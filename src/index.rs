use std::collections::{HashMap, HashSet};
use std::fmt;
use std::io::BufRead;

use crate::annotation::io::{AnnotationReader, AnnotationRecord, ParseError};
use crate::model::gene::Gene;
use crate::model::transcript::Transcript;
use crate::model::types::{GeneId, MatchClass, MatchHit, MatchOptions, TranscriptId};
#[allow(unused_imports)]
use crate::types::{RefBlock, SplicedRead, Strand};

// to serialize the data
use anyhow::{Context, Result, bail};
use flate2::read::GzDecoder;
use serde::{Deserialize, Serialize};
use std::fs::File;
use std::io::{BufReader, Read, Write};
use std::path::Path;

const MAGIC: &[u8; 4] = b"SPX1";
const VERSION_STR: &str = env!("CARGO_PKG_VERSION");

/// Configure which attribute keys are used to extract:
/// - gene stable identifier (used to intern -> GeneId)
/// - gene display names/aliases (stored in Gene.names)
/// - transcript stable identifier (used to intern -> TranscriptId)
/// - transcript display names/aliases (stored in Transcript.names)
/// - (GFF3) exon -> transcript linking keys (usually Parent)
///
/// Notes:
/// - We allow multiple keys per category; first present wins.
/// - For GFF3 Parent values, we split by ',' and treat each parent as a transcript ID.
#[derive(Debug, Clone)]
pub struct IdNameKeys {
    pub gene_id_keys: Vec<String>,
    pub gene_name_keys: Vec<String>,

    pub transcript_id_keys: Vec<String>,
    pub transcript_name_keys: Vec<String>,

    /// GFF3 exon->transcript linkage (most commonly: Parent)
    pub parent_keys: Vec<String>,

    /// Feature types that count as exon blocks (default: ["exon"])
    pub exon_feature_types: Vec<String>,
}

impl Default for IdNameKeys {
    fn default() -> Self {
        Self {
            // Common GTF + some common variants
            gene_id_keys: vec!["gene_id".into(), "gene".into(), "GeneID".into()],
            gene_name_keys: vec!["gene_name".into(), "Name".into(), "gene".into()],

            transcript_id_keys: vec!["transcript_id".into(), "transcript".into(), "ID".into()],
            transcript_name_keys: vec!["transcript_name".into(), "Name".into()],

            parent_keys: vec!["Parent".into()],
            exon_feature_types: vec!["exon".into()],
        }
    }
}

/// Per-chromosome bucket index: bin -> transcript ids.
///
/// This is a pre-filter only: it returns candidate transcript IDs that overlap bins.
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct ChrBuckets {
    pub bin_width: u32,
    pub bins: Vec<Vec<TranscriptId>>,
    pub max_end: u32,
}

impl ChrBuckets {
    pub fn new(bin_width: u32) -> Self {
        Self {
            bin_width,
            bins: Vec::new(),
            max_end: 0,
        }
    }

    fn ensure_len_for_end(&mut self, end0: u32) {
        self.max_end = self.max_end.max(end0);

        let need_bins =
            ((self.max_end as u64 + self.bin_width as u64 - 1) / self.bin_width as u64) as usize;
        if self.bins.len() < need_bins {
            self.bins.resize_with(need_bins, Vec::new);
        }
    }

    fn add_span(&mut self, tx_id: TranscriptId, start0: u32, end0: u32) {
        if end0 <= start0 {
            return;
        }

        self.ensure_len_for_end(end0);

        let b0 = (start0 / self.bin_width) as usize;
        let b1 = ((end0.saturating_sub(1)) / self.bin_width) as usize;

        for b in b0..=b1 {
            self.bins[b].push(tx_id);
        }
    }

    pub fn finalize_by_tx_start(&mut self, tx_span_start: &[u32], tx_span_end: &[u32]) {
        for bin in &mut self.bins {
            // Sort by (start, end, id) to make it deterministic.
            bin.sort_unstable_by(|a, b| {
                let sa = tx_span_start[*a];
                let sb = tx_span_start[*b];
                match sa.cmp(&sb) {
                    std::cmp::Ordering::Equal => {
                        let ea = tx_span_end[*a];
                        let eb = tx_span_end[*b];
                        match ea.cmp(&eb) {
                            std::cmp::Ordering::Equal => a.cmp(b),
                            other => other,
                        }
                    }
                    other => other,
                }
            });

            // Now adjacent duplicates are guaranteed adjacent (because equal ids compare equal)
            bin.dedup();
        }
    }
}

/// Best transcript matches (ties allowed).
///
/// This borrows the transcript from the index (zero-copy).
#[derive(Debug, Clone)]
pub struct TranscriptMatch<'a> {
    pub transcript_id: TranscriptId,
    pub transcript: &'a Transcript,
    pub hit: MatchHit,
}

/// Best gene matches (ties allowed).
///
/// `winning_transcripts` are the transcript(s) of this gene that achieved `best_hit`.
#[derive(Debug, Clone)]
pub struct GeneMatch<'a> {
    pub gene_id: GeneId,
    pub gene: &'a Gene,
    pub best_hit: MatchHit,
    pub winning_transcripts: Vec<&'a Transcript>,
}

/// The owning index type:
/// - chromosome dictionary (chr name -> chr_id)
/// - genes + transcripts
/// - per-chromosome buckets for fast candidate lookup
#[derive(Debug, Clone, Serialize, Deserialize)]
pub struct SpliceIndex {
    pub bin_width: u32,

    pub chr_names: Vec<String>,
    chr_to_id: HashMap<String, usize>,

    pub genes: Vec<Gene>,
    pub transcripts: Vec<Transcript>,
    transcript_by_name: HashMap<String, usize>,

    pub chr_buckets: Vec<ChrBuckets>,

    // Cached transcript spans, indexed by TranscriptId
    pub tx_span_start: Vec<u32>,
    pub tx_span_end: Vec<u32>,
}

/// Human-readable summary of the `SpliceIndex`.
///
/// The output is intended for quick inspection and debugging. It prints:
///
/// Global summary:
/// - Total number of genes
/// - Total number of transcripts
/// - Number of chromosomes indexed
/// - Global bin width in base pairs
///
/// Per-chromosome summary:
/// - Number of bins
/// - Number of unique genes represented on the chromosome
/// - Number of unique transcripts represented on the chromosome
/// - Mean number of genes per bin
/// - Mean number of transcripts per bin
///
/// Notes:
/// - Means are calculated over all bins, including empty ones.
/// - Gene counts per bin are derived from the transcripts stored in each bin.
/// - All statistics are computed on the fly during formatting and are not cached.
///
/// This implementation is meant for logging and diagnostics, not for
/// high-performance or machine-readable output.
impl fmt::Display for SpliceIndex {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        // ---- global summary ----
        let n_genes = self.genes.len();
        let n_txs = self.transcripts.len();
        let n_chrs = self.chr_names.len();

        writeln!(
            f,
            "SpliceIndex: {} genes, {} transcripts, {} chromosomes, bin_width={} bp",
            n_genes, n_txs, n_chrs, self.bin_width
        )?;

        if let Some(g) = self.genes.get(0) {
            writeln!(f, "Gene Names like: {}", g.names.join(", "))?;
        } else {
            writeln!(f, "Gene Names like: No genes detected")?;
        }

        // ---- per-chromosome stats ----
        for (i, chr_name) in self.chr_names.iter().enumerate() {
            let Some(chr) = self.chr_buckets.get(i) else {
                writeln!(f, "  - {}: <missing ChrBuckets>", chr_name)?;
                continue;
            };

            let nbins = chr.bins.len();

            let mut total_tx_hits: u64 = 0;
            let mut total_gene_hits: u64 = 0;

            let mut uniq_txs: HashSet<TranscriptId> = HashSet::new();
            let mut uniq_genes: HashSet<GeneId> = HashSet::new();

            for bin in chr.bins.iter() {
                total_tx_hits += bin.len() as u64;

                let mut bin_genes: HashSet<GeneId> = HashSet::new();

                for &tx_id in bin.iter() {
                    uniq_txs.insert(tx_id);

                    // assumes Transcript has a gene_id field
                    let gene_id = self.transcripts[tx_id].gene_id;
                    uniq_genes.insert(gene_id);
                    bin_genes.insert(gene_id);
                }

                total_gene_hits += bin_genes.len() as u64;
            }

            let mean_tx_per_bin = if nbins == 0 {
                0.0
            } else {
                total_tx_hits as f64 / nbins as f64
            };

            let mean_genes_per_bin = if nbins == 0 {
                0.0
            } else {
                total_gene_hits as f64 / nbins as f64
            };

            writeln!(
                f,
                "  - {}: bins={}, genes={}, transcripts={}, mean_genes/bin={:.3}, mean_tx/bin={:.3}",
                chr_name,
                nbins,
                uniq_genes.len(),
                uniq_txs.len(),
                mean_genes_per_bin,
                mean_tx_per_bin
            )?;
        }

        Ok(())
    }
}

impl SpliceIndex {
    pub fn new(bin_width: u32) -> Self {
        Self {
            bin_width,
            chr_names: Vec::new(),
            chr_to_id: HashMap::new(),
            genes: Vec::new(),
            transcripts: Vec::new(),
            transcript_by_name: HashMap::new(),
            chr_buckets: Vec::new(),
            tx_span_start: Vec::new(),
            tx_span_end: Vec::new(),
        }
    }

    /// sort the internal data structure by transcript start position.
    fn finalize(&mut self) {
        for cb in &mut self.chr_buckets {
            cb.finalize_by_tx_start(&self.tx_span_start, &self.tx_span_end);
        }
    }

    /// get all gene names from the index (id == gene_id)
    pub fn gene_names(&self) -> Vec<String> {
        self.genes
            .iter()
            .map(|g| g.primary_name().unwrap_or("NA").to_string())
            .collect()
    }

    /// get all transcript names from the index (id == transcript_id)
    pub fn transcript_names(&self) -> Vec<String> {
        self.transcripts
            .iter()
            .map(|t| t.primary_name().unwrap_or("NA").to_string())
            .collect()
    }

    /// Build a `SpliceIndex` from a GTF/GFF path.
    ///
    /// Chromosome order is derived from the annotation file itself:
    /// the first time we see a chromosome name in column 1, we append it to
    /// `chr_names` (stable "first-seen" order).
    ///
    /// This avoids requiring an external chromosome ordering for the common case.
    pub fn from_path<P: AsRef<Path>>(path: P, bin_width: u32, keys: IdNameKeys) -> Result<Self> {
        let path = path.as_ref();

        let is_gz = path
            .extension()
            .and_then(|s| s.to_str())
            .map(|ext| ext.eq_ignore_ascii_case("gz"))
            .unwrap_or(false);

        // Helper: open as BufRead (plain or gz)
        let open_bufread = || -> Result<Box<dyn BufRead>> {
            let f = File::open(path)
                .with_context(|| format!("open annotation file {}", path.display()))?;

            if is_gz {
                let gz = GzDecoder::new(f);
                Ok(Box::new(BufReader::new(gz)))
            } else {
                Ok(Box::new(BufReader::new(f)))
            }
        };

        // ----------------------------
        // Pass 1: discover chr order (first-seen in file)
        // ----------------------------
        let mut reader =
            open_bufread().with_context(|| format!("prepare reader for {}", path.display()))?;

        let mut chr_names: Vec<String> = Vec::new();
        let mut chr_to_id: HashMap<String, usize> = HashMap::new();

        let mut line = String::new();
        loop {
            line.clear();
            let n = reader
                .read_line(&mut line)
                .with_context(|| format!("read annotation {}", path.display()))?;
            if n == 0 {
                break;
            }

            let s = line.trim_end();
            if s.is_empty() || s.starts_with('#') {
                continue;
            }

            let chr = match s.split('\t').next() {
                Some(c) if !c.is_empty() => c,
                _ => continue,
            };

            if !chr_to_id.contains_key(chr) {
                let id = chr_names.len();
                chr_names.push(chr.to_string());
                chr_to_id.insert(chr.to_string(), id);
            }
        }
        // rebuild it fuzzy!
        let chr_to_id = Self::build_chr_map(&chr_names);

        // ----------------------------
        // Init index with discovered order
        // ----------------------------
        let mut idx = SpliceIndex {
            bin_width,

            chr_names,
            chr_to_id,

            genes: Vec::new(),
            transcripts: Vec::new(),
            transcript_by_name: HashMap::new(),

            chr_buckets: Vec::new(),

            tx_span_start: Vec::new(),
            tx_span_end: Vec::new(),
        };

        idx.chr_buckets = (0..idx.chr_names.len())
            .map(|_| ChrBuckets {
                bin_width,
                bins: Vec::new(),
                max_end: 0,
            })
            .collect();

        // ----------------------------
        // Pass 2: actual parse/build using your existing from_reader
        // ----------------------------
        let reader =
            open_bufread().with_context(|| format!("re-open reader for {}", path.display()))?;

        idx.from_reader(reader, keys)
            .with_context(|| format!("build splice index from {}", path.display()))
    }

    /// Build a fuzzy chromosome name map.
    ///
    /// Examples:
    /// - `chr1` maps also from `1`
    /// - `1` maps also from `chr1`
    /// - `MT`, `M`, `chrM` are treated as aliases
    pub fn build_chr_map(chr_names: &[String]) -> HashMap<String, usize> {
        let mut map = HashMap::with_capacity(chr_names.len() * 4);

        for (chr_id, name) in chr_names.iter().enumerate() {
            map.entry(name.clone()).or_insert(chr_id);

            if let Some(no_chr) = name.strip_prefix("chr") {
                map.entry(no_chr.to_string()).or_insert(chr_id);
            } else {
                map.entry(format!("chr{name}")).or_insert(chr_id);
            }

            match name.as_str() {
                "MT" => {
                    map.entry("chrM".to_string()).or_insert(chr_id);
                    map.entry("M".to_string()).or_insert(chr_id);
                }
                "chrM" => {
                    map.entry("MT".to_string()).or_insert(chr_id);
                    map.entry("M".to_string()).or_insert(chr_id);
                }
                "M" => {
                    map.entry("MT".to_string()).or_insert(chr_id);
                    map.entry("chrM".to_string()).or_insert(chr_id);
                }
                _ => {}
            }
        }

        map
    }

    pub fn gene_name(&self, gene_id: GeneId) -> Option<&str> {
        self.genes.get(gene_id).and_then(|g| g.primary_name())
    }

    pub fn transcript_by_name(&self, transcript_name: &str) -> Result<&Transcript, String> {
        let tx_id = self
            .transcript_by_name
            .get(transcript_name)
            .ok_or_else(|| format!("unknown transcript name {transcript_name:?}"))?;

        self.transcripts
            .get(*tx_id)
            .ok_or_else(|| format!("transcript id {tx_id} for {transcript_name:?} is out of bounds"))
    }



    pub fn transcript_name(&self, tx_id: TranscriptId) -> Option<&str> {
        self.transcripts.get(tx_id).and_then(|t| t.primary_name())
    }
    /// Build an index directly from a GTF/GFF3 reader.
    ///
    /// Workflow:
    /// 1) parse records
    /// 2) for exon features:
    ///    - extract gene key and transcript key(s)
    ///    - intern gene and transcript
    ///    - add exon block to transcript
    /// 3) finalize transcripts (sort/merge exons)
    /// 4) link transcripts into genes
    /// 5) build buckets
    ///
    /// # Example (minimal GTF via `BufRead`)
    /// ```
    /// use std::io::Cursor;
    /// use gtf_splice_index::index::{SpliceIndex, IdNameKeys};
    ///
    /// let gtf = "\
    /// chr1\tsrc\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n\
    /// chr1\tsrc\texon\t201\t250\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";\n";
    ///
    /// let idx = SpliceIndex::new(100)
    ///     .from_reader(Cursor::new(gtf.as_bytes()), IdNameKeys::default())
    ///     .unwrap();
    ///
    /// assert_eq!(idx.genes.len(), 1);
    /// assert_eq!(idx.transcripts.len(), 1);
    /// assert_eq!(idx.chr_names, vec!["chr1".to_string()]);
    /// ```
    pub fn from_reader<R: BufRead>(
        mut self,
        reader: R,
        keys: IdNameKeys,
    ) -> Result<Self, ParseError> {
        let mut transcript_by_name = HashMap::new();
        // stable key string -> internal id
        let mut gene_key_to_id: HashMap<String, GeneId> = HashMap::new();
        let mut tx_key_to_id: HashMap<String, TranscriptId> = HashMap::new();

        // transcript internal id -> gene internal id (best-effort)
        let mut tx_to_gene: HashMap<TranscriptId, GeneId> = HashMap::new();

        for rec in AnnotationReader::new(reader).records() {
            let rec = rec?;

            if !rec.is_exon_feature(&keys.exon_feature_types) {
                continue;
            }

            let chr_id = self.intern_chr(&rec.seqname);

            let gene_key = rec.pick_first_attr(&keys.gene_id_keys).ok_or_else(|| {
                ParseError::MalformedLine {
                    line_no: Some(rec.line_no),
                    expected: "gene_id attribute (configure via --gene-id-key)",
                    problem: "missing required gene id attribute".to_string(),
                    details: Some(format!("tried keys: {:?}", keys.gene_id_keys)),
                    line_preview: rec.line_preview.clone(),
                }
            })?;

            let tx_key_raw = rec
                .pick_first_attr(&keys.transcript_id_keys)
                .or_else(|| rec.pick_first_attr(&keys.parent_keys))
                .ok_or_else(|| {
                    let mut tried = keys.transcript_id_keys.clone();
                    tried.extend(keys.parent_keys.clone());
                    ParseError::MalformedLine {
                        line_no: Some(rec.line_no),
                        expected: "transcript_id attribute (or Parent for GFF3; configure via --transcript-id-key / --parent-key)",
                        problem: "missing required transcript id attribute".to_string(),
                        details: Some(format!("tried keys: {:?}", tried)),
                        line_preview: rec.line_preview.clone(),
                    }
                })?;

            // Parent can be comma-separated in GFF3; support multi-parent exons.
            let tx_keys: Vec<String> = split_gff3_parent_list(&tx_key_raw);

            // Intern gene (create once)
            let gene_id = self.intern_gene(&rec, &keys, &gene_key, &mut gene_key_to_id);

            // For each transcript key, intern transcript and add exon
            for tx_key in tx_keys {
                let tx_id =
                    self.intern_tx(&rec, &keys, chr_id, gene_id, &tx_key, &mut tx_key_to_id);

                tx_to_gene.entry(tx_id).or_insert(gene_id);

                self.transcripts[tx_id].add_exon(RefBlock {
                    start: rec.start0,
                    end: rec.end0,
                });
            }
        }

        // Finalize transcripts (sort/merge exons)

        self.tx_span_start.reserve(self.transcripts.len());
        self.tx_span_end.reserve(self.transcripts.len());

        for tx in &mut self.transcripts {
            let (start, end) = tx.finalize();
            self.tx_span_start.push(start);
            self.tx_span_end.push(end);
            transcript_by_name.insert( 
                tx.primary_name()
                .expect("Lib problems - transcript has no primary name!")
                .to_string(), tx.id 
            );
        }

        // Link transcripts into genes
        for (tx_id, gene_id) in tx_to_gene {
            self.genes[gene_id].add_transcript(tx_id);
        }
        for g in &mut self.genes {
            g.finalize();
        }

        // Build buckets
        self.build_buckets();
        self.transcript_by_name = transcript_by_name;
        Ok(self)
    }

    /// get the id for a chr name
    pub fn chr_id(&self, chr: &str) -> Option<usize> {
        self.chr_to_id.get(chr).copied()
    }

    /// Candidate transcripts for a span using bucket prefilter (union across bins).
    ///
    /// Returns a deduped Vec of TranscriptIds.
    pub fn candidates_for_span_union(
        &self,
        chr_id: usize,
        start0: u32,
        end0: u32,
    ) -> Vec<TranscriptId> {
        if chr_id >= self.chr_buckets.len() {
            return Vec::new();
        }

        let chr_bins = &self.chr_buckets[chr_id];

        if chr_bins.bins.is_empty() || end0 <= start0 {
            return Vec::new();
        }

        let first_bin = (start0 / chr_bins.bin_width) as usize;
        if first_bin >= chr_bins.bins.len() {
            return Vec::new();
        }

        let last_bin = ((end0.saturating_sub(1)) / chr_bins.bin_width) as usize;
        let last_bin = last_bin.min(chr_bins.bins.len() - 1);

        let mut candidates: HashSet<TranscriptId> = HashSet::new();

        for bin_idx in first_bin..=last_bin {
            let transcripts_in_bin = &chr_bins.bins[bin_idx];
            if transcripts_in_bin.is_empty() {
                continue;
            }

            // Binary search:
            // find first transcript where start >= end0
            let cutoff = Self::partition_point(
                transcripts_in_bin,
                |&tx_id| self.tx_span_start[tx_id] < end0,
            );

            // Only transcripts with start < end0 can overlap
            for &tx_id in &transcripts_in_bin[..cutoff] {
                // True overlap condition
                if self.tx_span_end[tx_id] > start0 {
                    candidates.insert(tx_id);
                }
            }
        }
        // Deduplicate across bins
        let ret = candidates.into_iter().collect();

        ret
    }

    #[inline]
    fn partition_point<T, F>(slice: &[T], mut pred: F) -> usize
    where
        F: FnMut(&T) -> bool,
    {
        let mut left = 0usize;
        let mut right = slice.len();

        while left < right {
            let mid = left + (right - left) / 2;

            if pred(&slice[mid]) {
                left = mid + 1;
            } else {
                right = mid;
            }
        }

        left
    }

    /// Convenience: candidates for a spliced read (based on its span).
    pub fn candidates_for_read_union(&self, read: &SplicedRead) -> Vec<TranscriptId> {
        let Some((s, e)) = read.span() else {
            return Vec::new();
        };
        self.candidates_for_span_union(read.chr_id, s, e)
    }

    /// Match a spliced read against transcripts:
    /// - buckets -> candidates
    /// - Transcript::match_spliced_read for each candidate
    ///
    /// Returns best transcript matches (ties allowed), deterministically ordered.
    pub fn match_transcripts<'a>(
        &'a self,
        read: &SplicedRead,
        opts: MatchOptions,
    ) -> Vec<TranscriptMatch<'a>> {
        read.assert_finalized();

        let candidates = self.candidates_for_read_union(read);
        if candidates.is_empty() {
            return Vec::new();
        }

        let mut scored: Vec<TranscriptMatch<'a>> = Vec::with_capacity(candidates.len());
        for tx_id in candidates {
            let tx = &self.transcripts[tx_id];
            let hit = tx.match_spliced_read(read, opts);
            scored.push(TranscriptMatch {
                transcript_id: tx_id,
                transcript: tx,
                hit,
            });
        }

        // Best class among candidates (requires MatchClass: Ord implemented in model/types.rs).
        let best_class = scored
            .iter()
            .map(|m| m.hit.class)
            .max()
            .unwrap_or(MatchClass::NoOverlap);

        // Keep ties by best class.
        scored.retain(|m| m.hit.class == best_class);

        // Deterministic: smaller total overhang first, then transcript_id.
        scored.sort_by_key(|m| (m.hit.overhang_5p_bp + m.hit.overhang_3p_bp, m.transcript_id));

        scored
    }

    /// Match a spliced read against genes by:
    /// - matching candidate transcripts
    /// - grouping by gene_id
    /// - picking the best transcript-hit per gene
    ///
    /// Returns best gene matches (ties allowed).
    pub fn match_genes<'a>(&'a self, read: &SplicedRead, opts: MatchOptions) -> Vec<GeneMatch<'a>> {
        read.assert_finalized();

        let candidates = self.candidates_for_read_union(read);
        if candidates.is_empty() {
            return Vec::new();
        }

        // Per gene: best class + best hit + winning tx ids
        let mut per_gene: HashMap<GeneId, (MatchClass, MatchHit, Vec<TranscriptId>)> =
            HashMap::new();

        for tx_id in candidates {
            let tx = &self.transcripts[tx_id];
            let hit = tx.match_spliced_read(read, opts);
            let gid: GeneId = tx.gene_id; // Transcript stores usize; treat as GeneId alias.

            match per_gene.get_mut(&gid) {
                None => {
                    per_gene.insert(gid, (hit.class, hit, vec![tx_id]));
                }
                Some((best_class, best_hit, winners)) => {
                    if hit.class > *best_class {
                        *best_class = hit.class;
                        *best_hit = hit;
                        winners.clear();
                        winners.push(tx_id);
                    } else if hit.class == *best_class {
                        winners.push(tx_id);

                        // Keep a deterministic representative best_hit:
                        // prefer smaller total overhang.
                        let cur_sum = best_hit.overhang_5p_bp + best_hit.overhang_3p_bp;
                        let new_sum = hit.overhang_5p_bp + hit.overhang_3p_bp;
                        if new_sum < cur_sum {
                            *best_hit = hit;
                        }
                    }
                }
            }
        }

        // Overall best gene class.
        let overall_best_class = per_gene
            .values()
            .map(|(c, _, _)| *c)
            .max()
            .unwrap_or(MatchClass::NoOverlap);

        // Collect tied best genes.
        let mut out: Vec<GeneMatch<'a>> = Vec::new();
        for (gid, (best_class, best_hit, mut txs)) in per_gene {
            if best_class != overall_best_class {
                continue;
            }

            txs.sort_unstable();
            txs.dedup();

            let gene = &self.genes[gid];
            let mut winning_transcripts: Vec<&'a Transcript> =
                txs.iter().map(|&tid| &self.transcripts[tid]).collect();

            // Deterministic order within gene: by tx_id
            winning_transcripts.sort_by_key(|t| t.id);

            out.push(GeneMatch {
                gene_id: gid,
                gene,
                best_hit,
                winning_transcripts,
            });
        }

        // Deterministic order across genes: smaller overhang sum first, then gene_id.
        out.sort_by_key(|g| {
            (
                g.best_hit.overhang_5p_bp + g.best_hit.overhang_3p_bp,
                g.gene_id,
            )
        });
        out
    }

    // -----------------------
    // Internal helpers
    // -----------------------

    fn intern_chr(&mut self, chr: &str) -> usize {
        if let Some(&id) = self.chr_to_id.get(chr) {
            return id;
        }
        let id = self.chr_names.len();
        self.chr_names.push(chr.to_string());
        self.chr_to_id.insert(chr.to_string(), id);
        self.chr_buckets.push(ChrBuckets::new(self.bin_width));
        id
    }

    fn intern_gene(
        &mut self,
        rec: &AnnotationRecord,
        keys: &IdNameKeys,
        gene_key: &str,
        gene_key_to_id: &mut HashMap<String, GeneId>,
    ) -> GeneId {
        if let Some(&gid) = gene_key_to_id.get(gene_key) {
            // Add any found name aliases
            for k in &keys.gene_name_keys {
                if let Some(v) = rec.attr(k) {
                    self.genes[gid].add_name(v);
                }
            }
            self.genes[gid].add_name(gene_key);
            return gid;
        }

        // Choose primary display name: first gene_name_keys present, else gene_key
        let primary = rec
            .pick_first_attr(&keys.gene_name_keys)
            .unwrap_or_else(|| gene_key.to_string());

        let gid = self.genes.len();
        self.genes.push(Gene::new(gid, primary));
        gene_key_to_id.insert(gene_key.to_string(), gid);

        // Store stable key as alias + any other display names
        self.genes[gid].add_name(gene_key);
        for k in &keys.gene_name_keys {
            if let Some(v) = rec.attr(k) {
                self.genes[gid].add_name(v);
            }
        }

        gid
    }

    fn intern_tx(
        &mut self,
        rec: &AnnotationRecord,
        keys: &IdNameKeys,
        chr_id: usize,
        gene_id: GeneId,
        tx_key: &str,
        tx_key_to_id: &mut HashMap<String, TranscriptId>,
    ) -> TranscriptId {
        if let Some(&tid) = tx_key_to_id.get(tx_key) {
            self.transcripts[tid].add_name(tx_key);
            for k in &keys.transcript_name_keys {
                if let Some(v) = rec.attr(k) {
                    self.transcripts[tid].add_name(v);
                }
            }
            return tid;
        }

        let primary = rec
            .pick_first_attr(&keys.transcript_name_keys)
            .unwrap_or_else(|| tx_key.to_string());

        let tid = self.transcripts.len();
        self.transcripts
            .push(Transcript::new(tid, gene_id, primary, chr_id, rec.strand));
        tx_key_to_id.insert(tx_key.to_string(), tid);

        self.transcripts[tid].add_name(tx_key);
        for k in &keys.transcript_name_keys {
            if let Some(v) = rec.attr(k) {
                self.transcripts[tid].add_name(v);
            }
        }

        tid
    }

    fn build_buckets(&mut self) {
        for (tx_id, tx) in self.transcripts.iter().enumerate() {
            let Some((s, e)) = tx.span() else {
                continue;
            };
            self.chr_buckets[tx.chr_id].add_span(tx_id, s, e);
        }
        self.finalize();
    }
    /// Serialize this index with a small header (magic + crate version) and a bincode payload.
    pub fn save(&self, path: impl AsRef<Path>) -> Result<()> {
        let path = path.as_ref();
        let mut f = File::create(path)?;

        // magic
        f.write_all(MAGIC)?;

        // version string from Cargo.toml
        let v = VERSION_STR.as_bytes();
        let len = v.len() as u16;
        f.write_all(&len.to_le_bytes())?;
        f.write_all(v)?;

        // payload
        let payload = bincode::serialize(self)?;
        f.write_all(&payload)?;

        Ok(())
    }

    /// Load an index written by `save()`. Rejects wrong file types and version mismatches.
    pub fn load(path: impl AsRef<Path>) -> Result<Self> {
        let path = path.as_ref();
        let mut f = File::open(path)?;

        // check magic
        let mut magic = [0u8; 4];
        f.read_exact(&mut magic)?;
        if &magic != MAGIC {
            bail!("Not a SpliceIndex file (bad magic)");
        }

        // read version string
        let mut len_buf = [0u8; 2];
        f.read_exact(&mut len_buf)?;
        let len = u16::from_le_bytes(len_buf) as usize;

        let mut ver_buf = vec![0u8; len];
        f.read_exact(&mut ver_buf)?;
        let file_version = std::str::from_utf8(&ver_buf)?;

        if file_version != VERSION_STR {
            bail!(
                "Index version mismatch: file={}, binary={}",
                file_version,
                VERSION_STR
            );
        }

        // read payload
        let mut payload = Vec::new();
        f.read_to_end(&mut payload)?;
        let idx: Self = bincode::deserialize(&payload)?;

        Ok(idx)
    }

    pub fn build_chr_map_fuzzy(&self) -> std::collections::HashMap<String, usize> {
        let mut map = std::collections::HashMap::new();

        for (i, n) in self.chr_names.iter().enumerate() {
            map.entry(n.clone()).or_insert(i);

            let no_chr = n.strip_prefix("chr").unwrap_or(n).to_string();
            map.entry(no_chr).or_insert(i);

            let with_chr = if n.starts_with("chr") {
                n.clone()
            } else {
                format!("chr{n}")
            };
            map.entry(with_chr).or_insert(i);

            if n == "MT" {
                map.entry("chrM".to_string()).or_insert(i);
                map.entry("M".to_string()).or_insert(i);
            }

            if n == "chrM" {
                map.entry("MT".to_string()).or_insert(i);
                map.entry("M".to_string()).or_insert(i);
            }
        }

        map
    }
}

/// Split Parent= list (GFF3) by commas; also trim whitespace.
fn split_gff3_parent_list(raw: &str) -> Vec<String> {
    let raw = raw.trim();
    if raw.contains(',') {
        raw.split(',')
            .map(|s| s.trim())
            .filter(|s| !s.is_empty())
            .map(|s| s.to_string())
            .collect()
    } else {
        vec![raw.to_string()]
    }
}

/*
#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn builds_index_from_minimal_gtf() {
        let gtf = "\
chr1\tsrc\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; gene_name \"Alpha\"; transcript_id \"T1\"; transcript_name \"TxA\";
chr1\tsrc\texon\t201\t250\t.\t+\t.\tgene_id \"G1\"; gene_name \"Alpha\"; transcript_id \"T1\"; transcript_name \"TxA\";
";
        let idx = SpliceIndex::new(100)
            .from_reader(Cursor::new(gtf.as_bytes()), IdNameKeys::default())
            .unwrap();

        assert_eq!(idx.chr_names, vec!["chr1".to_string()]);
        assert_eq!(idx.genes.len(), 1);
        assert_eq!(idx.transcripts.len(), 1);
        assert_eq!(idx.chr_buckets.len(), 1);
        assert!(idx.chr_buckets[0].bins.iter().any(|b| b.contains(&0)));
    }

    #[test]
    fn match_transcripts_picks_exact_chain() {
        let gtf = "\
chr1\tsrc\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\tsrc\texon\t201\t250\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\tsrc\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";
chr1\tsrc\texon\t201\t240\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T2\";
";
        let idx = SpliceIndex::new(50)
            .from_reader(Cursor::new(gtf.as_bytes()), IdNameKeys::default())
            .unwrap();

        // chr1 is first => chr_id 0
        let mut read = SplicedRead::new(
            0,
            Strand::Plus,
            vec![RefBlock::new(110, 150), RefBlock::new(200, 250)],
        );
        read.finalize();

        let best = idx.match_transcripts(&read, MatchOptions::default());
        assert!(!best.is_empty());
        assert_eq!(best[0].hit.class, MatchClass::ExactJunctionChain);
        assert_eq!(best.len(), 1);
        //panic!("Best vec: {:?}",best.iter().map(|b| b.hit.class ));
    }

    #[test]
    fn match_genes_returns_best_gene_and_winning_transcripts() {
        let gtf = "\
chr1\tsrc\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\tsrc\texon\t201\t250\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\tsrc\texon\t1001\t1050\t.\t+\t.\tgene_id \"G2\"; transcript_id \"T2\";
chr1\tsrc\texon\t1101\t1150\t.\t+\t.\tgene_id \"G2\"; transcript_id \"T2\";
";
        let idx = SpliceIndex::new(200)
            .from_reader(Cursor::new(gtf.as_bytes()), IdNameKeys::default())
            .unwrap();

        let mut read = SplicedRead::new(
            0,
            Strand::Plus,
            vec![RefBlock::new(110, 150), RefBlock::new(200, 250)],
        );
        read.finalize();

        let genes = idx.match_genes(&read, MatchOptions::default());
        assert_eq!(genes.len(), 1);
        assert_eq!(genes[0].best_hit.class, MatchClass::ExactJunctionChain);
        assert_eq!(genes[0].winning_transcripts.len(), 1);
    }

    #[test]
    fn gff3_parent_multi_value_creates_two_transcripts() {
        let gff = "\
chr2\tsrc\texon\t5\t20\t.\t-\t.\tParent=tx1,tx2;gene_id=G9;Name=GeneNice
";
        let mut keys = IdNameKeys::default();
        keys.transcript_id_keys = vec![]; // force Parent usage
        keys.parent_keys = vec!["Parent".into()];
        keys.gene_id_keys = vec!["gene_id".into()];
        keys.gene_name_keys = vec!["Name".into()];

        let idx = SpliceIndex::new(50)
            .from_reader(Cursor::new(gff.as_bytes()), keys)
            .unwrap();

        assert_eq!(idx.genes.len(), 1);
        assert_eq!(idx.transcripts.len(), 2);
        assert_eq!(idx.transcripts[0].strand, Strand::Minus);
        assert_eq!(idx.transcripts[1].strand, Strand::Minus);
        assert_eq!(idx.transcripts[0].exons().len(), 1);
        assert_eq!(idx.transcripts[0].exons()[0].start, 4);
        assert_eq!(idx.transcripts[0].exons()[0].end, 20);
    }
}
*/

#[cfg(test)]
mod tests {
    use super::*;

    use crate::model::types::{GeneId, TranscriptId};
    use std::collections::{HashMap, HashSet};
    use std::path::PathBuf;

    // If you already have a crate-local Result/Error type, keep using it.
    // Otherwise you can switch these tests to `anyhow::Result<()>` easily.
    type TestResult<T> = std::result::Result<T, Box<dyn std::error::Error>>;

    // -----------------------------
    // Helpers: stable temp path
    // -----------------------------
    fn tmp_path(name: &str) -> PathBuf {
        let mut p = std::env::temp_dir();
        // add a bit of uniqueness without external crates
        let pid = std::process::id();
        let t = std::time::SystemTime::now()
            .duration_since(std::time::UNIX_EPOCH)
            .unwrap()
            .as_nanos();
        p.push(format!("gtf_splice_index_test_{name}_{pid}_{t}.bin"));
        p
    }

    // -----------------------------
    // Helpers: ID constructors
    // -----------------------------
    // Adjust these if your GeneId/TranscriptId are not tuple structs.
    fn gid(i: usize) -> GeneId {
        i
    }
    fn tid(i: usize) -> TranscriptId {
        i
    }

    // -----------------------------
    // Helpers: create transcripts/genes (using your constructors)
    // -----------------------------
    fn make_gene(
        id: GeneId,
        primary: &str,
        aliases: &[&str],
        transcript_ids: &[TranscriptId],
    ) -> Gene {
        let mut g = Gene::new(id, primary);
        for &a in aliases {
            g.add_name(a);
        }
        for &tx in transcript_ids {
            g.add_transcript(tx);
        }
        g.finalize();
        g
    }

    fn make_transcript(
        id: TranscriptId,
        gene_id: usize,
        primary_name: &str,
        alias_names: &[&str],
        chr_id: usize,
        strand: Strand,
        exons: Vec<RefBlock>,
    ) -> Transcript {
        let mut tx = Transcript::new(id, gene_id, primary_name, chr_id, strand);
        for &a in alias_names {
            tx.add_name(a);
        }
        for b in exons {
            tx.add_exon(b);
        }
        tx.finalize(); // sorts + merges exons, sets finalized=true, returns span
        tx
    }

    fn span_of_exons(exons: &[RefBlock]) -> (u32, u32) {
        let mut s = u32::MAX;
        let mut e = 0u32;
        for b in exons {
            s = s.min(b.start);
            e = e.max(b.end);
        }
        (s, e)
    }

    // Fill buckets deterministically.
    // This assumes half-open coordinates [start,end) and bin assignment by overlap.
    fn rebuild_chr_buckets(index: &mut SpliceIndex) {
        index.chr_buckets.clear();

        // determine per-chromosome max_end
        let mut chr_max_end: Vec<u32> = vec![0; index.chr_names.len()];
        for tx in &index.transcripts {
            let (s, e) = span_of_exons(&tx.exons());
            let _ = s;
            chr_max_end[tx.chr_id] = chr_max_end[tx.chr_id].max(e);
        }

        // allocate bins and fill
        for chr_id in 0..index.chr_names.len() {
            let max_end = chr_max_end[chr_id];
            let nbins = if max_end == 0 {
                0
            } else {
                // ceiling division
                ((max_end as u64 + index.bin_width as u64 - 1) / index.bin_width as u64) as usize
            };
            let mut cb = ChrBuckets {
                bin_width: index.bin_width,
                bins: vec![Vec::<TranscriptId>::new(); nbins],
                max_end,
            };

            for tx in &index.transcripts {
                if tx.chr_id != chr_id {
                    continue;
                }
                let (start, end) = span_of_exons(&tx.exons());
                if end <= start {
                    continue;
                }
                let bw = index.bin_width;
                let b0 = (start / bw) as usize;
                let b1 = ((end - 1) / bw) as usize; // end-1 for half-open
                for b in b0..=b1 {
                    cb.bins[b].push(tx.id);
                }
            }

            // enforce deterministic contents: sort & dedup each bin

            index.chr_buckets.push(cb);
        }
        index.finalize();
    }

    fn rebuild_tx_spans(index: &mut SpliceIndex) {
        index.tx_span_start.clear();
        index.tx_span_end.clear();

        // We assume TranscriptId values are dense 0..N-1.
        let n = index.transcripts.len();
        index.tx_span_start.resize(n, 0);
        index.tx_span_end.resize(n, 0);

        for tx in &index.transcripts {
            let (s, e) = span_of_exons(&tx.exons());
            index.tx_span_start[tx.id] = s;
            index.tx_span_end[tx.id] = e;
        }
    }

    // Build a small but non-trivial index:
    //
    // bin_width = 100
    //
    // chr1:
    //   tx0: [ 10,  60) U [140, 180) => spans 10..180 => bins 0 and 1
    //   tx1: [200, 260)             => spans 200..260 => bins 2
    // chr2:
    //   tx2: [ 90, 110)             => spans 90..110  => bins 0 and 1 (crosses boundary)
    //
    fn build_index() -> SpliceIndex {
        let mut idx = SpliceIndex::new(100);

        // Because this is a unit test in the same module, we can fill private fields.
        idx.chr_names = vec!["chr1".to_string(), "chr2".to_string()];
        idx.chr_to_id = HashMap::from([("chr1".to_string(), 0usize), ("chr2".to_string(), 1usize)]);

        // transcripts
        let tx0 = make_transcript(
            tid(0),
            0, // gene index in idx.genes (below)
            "TX0",
            &["TX0", "TxZero"],
            0,
            Strand::Plus,
            vec![
                RefBlock { start: 10, end: 60 },
                RefBlock {
                    start: 140,
                    end: 180,
                },
            ],
        );

        let tx1 = make_transcript(
            tid(1),
            0,
            "TX1",
            &["TX1"],
            0,
            Strand::Plus,
            vec![RefBlock {
                start: 200,
                end: 260,
            }],
        );

        let tx2 = make_transcript(
            tid(2),
            1, // second gene
            "TX1",
            &["TX2"],
            1,
            Strand::Minus,
            vec![RefBlock {
                start: 90,
                end: 110,
            }],
        );

        // add this tx3
        let tx3 = make_transcript(
            tid(3),
            0, // same gene as tx0/tx1, or set 1 if you want a second gene on chr1
            "TX3",
            &["TX3"],
            0, // chr1
            Strand::Plus,
            vec![
                RefBlock {
                    start: 150,
                    end: 160,
                },
                RefBlock {
                    start: 210,
                    end: 240,
                },
            ],
        );

        idx.transcripts = vec![tx0, tx1, tx2, tx3];

        // genes
        let g0 = make_gene(gid(0), "G0", &["G0", "GeneZero"], &[tid(0), tid(1)]);
        let g1 = make_gene(gid(1), "G1", &["G1"], &[tid(2)]);
        idx.genes = vec![g0, g1];

        // build cached structures
        rebuild_tx_spans(&mut idx);
        rebuild_chr_buckets(&mut idx);

        idx
    }

    // -----------------------------
    // The "rigid structure" test
    // -----------------------------
    #[test]
    fn splice_index_internal_structure_is_consistent_and_deterministic() -> TestResult<()> {
        let idx = build_index();

        // ---- basic shape checks
        assert_eq!(idx.bin_width, 100);
        assert_eq!(idx.chr_names, vec!["chr1".to_string(), "chr2".to_string()]);
        assert_eq!(idx.genes.len(), 2);
        assert_eq!(idx.transcripts.len(), 4);
        assert_eq!(idx.chr_buckets.len(), 2);

        // ---- private map must match chr_names
        assert_eq!(idx.chr_to_id.get("chr1").copied(), Some(0));
        assert_eq!(idx.chr_to_id.get("chr2").copied(), Some(1));
        assert_eq!(idx.chr_to_id.len(), 2);

        // ---- tx span caches must be present and correct
        assert_eq!(idx.tx_span_start.len(), 4);
        assert_eq!(idx.tx_span_end.len(), 4);

        assert_eq!(idx.tx_span_start[0], 10);
        assert_eq!(idx.tx_span_end[0], 180);

        assert_eq!(idx.tx_span_start[1], 200);
        assert_eq!(idx.tx_span_end[1], 260);

        assert_eq!(idx.tx_span_start[2], 90);
        assert_eq!(idx.tx_span_end[2], 110);

        // ---- per-chromosome buckets must be consistent
        // chr1 max_end = 260 => nbins=3 (0..2)
        let cb1 = &idx.chr_buckets[0];
        assert_eq!(cb1.bin_width, 100);
        assert_eq!(cb1.max_end, 260);
        assert_eq!(cb1.bins.len(), 3);

        // chr1 bin membership (sorted, deduped)
        // bin0: tx0 (10..60)
        assert_eq!(cb1.bins[0], vec![tid(0)]);
        // bin1: tx0 (140..180)
        assert_eq!(cb1.bins[1], vec![tid(0), tid(3)]);
        // bin2: tx1 (200..260)
        assert_eq!(cb1.bins[2], vec![tid(3), tid(1)]);

        // chr2 max_end = 110 => nbins=2 (0..1)
        let cb2 = &idx.chr_buckets[1];
        assert_eq!(cb2.bin_width, 100);
        assert_eq!(cb2.max_end, 110);
        assert_eq!(cb2.bins.len(), 2);

        // chr2 bin membership: tx2 spans 90..110 crosses 100 boundary => bins 0 and 1
        assert_eq!(cb2.bins[0], vec![tid(2)]);
        assert_eq!(cb2.bins[1], vec![tid(2)]);

        // ---- candidates_for_span_union must union bins correctly and dedupe
        // chr1 span 0..150 touches bins 0 and 1 => tx0 only
        let mut c = idx.candidates_for_span_union(0, 0, 150);
        c.sort_by_key(|x| *x);
        assert_eq!(c, vec![tid(0)]);

        // chr1 span 150..250 touches bins 1 and 2 => tx0, tx1
        let mut c = idx.candidates_for_span_union(0, 150, 250);
        c.sort_by_key(|x| *x);
        assert_eq!(c, vec![tid(0), tid(1), tid(3)]);

        // chr2 span 95..96 touches only bin0 => tx2
        let mut c = idx.candidates_for_span_union(1, 95, 96);
        c.sort_by_key(|x| *x);
        assert_eq!(c, vec![tid(2)]);

        // ---- candidates_for_read_union should match span behavior across blocks
        // Adjust this section if your SplicedRead struct differs.
        let mut read_chr1 = SplicedRead::new(
            0,
            Strand::Plus,
            vec![
                RefBlock { start: 12, end: 20 },
                RefBlock {
                    start: 160,
                    end: 170,
                },
            ],
        );
        read_chr1.finalize();
        let mut c = idx.candidates_for_read_union(&read_chr1);
        c.sort_by_key(|x| *x);
        assert_eq!(c, vec![tid(0), tid(3)]);

        let mut read_chr1_wide = SplicedRead::new(
            0,
            Strand::Plus,
            vec![
                RefBlock {
                    start: 155,
                    end: 165,
                },
                RefBlock {
                    start: 205,
                    end: 210,
                },
            ],
        );
        read_chr1_wide.finalize();

        let mut c = idx.candidates_for_read_union(&read_chr1_wide);
        c.sort_by_key(|x| *x);
        assert_eq!(c, vec![tid(0), tid(1), tid(3)]);

        // ---- match_transcripts should pick the true best transcript
        // This assumes `MatchOptions` exists and your scoring prefers exon overlap and
        // correct splice structure. Adjust opts as needed.
        let opts = MatchOptions::default();

        let hits = idx.match_transcripts(&read_chr1, opts);
        assert!(!hits.is_empty(), "expected at least one transcript hit");
        // strongest expectation: tx0 must be a winner
        assert!(
            hits.iter().any(|h| h.transcript_id == tid(0)),
            "expected tx0 among winning transcript hits: {hits:?}"
        );
        // and no chr-mismatched tx2
        assert!(
            hits.iter().all(|h| h.transcript.chr_id == 0),
            "expected all hits to be on chr1"
        );

        // ---- match_genes should map to gene0 for chr1 reads
        let gene_hits = idx.match_genes(&read_chr1, opts);
        assert!(!gene_hits.is_empty(), "expected at least one gene hit");
        assert!(
            gene_hits.iter().any(|g| g.gene_id == gid(0)),
            "expected gene0 among winning gene hits: {gene_hits:?}"
        );
        assert!(
            gene_hits.iter().all(|g| g.gene.id == g.gene_id),
            "gene reference and gene_id must agree"
        );

        // winning transcripts for gene0 must be subset of {tx0, tx1}
        for g in &gene_hits {
            if g.gene_id == gid(0) {
                let s: HashSet<usize> = g.winning_transcripts.iter().map(|t| t.id).collect();
                assert!(
                    s.is_subset(&HashSet::from([0usize, 1usize])),
                    "unexpected winning transcript IDs for gene0: {s:?}"
                );
            }
        }

        // ---- Display should at least mention key global counts
        let s = format!("{idx}");
        // keep these weak enough to survive formatting changes, but still meaningful
        assert!(s.contains("genes") || s.contains("Genes"));
        assert!(s.contains("transcripts") || s.contains("Transcripts"));
        assert!(s.contains("chr") || s.contains("Chr") || s.contains("chrom"));

        Ok(())
    }

    // -----------------------------
    // Roundtrip test: save/load
    // -----------------------------
    #[test]
    fn splice_index_save_load_roundtrip_preserves_internal_structure() -> TestResult<()> {
        use std::io::Cursor;

        // -----------------------------
        // Inline GTF (2 genes, 3 transcripts)
        // -----------------------------
        let gtf = concat!(
            "chr1\tsource\texon\t10\t60\t.\t+\t.\tgene_id \"G0\"; transcript_id \"TX0\"; gene_name \"GeneZero\";\n",
            "chr1\tsource\texon\t140\t180\t.\t+\t.\tgene_id \"G0\"; transcript_id \"TX0\";\n",
            "chr1\tsource\texon\t200\t260\t.\t+\t.\tgene_id \"G0\"; transcript_id \"TX1\";\n",
            "chr2\tsource\texon\t90\t110\t.\t-\t.\tgene_id \"G1\"; transcript_id \"TX2\";\n",
        );

        // -----------------------------
        // Build index from text
        // -----------------------------
        let mut idx = SpliceIndex::new(50).from_reader(
            Cursor::new(gtf),
            IdNameKeys::default(), // your default keys
        )?;

        // IMPORTANT: finalize after build
        idx.finalize();

        let path = tmp_path("splice_index_inline_gtf");

        idx.save(&path)?;

        let loaded = SpliceIndex::load(&path)?;

        // ---- compare high-level invariants
        assert_eq!(loaded.bin_width, idx.bin_width);
        assert_eq!(loaded.chr_names, idx.chr_names);
        assert_eq!(loaded.genes.len(), idx.genes.len());
        assert_eq!(loaded.transcripts.len(), idx.transcripts.len());

        // ---- compare private dictionary too (unit-test: same module)
        assert_eq!(loaded.chr_to_id.len(), idx.chr_to_id.len());
        for (k, v) in &idx.chr_to_id {
            assert_eq!(loaded.chr_to_id.get(k).copied(), Some(*v));
        }

        // ---- compare cached spans
        assert_eq!(loaded.tx_span_start, idx.tx_span_start);
        assert_eq!(loaded.tx_span_end, idx.tx_span_end);

        // ---- compare buckets rigidly
        assert_eq!(loaded.chr_buckets.len(), idx.chr_buckets.len());
        for (a, b) in loaded.chr_buckets.iter().zip(idx.chr_buckets.iter()) {
            assert_eq!(a.bin_width, b.bin_width);
            assert_eq!(a.max_end, b.max_end);
            assert_eq!(a.bins.len(), b.bins.len());
            assert_eq!(a.bins, b.bins);
        }

        Ok(())
    }

    #[test]
    fn competitive_match_transcripts_and_genes_choose_correct_winner() -> TestResult<()> {
        use std::collections::HashMap;

        // -----------------------------
        // Build a competitive index:
        // chr1 has TX0 (G0) and TX3 (G2) overlapping the same bin(s)
        // chr2 has TX2 (G1) as a distractor on another chromosome
        // -----------------------------
        let mut idx = SpliceIndex::new(100);

        idx.chr_names = vec!["chr1".to_string(), "chr2".to_string()];
        idx.chr_to_id = HashMap::from([("chr1".to_string(), 0usize), ("chr2".to_string(), 1usize)]);

        // Gene layout:
        // genes[0] => G0 on chr1 with TX0 + TX1
        // genes[1] => G1 on chr2 with TX2
        // genes[2] => G2 on chr1 with TX3   (this is the competitor gene)
        let g0 = make_gene(gid(0), "G0", &["GeneZero"], &[]);
        let g1 = make_gene(gid(1), "G1", &[], &[]);
        let g2 = make_gene(gid(2), "G2", &[], &[]);
        idx.genes = vec![g0, g1, g2];

        // TX0 (gene index 0), chr1, two exons => junction (60,140)
        let tx0 = make_transcript(
            tid(0),
            0,
            "TX0",
            &["TxZero"],
            0,
            Strand::Plus,
            vec![RefBlock::new(10, 60), RefBlock::new(140, 180)],
        );

        // TX1 (gene index 0), chr1, single exon in bin2 (not needed for this test, but realistic)
        let tx1 = make_transcript(
            tid(1),
            0,
            "TX1",
            &[],
            0,
            Strand::Plus,
            vec![RefBlock::new(200, 260)],
        );

        // TX2 (gene index 1), chr2, crosses boundary (distractor)
        let tx2 = make_transcript(
            tid(2),
            1,
            "TX2",
            &[],
            1,
            Strand::Minus,
            vec![RefBlock::new(90, 110)],
        );

        // TX3 (gene index 2), chr1 competitor:
        // exons overlap read's 2nd block, but junction differs: (160,210)
        // span 150..240 => touches bin1 and bin2, so it will be a candidate for reads in bin1.
        let tx3 = make_transcript(
            tid(3),
            2,
            "TX3",
            &[],
            0,
            Strand::Plus,
            vec![RefBlock::new(150, 160), RefBlock::new(210, 240)],
        );

        idx.transcripts = vec![tx0, tx1, tx2, tx3];

        // Wire genes -> transcript_ids (important for match_genes)
        idx.genes[0].add_transcript(tid(0));
        idx.genes[0].add_transcript(tid(1));
        idx.genes[1].add_transcript(tid(2));
        idx.genes[2].add_transcript(tid(3));
        for g in &mut idx.genes {
            g.finalize();
        }

        // Build cached structures + buckets
        rebuild_tx_spans(&mut idx);
        rebuild_chr_buckets(&mut idx);

        // If your SpliceIndex has its own finalize() that sorts/dedups bins by tx start, call it here:
        // idx.finalize();

        // Sanity: prove we truly have a bin with 2 entries on chr1
        let cb1 = &idx.chr_buckets[0];
        assert_eq!(cb1.bins.len(), 3);
        // bin1 should contain TX0 (span 10..180) and TX3 (span 150..240)
        assert!(
            cb1.bins[1].contains(&tid(0)) && cb1.bins[1].contains(&tid(3)),
            "expected chr1 bin1 to contain TX0 and TX3, got: {:?}",
            cb1.bins[1]
        );

        // -----------------------------
        // Competitive read: exact for TX0, not exact for TX3
        // Read blocks match TX0 exons with the same junction (60,140)
        // -----------------------------
        let mut read_tx0 = SplicedRead::new(
            0,
            Strand::Plus,
            vec![RefBlock::new(12, 60), RefBlock::new(140, 170)],
        );
        read_tx0.finalize(); // REQUIRED for matching

        // Require exact junction chain and strand
        let mut opts = MatchOptions::default();
        opts.require_strand = true;
        opts.require_exact_junction_chain = true;
        opts.max_5p_overhang_bp = 50;
        opts.max_3p_overhang_bp = 50;

        // -----------------------------
        // 1) match_transcripts: TX0 must win (TX3 must not win under exact requirement)
        // -----------------------------
        let tx_hits = idx.match_transcripts(&read_tx0, opts);
        assert!(!tx_hits.is_empty(), "expected transcript hits");

        // Every returned hit should be ExactJunctionChain under this option set
        assert!(
            tx_hits
                .iter()
                .all(|h| h.hit.class == MatchClass::ExactJunctionChain),
            "expected only ExactJunctionChain hits, got: {tx_hits:?}"
        );

        // TX0 must be among winners
        assert!(
            tx_hits.iter().any(|h| h.transcript_id == tid(0)),
            "expected TX0 among winners, got: {tx_hits:?}"
        );

        // TX3 must NOT be a winner when exact chain is required (its junction differs)
        assert!(
            !tx_hits.iter().any(|h| h.transcript_id == tid(3)),
            "did not expect TX3 among winners under exact requirement, got: {tx_hits:?}"
        );

        // -----------------------------
        // 2) match_genes: G0 must win (and winning transcript must be TX0)
        // -----------------------------
        let gene_hits = idx.match_genes(&read_tx0, opts);
        assert!(!gene_hits.is_empty(), "expected gene hits");

        // Under exact requirement, best_hit should also be ExactJunctionChain
        assert!(
            gene_hits
                .iter()
                .all(|g| g.best_hit.class == MatchClass::ExactJunctionChain),
            "expected ExactJunctionChain gene best_hit, got: {gene_hits:?}"
        );

        // G0 (id 0) must be among winners
        assert!(
            gene_hits.iter().any(|g| g.gene_id == gid(0)),
            "expected G0 among gene winners, got: {gene_hits:?}"
        );

        // G2 (id 2) must NOT be among winners (because TX3 should not win)
        assert!(
            !gene_hits.iter().any(|g| g.gene_id == gid(2)),
            "did not expect G2 among winners under exact requirement, got: {gene_hits:?}"
        );

        // Winning transcripts for G0 must include TX0 (and should not include TX1 for this read)
        for g in &gene_hits {
            if g.gene_id == gid(0) {
                let ids: Vec<TranscriptId> = g.winning_transcripts.iter().map(|t| t.id).collect();
                assert!(
                    ids.contains(&tid(0)),
                    "expected TX0 to be a winning transcript for G0, got: {ids:?}"
                );
                assert!(
                    !ids.contains(&tid(1)),
                    "did not expect TX1 to be a winning transcript for this read, got: {ids:?}"
                );
            }
        }

        Ok(())
    }
}

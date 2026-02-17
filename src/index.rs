use std::collections::HashMap;
use std::io::BufRead;

use crate::annotation::io::{AnnotationReader, AnnotationRecord, ParseError};
use crate::model::gene::Gene;
use crate::model::transcript::Transcript;
use crate::model::types::{GeneId, MatchClass, MatchHit, MatchOptions, TranscriptId};
use crate::types::{RefBlock, SplicedRead, Strand};

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
#[derive(Debug, Clone)]
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

    fn finalize(&mut self) {
        for bin in &mut self.bins {
            bin.sort_unstable();
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
#[derive(Debug, Clone)]
pub struct SpliceIndex {
    pub bin_width: u32,

    pub chr_names: Vec<String>,
    chr_to_id: HashMap<String, usize>,

    pub genes: Vec<Gene>,
    pub transcripts: Vec<Transcript>,

    pub chr_buckets: Vec<ChrBuckets>,
}

impl SpliceIndex {
    pub fn new(bin_width: u32) -> Self {
        Self {
            bin_width,
            chr_names: Vec::new(),
            chr_to_id: HashMap::new(),
            genes: Vec::new(),
            transcripts: Vec::new(),
            chr_buckets: Vec::new(),
        }
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
        // stable key string -> internal id
        let mut gene_key_to_id: HashMap<String, GeneId> = HashMap::new();
        let mut tx_key_to_id: HashMap<String, TranscriptId> = HashMap::new();

        // transcript internal id -> gene internal id (best-effort)
        let mut tx_to_gene: HashMap<TranscriptId, GeneId> = HashMap::new();

        for rec in AnnotationReader::new(reader).records() {
            let rec = rec?;

            if !rec.is_exon_feature( &keys.exon_feature_types) {
                continue;
            }

            let chr_id = self.intern_chr(&rec.seqname);

            // Gene stable key (for interning)
            let gene_key = &rec.pick_first_attr( &keys.gene_id_keys)
                .unwrap_or_else(|| "<NO_GENE_ID>".to_string());

            // Transcript stable key(s) (for interning):
            // prefer transcript_id_keys, else use Parent keys (GFF3)
            let tx_key_raw = &rec.pick_first_attr( &keys.transcript_id_keys)
                .or_else(|| rec.pick_first_attr( &keys.parent_keys))
                .unwrap_or_else(|| "<NO_TX_ID>".to_string());

            // Parent can be comma-separated in GFF3; support multi-parent exons.
            let tx_keys: Vec<String> = split_gff3_parent_list(&tx_key_raw);

            // Intern gene (create once)
            let gene_id = self.intern_gene(&rec, &keys, &gene_key, &mut gene_key_to_id);

            // For each transcript key, intern transcript and add exon
            for tx_key in tx_keys {
                let tx_id = self.intern_tx(
                    &rec,
                    &keys,
                    chr_id,
                    gene_id,
                    &tx_key,
                    &mut tx_key_to_id,
                );

                tx_to_gene.entry(tx_id).or_insert(gene_id);

                self.transcripts[tx_id].add_exon(RefBlock {
                    start: rec.start0,
                    end: rec.end0,
                });
            }
        }

        // Finalize transcripts (sort/merge exons)
        for tx in &mut self.transcripts {
            tx.finalize();
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

        Ok(self)
    }

    /// Candidate transcripts for a span using bucket prefilter (union across bins).
    ///
    /// Returns a deduped Vec of TranscriptIds.
    pub fn candidates_for_span_union(&self, chr_id: usize, start0: u32, end0: u32) -> Vec<TranscriptId> {
        if chr_id >= self.chr_buckets.len() {
            return Vec::new();
        }
        let cb = &self.chr_buckets[chr_id];
        if cb.bins.is_empty() || end0 <= start0 {
            return Vec::new();
        }

        let b0 = (start0 / cb.bin_width) as usize;
        let b1 = ((end0.saturating_sub(1)) / cb.bin_width) as usize;

        if b0 >= cb.bins.len() {
            return Vec::new();
        }
        let b1 = b1.min(cb.bins.len() - 1);

        let mut out: Vec<TranscriptId> = Vec::new();
        for b in b0..=b1 {
            out.extend_from_slice(&cb.bins[b]);
        }
        out.sort_unstable();
        out.dedup();
        out
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
        let mut per_gene: HashMap<GeneId, (MatchClass, MatchHit, Vec<TranscriptId>)> = HashMap::new();

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
        out.sort_by_key(|g| (g.best_hit.overhang_5p_bp + g.best_hit.overhang_3p_bp, g.gene_id));
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
        let primary = rec.pick_first_attr(&keys.gene_name_keys).unwrap_or_else(|| gene_key.to_string());

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

        let primary = rec.pick_first_attr( &keys.transcript_name_keys).unwrap_or_else(|| tx_key.to_string());

        let tid = self.transcripts.len();
        self.transcripts.push(Transcript::new(
            tid,
            gene_id,
            primary,
            chr_id,
            rec.strand,
        ));
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
            let Some((s, e)) = tx.span() else { continue; };
            self.chr_buckets[tx.chr_id].add_span(tx_id, s, e);
        }
        for cb in &mut self.chr_buckets {
            cb.finalize();
        }
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

use std::io::{ BufRead, BufReader};
use std::path::Path;

use crate::index::{IdNameKeys, SpliceIndex};
use crate::annotation::io::ParseError;

/// High-level builder for creating a `SpliceIndex` from a GTF/GFF3 file.
///
/// This is the “annotation class” you asked for:
/// - parses whole file (optionally gzipped)
/// - configurable mapping of ID/NAME keys for gene + transcript
/// - builds genes + transcripts + chromosome dictionary + buckets
#[derive(Debug, Clone)]
pub struct AnnotationBuilder {
    pub bin_width: u32,
    pub keys: IdNameKeys,
}

impl AnnotationBuilder {
    /// Start with defaults that work reasonably for many GTF/GFF3 files.
    pub fn new(bin_width: u32) -> Self {
        Self {
            bin_width,
            keys: IdNameKeys::default(),
        }
    }

    /// Convenience: set a single key (or first-preference key) for gene id.
    pub fn gene_id_key(mut self, key: &str) -> Self {
        self.keys.gene_id_keys = vec![key.to_string()];
        self
    }

    /// Convenience: set preferred keys for gene display names/aliases.
    pub fn gene_name_keys(mut self, keys: &[&str]) -> Self {
        self.keys.gene_name_keys = keys.iter().map(|s| s.to_string()).collect();
        self
    }

    /// Convenience: set transcript id key(s).
    pub fn transcript_id_keys(mut self, keys: &[&str]) -> Self {
        self.keys.transcript_id_keys = keys.iter().map(|s| s.to_string()).collect();
        self
    }

    /// Convenience: set transcript name key(s).
    pub fn transcript_name_keys(mut self, keys: &[&str]) -> Self {
        self.keys.transcript_name_keys = keys.iter().map(|s| s.to_string()).collect();
        self
    }

    /// Convenience: set parent keys for GFF3 exon->transcript linking (usually ["Parent"]).
    pub fn parent_keys(mut self, keys: &[&str]) -> Self {
        self.keys.parent_keys = keys.iter().map(|s| s.to_string()).collect();
        self
    }

    /// Convenience: define what feature types count as exon blocks.
    /// Typical:
    /// - GTF: ["exon"]
    /// - GFF3: ["exon"] (sometimes also "CDS" depending on what you want)
    pub fn exon_feature_types(mut self, types: &[&str]) -> Self {
        self.keys.exon_feature_types = types.iter().map(|s| s.to_string()).collect();
        self
    }

    /// Build index from anything implementing `BufRead`.
    pub fn build_from_reader<R: BufRead>(&self, reader: R) -> Result<SpliceIndex, ParseError> {
        SpliceIndex::new(self.bin_width).from_reader(reader, self.keys.clone())
    }

    /// Build index from a file path.
    ///
    /// - If path ends with `.gz`, uses gzip decoder (requires the `flate2` dependency).
    /// - Otherwise reads as plain text.
    pub fn build_from_path<P: AsRef<Path>>(&self, path: P) -> Result<SpliceIndex, ParseError> {
        let path = path.as_ref();
        let file = std::fs::File::open(path).map_err(|e| ParseError::IoPath {
            path: path.display().to_string(),
            source: e,
        })?;

        // Decide gz vs plain by extension.
        let is_gz = path.extension().map(|e| e == "gz").unwrap_or(false);

        if is_gz {
            // NOTE: requires `flate2 = "1"` in Cargo.toml
            let decoder = flate2::read::GzDecoder::new(file);
            let reader = BufReader::new(decoder);
            self.build_from_reader(reader)
        } else {
            let reader = BufReader::new(file);
            self.build_from_reader(reader)
        }
    }
}

// -------------------- tests --------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn builder_gtf_default_keys_builds_index() {
        // Two exons of one transcript T1 belonging to gene G1.
        let gtf = "\
chr1\tsrc\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; gene_name \"Alpha\"; transcript_id \"T1\"; transcript_name \"TxA\";
chr1\tsrc\texon\t201\t250\t.\t+\t.\tgene_id \"G1\"; gene_name \"Alpha\"; transcript_id \"T1\"; transcript_name \"TxA\";
";

        let idx = AnnotationBuilder::new(100)
            .build_from_reader(Cursor::new(gtf.as_bytes()))
            .unwrap();

        assert_eq!(idx.chr_names, vec!["chr1".to_string()]);
        assert_eq!(idx.genes.len(), 1);
        assert_eq!(idx.transcripts.len(), 1);

        // We store multiple names/aliases; at minimum expect these to be present.
        assert!(idx.genes[0].names.iter().any(|n| n == "Alpha" || n == "G1"));
        assert!(idx.transcripts[0].names.iter().any(|n| n == "TxA" || n == "T1"));

        // Exons should be converted to 0-based half-open: 101..150 => [100,150)
        assert_eq!(idx.transcripts[0].exons()[0].start, 100);
        assert_eq!(idx.transcripts[0].exons()[0].end, 150);

        // Buckets exist
        assert_eq!(idx.chr_buckets.len(), 1);
        assert!(!idx.chr_buckets[0].bins.is_empty());
    }

    #[test]
    fn builder_gff3_parent_linking_builds_index() {
        // Minimal GFF3-like exon lines using Parent for transcript ID.
        // Here we intentionally do NOT include "transcript_id" etc.
        let gff = "\
chr2\tsrc\texon\t5\t20\t.\t-\t.\tParent=tx1;gene_id=G9;Name=GeneNice
chr2\tsrc\texon\t30\t40\t.\t-\t.\tParent=tx1;gene_id=G9;Name=GeneNice
";

        // Force transcript ID to come from Parent.
        let builder = AnnotationBuilder::new(50)
            .transcript_id_keys(&[]) // empty => will fall back to parent_keys in the index build
            .parent_keys(&["Parent"])
            .gene_id_key("gene_id")
            .gene_name_keys(&["Name"])
            .exon_feature_types(&["exon"]);

        let idx = builder
            .build_from_reader(Cursor::new(gff.as_bytes()))
            .unwrap();

        assert_eq!(idx.chr_names, vec!["chr2".to_string()]);
        assert_eq!(idx.genes.len(), 1);
        assert_eq!(idx.transcripts.len(), 1);

        assert_eq!(idx.transcripts[0].strand, crate::types::Strand::Minus);

        // GFF3 coordinates 5..20 => [4,20)
        assert_eq!(idx.transcripts[0].exons()[0].start, 4);
        assert_eq!(idx.transcripts[0].exons()[0].end, 20);

        // Gene/transcript names should include configured values
        assert!(idx.genes[0].names.iter().any(|n| n == "GeneNice" || n == "G9"));
        assert!(idx.transcripts[0].names.iter().any(|n| n == "tx1"));
    }

    #[test]
    fn builder_respects_exon_feature_types_filter() {
        // One exon and one CDS; we only want exon.
        let gtf = "\
chr1\tsrc\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
chr1\tsrc\tCDS\t201\t250\t.\t+\t0\tgene_id \"G1\"; transcript_id \"T1\";
";

        let builder = AnnotationBuilder::new(100).exon_feature_types(&["exon"]);

        let idx = builder
            .build_from_reader(Cursor::new(gtf.as_bytes()))
            .unwrap();

        assert_eq!(idx.transcripts.len(), 1);
        // Only the exon should be used => one block after finalize.
        assert_eq!(idx.transcripts[0].exons().len(), 1);
        assert_eq!(idx.transcripts[0].exons()[0].start, 100);
        assert_eq!(idx.transcripts[0].exons()[0].end, 150);
    }
}
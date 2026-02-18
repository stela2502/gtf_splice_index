use crate::model::types::{TranscriptId, GeneId};
use serde::{Serialize, Deserialize};

/// Gene model: stores one or more names/aliases and transcript ids.
///
/// Notes:
/// - `names[0]` is treated as the primary name (if present).
/// - additional names are aliases (deduped).
#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Gene {
    pub id: GeneId,
    pub names: Vec<String>,
    transcript_ids: Vec<TranscriptId>,
}

impl Gene {
    pub fn new(id: GeneId, primary_name: impl Into<String>) -> Self {
        Self {
            id,
            names: vec![primary_name.into()],
            transcript_ids: Vec::new(),
        }
    }

    /// Add an alias/alternative name (deduped).
    pub fn add_name(&mut self, name: &str) {
        let name = name.trim();
        if name.is_empty() {
            return;
        }
        if !self.names.iter().any(|n| n == name) {
            self.names.push(name.to_string());
        }
    }

    /// Primary name (if any).
    pub fn primary_name(&self) -> Option<&str> {
        self.names.first().map(|s| s.as_str())
    }

    pub fn add_transcript(&mut self, tx_id: TranscriptId) {
        self.transcript_ids.push(tx_id);
    }

    pub fn transcript_ids(&self) -> &[TranscriptId] {
        &self.transcript_ids
    }

    /// Sort transcript IDs and remove duplicates.
    pub fn finalize(&mut self) {
        self.transcript_ids.sort_unstable();
        self.transcript_ids.dedup();
    }

}



#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn names_are_deduped_and_primary_kept() {
        let mut g = Gene::new(0, "G1");
        g.add_name("G1");
        g.add_name("GeneSymbol");
        g.add_name("GeneSymbol");
        assert_eq!(g.names, vec!["G1".to_string(), "GeneSymbol".to_string()]);
        assert_eq!(g.primary_name(), Some("G1"));
    }

    #[test]
    fn transcripts_finalize_dedups() {
        let mut g = Gene::new(0, "G1");
        g.add_transcript(3);
        g.add_transcript(1);
        g.add_transcript(3);
        g.finalize();
        assert_eq!(g.transcript_ids(), &[1, 3]);
    }

    #[test]
    fn add_and_finalize_sorts_and_dedups() {
        let mut g = Gene::new(0, "G1");
        g.add_transcript(10);
        g.add_transcript(3);
        g.add_transcript(10);
        g.add_transcript(7);

        // before finalize: unsorted, may have duplicates
        assert_eq!(g.transcript_ids().len(), 4);

        g.finalize();

        assert_eq!(g.transcript_ids(), &[3, 7, 10]);
    }

    #[test]
    fn contains_transcript_works_even_without_finalize() {
        let mut g = Gene::new(1, "G2");
        g.add_transcript(5);
        g.add_transcript(2);

        assert_eq!(g.transcript_ids(), &[5,2]);

        // If you want stable ordering, finalize().
        g.finalize();
        assert_eq!(g.transcript_ids(), &[2, 5]);

    }
}


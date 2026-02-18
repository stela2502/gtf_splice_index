use serde::{Serialize, Deserialize};

/// Genomic strand/orientation.
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize )]
pub enum Strand {
    Plus,
    Minus,
    Unknown,
}

impl Strand {
    #[inline]
    pub fn is_compatible_with(self, other: Strand) -> bool {
        // "Unknown" is treated as compatible with either.
        self == Strand::Unknown || other == Strand::Unknown || self == other
    }
}

/// A contiguous genomic interval.
/// Coordinates are 0-based, half-open: [start, end)
#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize )]
pub struct RefBlock {
    pub start: u32,
    pub end: u32,
}

impl RefBlock {
    /// Create a new block. Panics if start >= end.
    pub fn new(start: u32, end: u32) -> Self {
        assert!(start < end, "RefBlock requires start < end");
        Self { start, end }
    }

    #[inline]
    pub fn len(self) -> u32 {
        self.end - self.start
    }

    #[inline]
    pub fn overlaps(self, other: RefBlock) -> bool {
        self.start < other.end && other.start < self.end
    }

    #[inline]
    pub fn contains(self, other: RefBlock) -> bool {
        self.start <= other.start && other.end <= self.end
    }
    
    /// Compute splice junctions from an ordered list of blocks.
    ///
    /// Junction i is: (blocks[i].end, blocks[i+1].start),
    /// except that small gaps (<= allowed_gap_size) are treated as sequencing/alignment artifacts
    /// and do NOT produce junctions.
    ///
    /// For unspliced reads (len < 2), this returns an empty vec.
    pub fn junctions_from_blocks(blocks: &[RefBlock], allowed_gap_size: u32) -> Vec<(u32, u32)> {
        let mut out = Vec::new();
        if blocks.len() < 2 {
            return out;
        }

        for w in blocks.windows(2) {
            let a = w[0];
            let b = w[1];

            // Overlapping or touching blocks: not a junction
            if b.start <= a.end {
                continue;
            }

            let gap = b.start - a.end;

            // Small gap inside an exon (sequencing/alignment artifact): ignore
            if gap <= allowed_gap_size {
                continue;
            }

            // Real splice junction
            out.push((a.end, b.start));
        }

        out
    }

}



/// A spliced read/query in genomic coordinates (BAM-independent).
///
/// Design:
/// - `chr_id` (usize) is the pre-mapped chromosome index used by the Index.
/// - `strand` is the read orientation (or Unknown).
/// - `blocks` are genomic blocks (0-based, half-open), usually exonic align blocks.
///
/// You can create it with unsorted blocks and then call `finalize()`
/// to sort and merge/clean them for stable matching.
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplicedRead {
    pub chr_id: usize,
    pub strand: Strand,
    pub blocks: Vec<RefBlock>,
    finalized: bool,
}

impl SplicedRead {
    pub fn new(chr_id: usize, strand: Strand, blocks: Vec<RefBlock>) -> Self {
        Self {
            chr_id,
            strand,
            blocks,
            finalized: false,
        }
    }

    /// Sort blocks by start and merge overlaps/adjacent.
    ///
    /// This makes matching stable and guards against messy upstream block creation.
    pub fn finalize(&mut self) {
        if self.blocks.is_empty() {
            self.finalized = true;
            return;
        }

        self.blocks.sort_by_key(|b| (b.start, b.end));

        let mut merged: Vec<RefBlock> = Vec::with_capacity(self.blocks.len());
        let mut cur = self.blocks[0];

        for &b in &self.blocks[1..] {
            // Merge if overlapping or adjacent.
            if b.start <= cur.end {
                cur.end = cur.end.max(b.end);
            } else {
                merged.push(cur);
                cur = b;
            }
        }
        merged.push(cur);

        self.blocks = merged;
        self.finalized = true;
    }

    /// Panics if not finalized (development-time safety).
    #[inline]
    pub fn assert_finalized(&self) {
        assert!(
            self.finalized,
            "SplicedRead must be finalized() before use in matching"
        );
    }

    /// Span of the read blocks (min start, max end).
    pub fn span(&self) -> Option<(u32, u32)> {
        if self.blocks.is_empty() {
            return None;
        }
        Some((self.blocks.first().unwrap().start, self.blocks.last().unwrap().end))
    }

    /// Junctions implied by the read blocks.
    pub fn junctions(&self) -> Vec<(u32, u32)> {
        RefBlock::junctions_from_blocks(&self.blocks, 0)
    }

}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn refblock_junctions() {
        let blocks = vec![RefBlock::new(10, 20), RefBlock::new(30, 40), RefBlock::new(50, 60)];
        assert_eq!(RefBlock::junctions_from_blocks(&blocks, 0), vec![(20, 30), (40, 50)]);
    }

    #[test]
    fn spliced_read_finalize_sorts_and_merges() {
        let mut sr = SplicedRead::new(
            1,
            Strand::Plus,
            vec![
                RefBlock::new(100, 120),
                RefBlock::new(110, 130), // overlap with first
                RefBlock::new(200, 210),
                RefBlock::new(210, 220), // adjacent
            ],
        );

        sr.finalize();
        assert_eq!(sr.blocks, vec![RefBlock::new(100, 130), RefBlock::new(200, 220)]);
        assert_eq!(sr.junctions(), vec![(130, 200)]);
        assert_eq!(sr.span(), Some((100, 220)));
    }
}


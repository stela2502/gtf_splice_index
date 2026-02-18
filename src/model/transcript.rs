use crate::model::types::{MatchClass, MatchHit, MatchOptions, TranscriptId};
use crate::types::{RefBlock, SplicedRead, Strand};
use serde::{Serialize, Deserialize};


#[derive(Debug, Clone, PartialEq, Eq, Serialize, Deserialize)]
pub struct Transcript {
    pub id: TranscriptId,
    pub gene_id: usize,
    pub names: Vec<String>,
    pub chr_id: usize,
    pub strand: Strand,
    exons: Vec<RefBlock>,
    finalized: bool,
}

impl Transcript {
    pub fn new(
        id: TranscriptId,
        gene_id: usize,
        primary_name: impl Into<String>,
        chr_id: usize,
        strand: Strand,
    ) -> Self {
        Self {
            id,
            gene_id,
            names: vec![primary_name.into()],
            chr_id,
            strand,
            exons: Vec::new(),
            finalized: false,
        }
    }

    pub fn add_name(&mut self, name: &str) {
        let name = name.trim();
        if name.is_empty() {
            return;
        }
        if !self.names.iter().any(|n| n == name) {
            self.names.push(name.to_string());
        }
    }

    pub fn primary_name(&self) -> Option<&str> {
        self.names.first().map(|s| s.as_str())
    }

    pub fn add_exon(&mut self, block: RefBlock) {
        self.exons.push(block);
        self.finalized = false;
    }

    pub fn exons(&self) -> &[RefBlock] {
        &self.exons
    }

    fn introns(&self) -> Vec<RefBlock> {
        let mut out = Vec::new();
        if self.exons.len() < 2 {
            return out;
        }
        for w in self.exons.windows(2) {
            let a = w[0];
            let b = w[1];
            // intron is [a.end, b.start)
            if a.end < b.start {
                out.push(RefBlock::new(a.end, b.start));
            }
        }
        out
    }
    /// sorts the transcripts exons and returns (total start: u32, total end: u32)
    pub fn finalize(&mut self) -> (u32, u32) {
        if self.exons.is_empty() {
            self.finalized = true;
            return (0,0);
        }

        self.exons.sort_by_key(|b| (b.start, b.end));

        let mut merged: Vec<RefBlock> = Vec::with_capacity(self.exons.len());
        let mut cur = self.exons[0];

        for &b in &self.exons[1..] {
            if b.start <= cur.end {
                cur.end = cur.end.max(b.end);
            } else {
                merged.push(cur);
                cur = b;
            }
        }
        merged.push(cur);

        self.exons = merged;
        self.finalized = true;

        let start = self.exons.first().unwrap().start;
        let end = self.exons.last().unwrap().end;

        (start, end)
    }

    pub fn span(&self) -> Option<(u32, u32)> {
        if self.exons.is_empty() {
            return None;
        }
        Some((self.exons.first().unwrap().start, self.exons.last().unwrap().end))
    }

    pub fn junctions(&self) -> Vec<(u32, u32)> {
        RefBlock::junctions_from_blocks(&self.exons, 0)
    }

    pub fn match_spliced_read(&self, read: &SplicedRead, opts: MatchOptions) -> MatchHit {
        read.assert_finalized();
        self.match_read_blocks(read.chr_id, read.strand, &read.blocks, opts)
    }
    fn overlaps(a0: u32, a1: u32, b0: u32, b1: u32) -> bool {
        a0 < b1 && b0 < a1
    }


    fn match_read_blocks(
        &self,
        read_chr_id: usize,
        read_strand: Strand,
        read_blocks: &[RefBlock],
        opts: MatchOptions,
    ) -> MatchHit {
        if !self.finalized {
            panic!("Transcript::match_read_blocks called before finalize()");
        }

        let mut over5 = 0u32;
        let mut over3 = 0u32;

        if read_blocks.is_empty() {
            return MatchHit::new(MatchClass::NoOverlap, 0, 0);
        }

        if self.chr_id != read_chr_id {
            return MatchHit::new(MatchClass::NoOverlap, 0, 0);
        }
        // Read span (blocks must be sorted + non-empty)
        let r0 = read_blocks[0].start;
        let r1 = read_blocks[read_blocks.len() - 1].end;

        let Some((t0, t1)) = self.span() else {
            return MatchHit::new(MatchClass::NoOverlap, 0, 0);
        };

        let read_span = RefBlock { start: r0, end: r1 };
        let tx_span = RefBlock { start: t0, end: t1 };

        if !read_span.overlaps(tx_span) {
            return MatchHit::new(MatchClass::NoOverlap, 0, 0);
        }

        if opts.require_strand && !self.strand.is_compatible_with(read_strand) {
            if let Some((o5, o3)) = self.compute_overhangs_strand_aware(read_blocks) {
                over5 = o5;
                over3 = o3;
            }
            return MatchHit::new(MatchClass::StrandMismatch, over5, over3);
        }
        
        if let Some((o5, o3)) = self.compute_overhangs_strand_aware(read_blocks) {
            over5 = o5;
            over3 = o3;
        }

        if over5 > opts.max_5p_overhang_bp || over3 > opts.max_3p_overhang_bp {
            return MatchHit::new(MatchClass::OverhangTooLarge, over5, over3);
        }

        if !self.blocks_fit_exons_allowing_end_overhang(read_blocks, 10) {
            return MatchHit::new(MatchClass::Intronic, over5, over3);
        }

        let read_junctions = RefBlock::junctions_from_blocks(read_blocks, opts.allowed_intronic_gap_size);
        let tx_junctions = self.junctions();

        if opts.require_exact_junction_chain {
            let class = if read_junctions == tx_junctions {
                MatchClass::ExactJunctionChain
            } else {
                MatchClass::JunctionMismatch
            };
            return MatchHit::new(class, over5, over3);
        }

        let class = if read_junctions.iter().all(|j| tx_junctions.contains(j)) {
            if read_junctions == tx_junctions {
                MatchClass::ExactJunctionChain
            } else {
                MatchClass::Compatible
            }
        } else {
            MatchClass::JunctionMismatch
        };

        MatchHit::new(class, over5, over3)
    }

    fn compute_overhangs_strand_aware(&self, read_blocks: &[RefBlock]) -> Option<(u32, u32)> {
        let (t0, t1) = self.span()?;
        let r0 = read_blocks.first()?.start;
        let r1 = read_blocks.last()?.end;

        let left_over = t0.saturating_sub(r0);
        let right_over = r1.saturating_sub(t1);

        match self.strand {
            Strand::Plus => Some((left_over, right_over)),
            Strand::Minus => Some((right_over, left_over)),
            Strand::Unknown => Some((left_over, right_over)),
        }
    }

    fn blocks_fit_exons_allowing_end_overhang(
        &self,
        read_blocks: &[RefBlock],
        max_in_exon_gap_bp: u32,
    ) -> bool {
        let mut exon_idx = 0usize;
        let num_blocks = read_blocks.len();

        // Track previous block and which exon it matched,
        // so we can validate gaps inside an exon.
        let mut prev_block_end: Option<u32> = None;
        let mut prev_exon_idx: Option<usize> = None;

        for (block_idx, &read_block) in read_blocks.iter().enumerate() {
            // Move exon pointer until we reach an exon whose end is beyond the block start.
            while exon_idx < self.exons.len() && self.exons[exon_idx].end <= read_block.start {
                exon_idx += 1;
            }
            if exon_idx == self.exons.len() {
                return false;
            }

            let exon = self.exons[exon_idx];

            // Block must overlap the exon it is assigned to.
            if !read_block.overlaps(exon) {
                return false;
            }

            let is_first_block = block_idx == 0;
            let is_last_block = block_idx + 1 == num_blocks;

            // For ANY block: it must not have intronic sequence on the "wrong side".
            // If it's the first block: allow 5' overhang, but NOT intronic 3' extension past exon end.
            // If it's the last block: allow 3' overhang, but NOT intronic 5' extension before exon start.
            // If it's an interior block: must be fully contained (you already do that).
            if is_first_block {
                if read_block.end > exon.end {
                    return false; // would include intronic bases after exon end
                }
            }
            if is_last_block {
                if read_block.start < exon.start {
                    return false; // would include intronic bases before exon start
                }
            }
            // Interior blocks must be fully inside the exon (no end-overhang allowed there).
            if !is_first_block && !is_last_block && !exon.contains(read_block) {
                return false;
            }

            // -----------------------------
            // NEW: allow small gaps inside the same exon
            // -----------------------------
            if let (Some(prev_end), Some(prev_ei)) = (prev_block_end, prev_exon_idx) {
                if read_block.start > prev_end {
                    let gap = read_block.start - prev_end;

                    // If two consecutive blocks map to the SAME exon, the gap must be small
                    // and fully internal to that exon. Otherwise it's a structural mismatch.
                    if prev_ei == exon_idx {
                        if gap > max_in_exon_gap_bp {
                            return false;
                        }

                        // Gap endpoints must lie within the exon bounds.
                        // (We allow equality at exon end/start since coordinates are half-open.)
                        let prev_end_in_exon = prev_end >= exon.start && prev_end <= exon.end;
                        let next_start_in_exon = read_block.start >= exon.start && read_block.start <= exon.end;

                        if !prev_end_in_exon || !next_start_in_exon {
                            return false;
                        }
                    }
                }
            }

            prev_block_end = Some(read_block.end);
            prev_exon_idx = Some(exon_idx);
        }

        true
    }
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::model::types::MatchOptions;

    // ---------- small helpers ----------
    fn opts() -> MatchOptions {
        MatchOptions {
            max_5p_overhang_bp: 15,
            max_3p_overhang_bp: 15,
            require_strand: false,
            require_exact_junction_chain: false,
            ..Default::default()
        }
    }

    fn tx_two_exons(chr_id: usize, strand: Strand) -> Transcript {
        let mut tx = Transcript::new(0, 0, "T", chr_id, strand);
        tx.add_exon(RefBlock::new(100, 150));
        tx.add_exon(RefBlock::new(200, 250));
        tx.finalize();
        tx
    }

    fn tx_three_exons(chr_id: usize, strand: Strand) -> Transcript {
        let mut tx = Transcript::new(0, 0, "T", chr_id, strand);
        tx.add_exon(RefBlock::new(100, 150));
        tx.add_exon(RefBlock::new(200, 250));
        tx.add_exon(RefBlock::new(300, 350));
        tx.finalize();
        tx
    }

    fn read(chr_id: usize, strand: Strand, blocks: Vec<RefBlock>) -> SplicedRead {
        let mut r = SplicedRead::new(chr_id, strand, blocks);
        r.finalize();
        r
    }

    // ---------- your existing tests ----------
    #[test]
    fn transcript_names_dedup() {
        let mut t = Transcript::new(0, 0, "T1", 0, Strand::Plus);
        t.add_name("T1");
        t.add_name("ENST0000");
        t.add_name("ENST0000");
        assert_eq!(t.names, vec!["T1".to_string(), "ENST0000".to_string()]);
        assert_eq!(t.primary_name(), Some("T1"));
    }

    #[test]
    fn exact_chain_and_overhang_reported() {
        let tx = tx_two_exons(1, Strand::Plus);

        let read = read(
            1,
            Strand::Plus,
            vec![RefBlock::new(90, 150), RefBlock::new(200, 260)],
        );

        let hit = tx.match_spliced_read(
            &read,
            MatchOptions { max_5p_overhang_bp: 15, max_3p_overhang_bp: 15, ..Default::default() },
        );

        assert_eq!(hit.class, MatchClass::ExactJunctionChain);
        assert_eq!(hit.overhang_5p_bp, 10);
        assert_eq!(hit.overhang_3p_bp, 10);
    }

    #[test]
    fn exact_chain_with_sequenceing_errors() {
        let tx = tx_two_exons(1, Strand::Plus);

        let read = read(
            1,
            Strand::Plus,
            vec![
                RefBlock::new(90, 150),
                RefBlock::new(200, 229),
                RefBlock::new(230, 260),
            ],
        );

        let hit = tx.match_spliced_read(
            &read,
            MatchOptions { max_5p_overhang_bp: 15, max_3p_overhang_bp: 15, allowed_intronic_gap_size:10, ..Default::default() },
        );

        assert_eq!(hit.class, MatchClass::ExactJunctionChain);
        assert_eq!(hit.overhang_5p_bp, 10);
        assert_eq!(hit.overhang_3p_bp, 10);
    }

    // ---------- new comprehensive coverage tests ----------

    #[test]
    fn matchclass_nooverlap_same_chr_far_away() {
        let tx = tx_two_exons(1, Strand::Plus);

        let read = read(1, Strand::Plus, vec![RefBlock::new(1000, 1050)]);

        let hit = tx.match_spliced_read(&read, opts());

        assert_eq!(hit.class, MatchClass::NoOverlap);
        assert_eq!(hit.overhang_5p_bp, 0);
        assert_eq!(hit.overhang_3p_bp, 0);
    }

    #[test]
    fn matchclass_strand_mismatch_when_required() {
        let tx = tx_two_exons(1, Strand::Plus);

        let read = read(
            1,
            Strand::Minus,
            vec![RefBlock::new(100, 150), RefBlock::new(200, 250)],
        );

        let mut o = opts();
        o.require_strand = true;

        let hit = tx.match_spliced_read(&read, o);

        assert_eq!(hit.class, MatchClass::StrandMismatch);
        // overhangs should still be reported (typically 0 here)
        assert_eq!(hit.overhang_5p_bp, 0);
        assert_eq!(hit.overhang_3p_bp, 0);
    }

    #[test]
    fn matchclass_overhang_too_large() {
        let tx = tx_two_exons(1, Strand::Plus);

        // 5' overhang = 30 (start 70 vs tx start 100)
        let read = read(
            1,
            Strand::Plus,
            vec![RefBlock::new(70, 150), RefBlock::new(200, 250)],
        );

        let hit = tx.match_spliced_read(&read, opts());

        assert_eq!(hit.class, MatchClass::OverhangTooLarge);
        assert_eq!(hit.overhang_5p_bp, 30);
        assert_eq!(hit.overhang_3p_bp, 0);
    }

    #[test]
    fn matchclass_intronic_due_to_intron_overlap() {
        let tx = tx_two_exons(1, Strand::Plus);

        // This second block overlaps the intron region [150,200) (and also a bit into exon2 start),
        // which should fail blocks_fit_exons_allowing_end_overhang(...) and be classified as Intronic.
        let read = read(
            1,
            Strand::Plus,
            vec![RefBlock::new(100, 150), RefBlock::new(170, 210)],
        );

        let hit = tx.match_spliced_read(&read, opts());

        assert_eq!(hit.class, MatchClass::Intronic);
    }

    #[test]
    fn matchclass_junction_mismatch_skipping_middle_exon() {
        let tx = tx_three_exons(1, Strand::Plus);

        // Read uses exon1 + exon3, skipping exon2.
        // Blocks still fit exons, but junction (150 -> 300) is not in tx junction set.
        let read = read(
            1,
            Strand::Plus,
            vec![RefBlock::new(100, 150), RefBlock::new(300, 350)],
        );

        let hit = tx.match_spliced_read(&read, opts());

        assert_eq!(hit.class, MatchClass::JunctionMismatch);
    }

    #[test]
    fn matchclass_compatible_when_junctions_are_subset() {
        let tx = tx_three_exons(1, Strand::Plus);

        // Read uses exon1 + exon2 only (junction subset, not equal to tx's full chain).
        let read = read(
            1,
            Strand::Plus,
            vec![RefBlock::new(100, 150), RefBlock::new(200, 250)],
        );

        let hit = tx.match_spliced_read(&read, opts());

        assert_eq!(hit.class, MatchClass::Compatible);
    }

    #[test]
    fn require_exact_junction_chain_rejects_subset_as_junction_mismatch() {
        let tx = tx_three_exons(1, Strand::Plus);

        // Same read as Compatible case, but now exact is required.
        let read = read(
            1,
            Strand::Plus,
            vec![RefBlock::new(100, 150), RefBlock::new(200, 250)],
        );

        let mut o = opts();
        o.require_exact_junction_chain = true;

        let hit = tx.match_spliced_read(&read, o);

        assert_eq!(hit.class, MatchClass::JunctionMismatch);
    }
}

use std::cmp::Ordering;
use std::fmt;

/// Internal numeric IDs (indexes into Vecs).
pub type GeneId = usize;
pub type TranscriptId = usize;

/// Classification of how a spliced read matches a transcript model.
///
/// Overhangs are *not* part of this enum; they are always reported via `MatchHit`.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum MatchClass {
    /// Read blocks/junctions are fully consistent with the transcript model.
    Compatible,

    /// Read’s splice junction chain exactly equals the transcript’s junction chain.
    ExactJunctionChain,

    /// Read overlaps transcript span but places sequence in intronic regions.
    Intronic,

    /// Read overlaps transcript span but contains junctions not in the transcript.
    JunctionMismatch,

    /// Read exceeds allowed 5′ or 3′ overhang tolerance.
    OverhangTooLarge,

    /// Read does not overlap transcript span at all.
    NoOverlap,

    /// Strand was incompatible (unless unknown is allowed).
    StrandMismatch,
}


impl fmt::Display for MatchClass {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        let s = match self {
            MatchClass::Compatible => "Compatible",
            MatchClass::ExactJunctionChain => "ExactJunctionChain",
            MatchClass::Intronic => "Intronic",
            MatchClass::JunctionMismatch => "JunctionMismatch",
            MatchClass::OverhangTooLarge => "OverhangTooLarge",
            MatchClass::NoOverlap => "NoOverlap",
            MatchClass::StrandMismatch => "StrandMismatch",
        };
        write!(f, "{s}")
    }
}

impl MatchClass {
    /// Numeric ranking used for comparisons.
    /// Higher is better.
    pub fn rank(self) -> u8 {
        match self {
            MatchClass::ExactJunctionChain => 6,
            MatchClass::Compatible => 5,
            MatchClass::JunctionMismatch => 3,
            MatchClass::Intronic => 2,
            MatchClass::OverhangTooLarge => 1,
            MatchClass::StrandMismatch => 0,
            MatchClass::NoOverlap => 0,
        }
    }
}

impl Ord for MatchClass {
    fn cmp(&self, other: &Self) -> Ordering {
        self.rank().cmp(&other.rank())
    }
}

impl PartialOrd for MatchClass {
    fn partial_cmp(&self, other: &Self) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

/// Convenience category for overhang reporting.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum OverhangClass {
    None,
    Overhang5p,
    Overhang3p,
    Both,
}

/// Always returned from matching:
/// - `class`: the core classification
/// - `overhang_5p_bp` and `overhang_3p_bp`: always filled (0 means none)
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MatchHit {
    pub class: MatchClass,
    pub overhang_5p_bp: u32,
    pub overhang_3p_bp: u32,
}

impl MatchHit {
    pub fn new(class: MatchClass, overhang_5p_bp: u32, overhang_3p_bp: u32) -> Self {
        Self {
            class,
            overhang_5p_bp,
            overhang_3p_bp,
        }
    }

    pub fn overhang_class(&self) -> OverhangClass {
        match (self.overhang_5p_bp > 0, self.overhang_3p_bp > 0) {
            (false, false) => OverhangClass::None,
            (true, false) => OverhangClass::Overhang5p,
            (false, true) => OverhangClass::Overhang3p,
            (true, true) => OverhangClass::Both,
        }
    }
}

/// Options controlling what "match" means.
///
/// IMPORTANT:
/// - Overhangs are always *measured* and returned.
/// - These options only control whether the match is still considered valid
///   when overhangs exist.
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct MatchOptions {
    /// If true, require read blocks to be on a compatible strand.
    pub require_strand: bool,

    /// If true, require the read to have the exact same splice junction chain as the transcript.
    pub require_exact_junction_chain: bool,

    /// Maximum allowed 5′ overhang (bp). If exceeded -> OverhangTooLarge.
    pub max_5p_overhang_bp: u32,

    /// Maximum allowed 3′ overhang (bp). If exceeded -> OverhangTooLarge.
    pub max_3p_overhang_bp: u32,

    /// allowed sequenceing error gap. If exceed -> JunctionMismatch
    pub allowed_intronic_gap_size: u32,
}

impl Default for MatchOptions {
    fn default() -> Self {
        Self {
            require_strand: true,
            require_exact_junction_chain: false,
            max_5p_overhang_bp: 0,
            max_3p_overhang_bp: 0,
            allowed_intronic_gap_size: 0,
        }
    }
}


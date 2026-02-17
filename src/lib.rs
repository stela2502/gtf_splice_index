//! gtf_splice_index
//!
//! Pure-Rust, BAM-independent splice matching library.
//! This crate models transcripts/genes and matches spliced reads represented
//! as genomic blocks (0-based, half-open).

pub mod types;
pub mod model;
pub mod annotation;
pub mod index;

pub use index::{IdNameKeys, SpliceIndex};

pub use annotation::AnnotationBuilder;

pub use types::{RefBlock, SplicedRead, Strand, };

pub use model::transcript::{Transcript};
pub use model::gene::{Gene};
pub use model::types::{GeneId, TranscriptId};

// Re-export shared model matching types at crate root.
pub use model::{MatchClass, MatchOptions, OverhangClass};

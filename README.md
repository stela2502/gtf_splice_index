# gtf_splice_index

A flexible, streaming **GTF/GFF3 parser** plus a **splice-aware transcript index**
for fast matching of spliced reads against transcript models.

This crate is designed to sit between:

- **Alignment layer** (e.g. BAM → spliced blocks)
- **Annotation layer** (GTF/GFF transcripts)

and provide fast, deterministic transcript matching.

Coordinates are **0-based, half-open**: `[start, end)`.

---

## Main workflow

1. Build a `SpliceIndex` from a GTF/GFF file.
2. Convert an aligned read into a `SplicedRead`.
3. Call `match_transcripts()` to obtain best transcript matches.

This is the intended **public entry point** of the crate.

---

## Core exported types

From the crate root:

```rust
pub use types::{RefBlock, SplicedRead, Strand};
```

### `SplicedRead`

```rust
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct SplicedRead {
    pub chr_id: usize,
    pub strand: Strand,
    pub blocks: Vec<RefBlock>,
    finalized: bool,
}
```

A `SplicedRead` is simply a list of aligned reference blocks.

Typically these come from a BAM CIGAR string.

---

## Building the index from a GTF/GFF

High-level API:

```rust
use gtf_splice_index::AnnotationBuilder;

// 1 Mb bins: good default for mammalian genomes
let index = AnnotationBuilder::new(1_000_000)
    .build_from_path("genes.gtf")
    .unwrap();
```

### Bin size parameter

The value passed to `AnnotationBuilder::new()` controls the **genomic bin size**
used for transcript bucketing.

Recommended values:

| Genome type | Bin size |
|-------------|----------|
| Human / mouse | 1,000,000 |
| Fly / yeast | 100,000 |
| Bacteria | 10,000–50,000 |

For most users: **1,000,000 is a safe default**.

---

## Streaming parser (advanced usage)

If you need low-level access to annotation records:

```rust
use std::fs::File;
use std::io::BufReader;
use gtf_splice_index::annotation::io::AnnotationReader;

let file = File::open("genes.gtf").unwrap();
let reader = BufReader::new(file);

let rdr = AnnotationReader::new(reader);
for rec in rdr.records() {
    let rec = rec.unwrap();
    println!("{} {}-{}", rec.seqname, rec.start0, rec.end0);
}
```

Most users should use `AnnotationBuilder` instead.

---

## MatchOptions: controlling what “match” means

Matching behavior is configured using `MatchOptions`.

```rust
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

    /// Allowed sequencing micro-gap inside exons.
    /// Gaps ≤ this size are ignored when computing junctions.
    /// Larger gaps become real junctions and may lead to JunctionMismatch.
    pub allowed_intronic_gap_size: u32,
}
```

### Semantics

- Overhangs are **always measured** and returned.
- Options only control whether a read is considered valid.
- Small internal gaps can be tolerated to avoid false junction mismatches.

---

## Matching a read against transcripts

The main API:

```rust
let hits = index.match_transcripts(&read, opts);
```

### Complete example

```rust
use gtf_splice_index::{
    AnnotationBuilder,
    RefBlock,
    SplicedRead,
    Strand,
    MatchOptions,
};

// 1) Build index
let index = AnnotationBuilder::new(1_000_000)
    .build_from_path("genes.gtf")
    .unwrap();

// 2) Build spliced read
let mut read = SplicedRead::new(
    1,
    Strand::Plus,
    vec![
        RefBlock::new(100, 150),
        RefBlock::new(200, 250),
    ],
);
read.finalize();

// 3) Configure matching
let opts = MatchOptions {
    require_strand: true,
    require_exact_junction_chain: false,
    max_5p_overhang_bp: 10,
    max_3p_overhang_bp: 10,
    allowed_intronic_gap_size: 5,
};

// 4) Match
let hits = index.match_transcripts(&read, opts);

for hit in hits {
    println!(
        "transcript={} class={:?} over5={} over3={}",
        hit.transcript.primary_name().unwrap_or("unknown"),
        hit.hit.class,
        hit.hit.overhang_5p_bp,
        hit.hit.overhang_3p_bp,
    );
}
```

---

## Match classes (conceptual overview)

Typical classifications:

- `ExactJunctionChain`  
  Read junction chain equals transcript junction chain.

- `Compatible`  
  Read junctions are a subset of transcript junctions.

- `JunctionMismatch`  
  Read contains junction(s) not present in transcript.

- `Intronic`  
  Read overlaps transcript span but includes intronic sequence.

- `OverhangTooLarge`  
  End overhang exceeds configured threshold.

- `StrandMismatch`  
  Strand incompatible (if required).

- `NoOverlap`  
  No genomic overlap with transcript.

The returned `MatchHit` always includes measured 5′ and 3′ overhangs.

---

## Typical integration pattern

1. Parse BAM record → build `SplicedRead`
2. Call `index.match_transcripts()`
3. Select best hit(s) by match class and overhangs
4. Assign transcript or gene label

---

## Performance characteristics

- Streaming annotation parser
- Binned transcript index for fast candidate lookup
- Junction-based transcript filtering
- Deterministic matching results

---

## License

Add your license here (e.g. MIT / Apache-2.0).

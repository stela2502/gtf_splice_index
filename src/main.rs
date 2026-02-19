use std::path::PathBuf;

use anyhow::{Context, Result};
use clap::{Args, Parser, Subcommand};

use gtf_splice_index::{IdNameKeys, SpliceIndex}; // <-- adjust crate path/module as needed

/// Build, inspect, or serialize a splice index.
#[derive(Parser, Debug)]
#[command(name = "splice-index")]
#[command(author, version, about)]
struct Cli {
    #[command(subcommand)]
    cmd: Command,
}

#[derive(Subcommand, Debug)]
enum Command {
    /// Build an index from a GTF/GFF annotation and write it to disk
    Build(BuildArgs),

    /// Load an index from disk and print summary stats
    Stats(StatsArgs),
}

#[derive(Args, Debug)]
struct StatsArgs {
    /// Serialized index file
    #[arg(long, short)]
    index: PathBuf,
}

#[derive(Args, Debug)]
struct BuildArgs {
    /// Input annotation file (.gtf/.gff/.gff3)
    #[arg(long, short)]
    annotation: PathBuf,

    /// Bin width in base pairs
    #[arg(long, short, default_value_t = 1_000_000)]
    bin_width: u32,

    /// Output serialized index file
    #[arg(long, short)]
    index: PathBuf,

    // -------------------------
    // Attribute key options
    // -------------------------

    /// Attribute keys to use for gene ID (repeatable).
    /// Default (GTF-safe): gene_id
    #[arg(
        long = "gene-id-key",
        value_name = "KEY",
        num_args = 1..,
        default_values_t = vec!["gene_id".to_string()]
    )]
    gene_id_keys: Vec<String>,

    /// Attribute keys to use for gene name (repeatable).
    /// Default (GTF-safe): gene_name
    #[arg(
        long = "gene-name-key",
        value_name = "KEY",
        num_args = 1..,
        default_values_t = vec!["gene_name".to_string()]
    )]
    gene_name_keys: Vec<String>,

    /// Attribute keys to use for transcript ID (repeatable).
    /// Default (GTF-safe): transcript_id
    #[arg(
        long = "transcript-id-key",
        value_name = "KEY",
        num_args = 1..,
        default_values_t = vec!["transcript_id".to_string()]
    )]
    transcript_id_keys: Vec<String>,

    /// Attribute keys to use for transcript name (repeatable).
    /// Default (GTF-safe): transcript_name
    #[arg(
        long = "transcript-name-key",
        value_name = "KEY",
        num_args = 1..,
        default_values_t = vec!["transcript_name".to_string()]
    )]
    transcript_name_keys: Vec<String>,

    /// GFF3 exon->transcript linkage keys (repeatable).
    /// Default: Parent
    #[arg(
        long = "parent-key",
        value_name = "KEY",
        num_args = 1..,
        default_values_t = vec!["Parent".to_string()]
    )]
    parent_keys: Vec<String>,

    /// Feature types that count as exon blocks (repeatable).
    /// Default: exon
    #[arg(
        long = "exon-feature-type",
        value_name = "TYPE",
        num_args = 1..,
        default_values_t = vec!["exon".to_string()]
    )]
    exon_feature_types: Vec<String>,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    match cli.cmd {
        Command::Build(args) => {
            let keys = IdNameKeys {
                gene_id_keys: args.gene_id_keys,
                gene_name_keys: args.gene_name_keys,
                transcript_id_keys: args.transcript_id_keys,
                transcript_name_keys: args.transcript_name_keys,
                parent_keys: args.parent_keys,
                exon_feature_types: args.exon_feature_types,
            };

            let idx = SpliceIndex::from_path(&args.annotation, args.bin_width, keys)
                .with_context(|| format!("building index from {}", args.annotation.display()))?;

            println!("{idx}");

            idx.save(&args.index)
                .with_context(|| format!("writing index to {}", args.index.display()))?;

            eprintln!("Index written to {}", args.index.display());
        }

        Command::Stats(args) => {
            let idx = SpliceIndex::load(&args.index)
                .with_context(|| format!("reading index {}", args.index.display()))?;
            println!("{idx}");
        }
    }

    Ok(())
}

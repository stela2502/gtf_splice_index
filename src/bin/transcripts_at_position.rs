use anyhow::{Context, Result, anyhow};
use clap::Parser;
use gtf_splice_index::SpliceIndex;
use gtf_splice_index::IdNameKeys;

#[derive(Debug, Parser)]
#[command(
    author,
    version,
    about = "Report transcripts whose exons cover genomic positions"
)]
struct Cli {

    /// Build index from GTF/GFF annotation
    #[arg(long, conflicts_with = "index", required_unless_present = "index")]
    annotation: Option<std::path::PathBuf>,

    /// Load prebuilt binary splice index
    #[arg(long, conflicts_with = "annotation", required_unless_present = "annotation")]
    index: Option<std::path::PathBuf>,

    /// Bin width used when building the index
    #[arg(long, default_value_t = 1_000_000)]
    bin_width: u32,

    /// Chromosome/contig name, e.g. chr17
    #[arg(short, long)]
    chr: String,

    /// 1-based genomic position
    #[arg(short, long)]
    pos: u32,
}

fn main() -> Result<()> {
    let cli = Cli::parse();

    if cli.pos == 0 {
        return Err(anyhow!("--pos must be 1-based and > 0"));
    }

    let pos0 = cli.pos - 1;

    let index = if let Some(path) = &cli.index {
        SpliceIndex::load(path)
            .with_context(|| format!("failed to load binary index: {}", path.display()))?
    } else {
        SpliceIndex::from_path(
            cli.annotation.as_ref().unwrap(),
            cli.bin_width,
            IdNameKeys::default(),
        )
        .with_context(|| {
            format!(
                "failed to build splice index from annotation: {}",
                cli.annotation.as_ref().unwrap().display()
            )
        })?
    };

    println!("Index loaded");

    let chr_id = index
        .chr_id(&cli.chr)
        .with_context(|| format!("chromosome not found in index: {}", cli.chr))?;

    let tx_ids = index.candidates_for_span_union(chr_id, pos0, pos0 + 1);

    println!(
        "query_chr\tquery_pos1\tgene_id\tgene_names\ttranscript_id\ttranscript_names\tstrand\texon_number\texon_start1\texon_end1"
    );

    for tx_id in tx_ids {
        let tx = &index.transcripts[tx_id as usize];

        for (exon_idx, exon) in tx.exons().iter().enumerate() {
            if exon.start <= pos0 && pos0 < exon.end {
                let gene = &index.genes[tx.gene_id];

                println!(
                    "{}\t{}\t{}\t{}\t{}\t{}\t{:?}\t{}\t{}\t{}",
                    cli.chr,
                    cli.pos,
                    tx.gene_id,
                    gene.names.join(","),
                    tx_id,
                    tx.names.join(","),
                    tx.strand,
                    exon_idx + 1,
                    exon.start + 1,
                    exon.end,
                );
            }
        }
    }

    Ok(())
}
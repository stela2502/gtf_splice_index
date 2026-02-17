use std::collections::HashMap;
use std::io::BufRead;

use crate::types::Strand;

/// File dialect detected from attribute syntax.
///
/// - GFF3 typically uses: key=value;key2=value2
/// - GTF typically uses: key "value"; key2 "value2";
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub enum Dialect {
    Gff3,
    Gtf,
    Unknown,
}

/// A single parsed record line from GTF/GFF3.
///
/// Coordinates:
/// - `start0` is 0-based start
/// - `end0` is 0-based end (half-open)
#[derive(Debug, Clone, PartialEq)]
pub struct AnnotationRecord {
    pub seqname: String,      // chromosome / contig
    pub source: String,       // column 2
    pub feature_type: String, // column 3
    pub start0: u32,          // 0-based start
    pub end0: u32,            // 0-based end (half-open)
    pub score: Option<f32>,   // '.' => None
    pub strand: Strand,       // + / - / . / ?
    pub phase: Option<u8>,    // '.' => None, else 0/1/2
    pub attrs: HashMap<String, String>,
    pub dialect: Dialect,
}

impl AnnotationRecord {
    /// Convenience: get an attribute value.
    pub fn attr(&self, key: &str) -> Option<&str> {
        self.attrs.get(key).map(|s| s.as_str())
    }


    pub fn is_exon_feature(&self, exon_types: &[String]) -> bool {
        exon_types.iter().any(|t| t == &self.feature_type)
    }

    pub fn pick_first_attr(&self, keys: &[String]) -> Option<String> {
        for k in keys {
            if let Some(v) = self.attr(k) {
                let v = v.trim();
                if !v.is_empty() {
                    return Some(v.to_string());
                }
            }
        }
        None
    }
}

/// Parsing errors for GTF/GFF3.
#[derive(Debug)]
pub enum ParseError {
    IoPath { path: String, source: std::io::Error },
    MalformedLine { line: String },
    BadCoordinates { line: String },
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::IoPath { path, source } => {
                write!(f, "I/O error while reading '{}': {}", path, source)
            }
            ParseError::MalformedLine { line } => write!(f, "Malformed GTF/GFF line: {}", line),
            ParseError::BadCoordinates { line } => write!(f, "Bad coordinates in line: {}", line),
        }
    }
}

impl std::error::Error for ParseError {
    fn source(&self) -> Option<&(dyn std::error::Error + 'static)> {
        match self {
            ParseError::IoPath { source, .. } => Some(source),
            _ => None,
        }
    }
}


/// Low-level streaming parser for GTF/GFF3 files.
///
/// Most users should **not** use this directly.  
/// Instead, use [`crate::annotation::AnnotationBuilder`] to build a full
/// `SpliceIndex` from a file in one step.
///
/// # Example (high-level usage)
/// ```no_run
/// use gtf_splice_index::AnnotationBuilder;
///
/// let index = AnnotationBuilder::new(100)
///     .build_from_path("genes.gtf")
///     .unwrap();
/// ```
///
/// # Example (streaming low-level usage)
/// ```no_run
/// use std::fs::File;
/// use std::io::BufReader;
/// use gtf_splice_index::annotation::io::AnnotationReader;
///
/// let file = File::open("genes.gtf").unwrap();
/// let reader = BufReader::new(file);
///
/// let rdr = AnnotationReader::new(reader);
/// for rec in rdr.records() {
///     let rec = rec.unwrap();
///     println!("{} {}-{}", rec.seqname, rec.start0, rec.end0);
/// }
/// ```
pub struct AnnotationReader<R: BufRead> {
    reader: R,
    buf: String,
}

impl<R: BufRead> AnnotationReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            buf: String::new(),
        }
    }

    /// Returns an iterator over parsed records.
    ///
    /// - Skips blank lines
    /// - Skips comment lines starting with '#'
    pub fn records(mut self) -> impl Iterator<Item = Result<AnnotationRecord, ParseError>> {
        std::iter::from_fn(move || loop {
            self.buf.clear();
            match self.reader.read_line(&mut self.buf) {
                Ok(0) => return None,
                Ok(_) => {}
                Err(e) => {
                    return Some(Err(ParseError::IoPath {
                        path: "<reader>".to_string(),
                        source: e,
                    }))
                },
            }

            let line = self.buf.trim_end_matches(&['\n', '\r'][..]);
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            return Some(parse_record_line(line));
        })
    }
}

/// Parse a single non-comment line into an `AnnotationRecord`.
pub fn parse_record_line(line: &str) -> Result<AnnotationRecord, ParseError> {
    // GTF/GFF have 9 tab-separated columns:
    // seqname source feature start end score strand phase attributes
    let mut it = line.split('\t');

    let seqname = it.next().ok_or_else(|| ParseError::MalformedLine {
        line: line.to_string(),
    })?;
    let source = it.next().ok_or_else(|| ParseError::MalformedLine {
        line: line.to_string(),
    })?;
    let feature_type = it.next().ok_or_else(|| ParseError::MalformedLine {
        line: line.to_string(),
    })?;
    let start_s = it.next().ok_or_else(|| ParseError::MalformedLine {
        line: line.to_string(),
    })?;
    let end_s = it.next().ok_or_else(|| ParseError::MalformedLine {
        line: line.to_string(),
    })?;
    let score_s = it.next().ok_or_else(|| ParseError::MalformedLine {
        line: line.to_string(),
    })?;
    let strand_s = it.next().ok_or_else(|| ParseError::MalformedLine {
        line: line.to_string(),
    })?;
    let phase_s = it.next().ok_or_else(|| ParseError::MalformedLine {
        line: line.to_string(),
    })?;
    let attrs_s = it.next().ok_or_else(|| ParseError::MalformedLine {
        line: line.to_string(),
    })?;

    // no extra columns
    if it.next().is_some() {
        return Err(ParseError::MalformedLine {
            line: line.to_string(),
        });
    }

    // Coordinates: input is 1-based inclusive; convert to 0-based half-open [start-1, end)
    let start_1: u64 = start_s.parse().map_err(|_| ParseError::BadCoordinates {
        line: line.to_string(),
    })?;
    let end_1: u64 = end_s.parse().map_err(|_| ParseError::BadCoordinates {
        line: line.to_string(),
    })?;

    if start_1 == 0 || end_1 == 0 || end_1 < start_1 {
        return Err(ParseError::BadCoordinates {
            line: line.to_string(),
        });
    }

    let start0 = (start_1 - 1) as u32;
    let end0 = end_1 as u32; // inclusive -> half-open end is end_1

    let score = if score_s == "." {
        None
    } else {
        Some(score_s.parse::<f32>().map_err(|_| ParseError::MalformedLine {
            line: line.to_string(),
        })?)
    };

    let strand = match strand_s {
        "+" => Strand::Plus,
        "-" => Strand::Minus,
        "." | "?" => Strand::Unknown,
        _ => {
            return Err(ParseError::MalformedLine {
                line: line.to_string(),
            })
        }
    };

    let phase = if phase_s == "." {
        None
    } else {
        let p: u8 = phase_s.parse().map_err(|_| ParseError::MalformedLine {
            line: line.to_string(),
        })?;
        if p > 2 {
            return Err(ParseError::MalformedLine {
                line: line.to_string(),
            });
        }
        Some(p)
    };

    let (dialect, attrs) = parse_attributes(attrs_s);

    Ok(AnnotationRecord {
        seqname: seqname.to_string(),
        source: source.to_string(),
        feature_type: feature_type.to_string(),
        start0,
        end0,
        score,
        strand,
        phase,
        attrs,
        dialect,
    })
}

/// Parse the attributes field for either GFF3 or GTF.
///
/// Returns (Dialect, map).
///
/// Heuristics:
/// - If it contains '=' => treat as GFF3
/// - Else if it contains quotes => treat as GTF
/// - Else Unknown, but parse best-effort
pub fn parse_attributes(s: &str) -> (Dialect, HashMap<String, String>) {
    let s = s.trim();

    let dialect = if s.contains('=') {
        Dialect::Gff3
    } else if s.contains('"') {
        Dialect::Gtf
    } else {
        Dialect::Unknown
    };

    let mut map = HashMap::new();

    match dialect {
        Dialect::Gff3 => {
            for part in s.split(';') {
                let part = part.trim();
                if part.is_empty() {
                    continue;
                }
                let mut kv = part.splitn(2, '=');
                let k = kv.next().unwrap().trim();
                let v = kv.next().unwrap_or("").trim();
                if !k.is_empty() {
                    map.insert(k.to_string(), unquote(v));
                }
            }
        }
        Dialect::Gtf => {
            // Split by ';' then parse key + rest (quoted or not)
            for part in s.split(';') {
                let part = part.trim();
                if part.is_empty() {
                    continue;
                }
                let mut it = part.splitn(2, char::is_whitespace);
                let key = it.next().unwrap().trim();
                let rest = it.next().unwrap_or("").trim();

                if key.is_empty() {
                    continue;
                }
                let value = unquote(rest);
                if !value.is_empty() {
                    map.insert(key.to_string(), value);
                }
            }
        }
        Dialect::Unknown => {
            for part in s.split(';') {
                let part = part.trim();
                if part.is_empty() {
                    continue;
                }
                if part.contains('=') {
                    let mut kv = part.splitn(2, '=');
                    let k = kv.next().unwrap().trim();
                    let v = kv.next().unwrap_or("").trim();
                    if !k.is_empty() {
                        map.insert(k.to_string(), unquote(v));
                    }
                } else {
                    let mut it = part.splitn(2, char::is_whitespace);
                    let key = it.next().unwrap().trim();
                    let rest = it.next().unwrap_or("").trim();
                    if !key.is_empty() && !rest.is_empty() {
                        map.insert(key.to_string(), unquote(rest));
                    }
                }
            }
        }
    }

    (dialect, map)
}

fn unquote(v: &str) -> String {
    let v = v.trim();
    let v = v.strip_prefix('"').unwrap_or(v);
    let v = v.strip_suffix('"').unwrap_or(v);
    v.to_string()
}

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn parse_gtf_line() {
        let line = "chr1\tsrc\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_number \"1\";";
        let rec = parse_record_line(line).unwrap();

        assert_eq!(rec.dialect, Dialect::Gtf);
        assert_eq!(rec.seqname, "chr1");
        assert_eq!(rec.feature_type, "exon");
        // 101..150 inclusive -> [100,150)
        assert_eq!(rec.start0, 100);
        assert_eq!(rec.end0, 150);
        assert_eq!(rec.strand, Strand::Plus);

        assert_eq!(rec.attr("gene_id"), Some("G1"));
        assert_eq!(rec.attr("transcript_id"), Some("T1"));
        assert_eq!(rec.attr("exon_number"), Some("1"));
    }

    #[test]
    fn parse_gff3_line() {
        let line = "chr2\tsrc\texon\t5\t20\t.\t-\t.\tID=ex1;Parent=tx1;gene_id=G9";
        let rec = parse_record_line(line).unwrap();

        assert_eq!(rec.dialect, Dialect::Gff3);
        assert_eq!(rec.seqname, "chr2");
        assert_eq!(rec.feature_type, "exon");
        // 5..20 inclusive -> [4,20)
        assert_eq!(rec.start0, 4);
        assert_eq!(rec.end0, 20);
        assert_eq!(rec.strand, Strand::Minus);

        assert_eq!(rec.attr("ID"), Some("ex1"));
        assert_eq!(rec.attr("Parent"), Some("tx1"));
        assert_eq!(rec.attr("gene_id"), Some("G9"));
    }

    #[test]
    fn streaming_reader_skips_comments_and_blank_lines() {
        let data = "\
#comment
chr1\tsrc\texon\t1\t2\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
\n
chr1\tsrc\texon\t3\t4\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\";
";
        let cur = Cursor::new(data.as_bytes());
        let reader = AnnotationReader::new(cur);

        let recs: Vec<_> = reader.records().collect::<Result<Vec<_>, _>>().unwrap();
        assert_eq!(recs.len(), 2);
        assert_eq!(recs[0].start0, 0);
        assert_eq!(recs[0].end0, 2);
        assert_eq!(recs[1].start0, 2);
        assert_eq!(recs[1].end0, 4);
    }
}

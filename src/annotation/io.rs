// io.rs
//
// Low-level streaming parser for GTF/GFF3 files.
// - Produces `AnnotationRecord`s
// - Provides *useful* errors (what we tried, what failed, preview of the line)
//
// This module is intended to be used by higher-level builders (e.g. AnnotationBuilder / SpliceIndex).

use std::collections::HashMap;
use std::io::BufRead;

use crate::types::Strand;

// ----------------------------
// Dialect
// ----------------------------

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

// ----------------------------
// ParseError
// ----------------------------

/// Parsing errors for GTF/GFF3.
///
/// NOTE:
/// - We keep this focused on I/O + format-level issues.
/// - Missing *required* attributes (gene_id/transcript_id) should be handled where you
///   require them (i.e. in the index builder), because that's a policy choice.
#[derive(Debug)]
pub enum ParseError {
    IoPath { path: String, source: std::io::Error },

    /// "We tried to parse X, but failed because Y."
    MalformedLine {
        line_no: Option<usize>,
        expected: &'static str,
        problem: String,
        details: Option<String>,
        line_preview: String,
    },

    BadCoordinates {
        line_no: Option<usize>,
        problem: String,
        line_preview: String,
    },
}

impl std::fmt::Display for ParseError {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        match self {
            ParseError::IoPath { path, source } => {
                write!(f, "I/O error while reading '{}': {}", path, source)
            }

            ParseError::MalformedLine {
                line_no,
                expected,
                problem,
                details,
                line_preview,
            } => {
                if let Some(n) = line_no {
                    write!(f, "Malformed annotation line (line {}): ", n)?;
                } else {
                    write!(f, "Malformed annotation line: ")?;
                }
                write!(f, "tried to parse {}, but failed: {}", expected, problem)?;
                if let Some(d) = details {
                    write!(f, " ({})", d)?;
                }
                write!(f, "\n  line: {}", line_preview)
            }

            ParseError::BadCoordinates {
                line_no,
                problem,
                line_preview,
            } => {
                if let Some(n) = line_no {
                    write!(f, "Bad coordinates (line {}): {}", n, problem)?;
                } else {
                    write!(f, "Bad coordinates: {}", problem)?;
                }
                write!(f, "\n  line: {}", line_preview)
            }
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

// ----------------------------
// AnnotationRecord
// ----------------------------

/// A single parsed record line from GTF/GFF3.
///
/// Coordinates:
/// - `start0` is 0-based start
/// - `end0` is 0-based end (half-open)
#[derive(Debug, Clone, PartialEq)]
pub struct AnnotationRecord {
    /// 1-based line number in the input (best-effort from the reader).
    pub line_no: usize,

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

    /// Short preview of the original line for better error messages/debugging.
    pub line_preview: String,
}

impl AnnotationRecord {
    /// Convenience: get an attribute value.
    pub fn attr(&self, key: &str) -> Option<&str> {
        self.attrs.get(key).map(|s| s.as_str())
    }

    /// Feature-type filter helper.
    pub fn is_exon_feature(&self, exon_types: &[String]) -> bool {
        exon_types.iter().any(|t| t == &self.feature_type)
    }

    /// Return the first non-empty attribute among `keys`.
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

    /// Deterministic preview of current attributes (useful for structured errors).
    pub fn attrs_preview(&self, max: usize) -> String {
        preview_attrs(&self.attrs, max)
    }
}

// ----------------------------
// AnnotationReader
// ----------------------------

/// Low-level streaming parser for GTF/GFF3 files.
///
/// Most users should not use this directly; prefer your higher-level builder.
///
/// Behavior:
/// - Skips blank lines
/// - Skips comment lines starting with '#'
/// - Tracks 1-based line numbers
pub struct AnnotationReader<R: BufRead> {
    reader: R,
    buf: String,
    line_no: usize,
}

impl<R: BufRead> AnnotationReader<R> {
    pub fn new(reader: R) -> Self {
        Self {
            reader,
            buf: String::new(),
            line_no: 0,
        }
    }

    /// Returns an iterator over parsed records.
    pub fn records(mut self) -> impl Iterator<Item = Result<AnnotationRecord, ParseError>> {
        std::iter::from_fn(move || loop {
            self.buf.clear();
            self.line_no += 1;

            match self.reader.read_line(&mut self.buf) {
                Ok(0) => return None,
                Ok(_) => {}
                Err(e) => {
                    return Some(Err(ParseError::IoPath {
                        path: "<reader>".to_string(),
                        source: e,
                    }))
                }
            }

            let line = self.buf.trim_end_matches(&['\n', '\r'][..]);
            if line.is_empty() || line.starts_with('#') {
                continue;
            }

            return Some(parse_record_line_with_lineno(line, self.line_no));
        })
    }
}

// ----------------------------
// Line parsing
// ----------------------------

/// Parse a single non-comment line into an `AnnotationRecord` with a line number.
pub fn parse_record_line_with_lineno(line: &str, line_no: usize) -> Result<AnnotationRecord, ParseError> {
    let line_preview = preview_line(line, 320);

    // GTF/GFF have 9 tab-separated columns:
    // seqname source feature start end score strand phase attributes
    let ncols = line.split('\t').count();
    if ncols != 9 {
        return Err(ParseError::MalformedLine {
            line_no: Some(line_no),
            expected: "GTF/GFF3: exactly 9 TAB-separated columns",
            problem: "wrong number of columns".to_string(),
            details: Some(format!("expected 9, found {}", ncols)),
            line_preview,
        });
    }

    let mut it = line.split('\t');

    let seqname = it.next().unwrap();
    let source = it.next().unwrap();
    let feature_type = it.next().unwrap();
    let start_s = it.next().unwrap();
    let end_s = it.next().unwrap();
    let score_s = it.next().unwrap();
    let strand_s = it.next().unwrap();
    let phase_s = it.next().unwrap();
    let attrs_s = it.next().unwrap();

    // Coordinates: input is 1-based inclusive; convert to 0-based half-open [start-1, end)
    let start_1: u64 = start_s.parse().map_err(|_| ParseError::BadCoordinates {
        line_no: Some(line_no),
        problem: format!("start '{}' is not an integer", start_s),
        line_preview: preview_line(line, 320),
    })?;

    let end_1: u64 = end_s.parse().map_err(|_| ParseError::BadCoordinates {
        line_no: Some(line_no),
        problem: format!("end '{}' is not an integer", end_s),
        line_preview: preview_line(line, 320),
    })?;

    if start_1 == 0 || end_1 == 0 || end_1 < start_1 {
        return Err(ParseError::BadCoordinates {
            line_no: Some(line_no),
            problem: format!("expected 1-based inclusive coordinates with end >= start; got {}..{}", start_1, end_1),
            line_preview,
        });
    }

    let start0 = (start_1 - 1) as u32;
    let end0 = end_1 as u32;

    let score = if score_s == "." {
        None
    } else {
        Some(score_s.parse::<f32>().map_err(|_| ParseError::MalformedLine {
            line_no: Some(line_no),
            expected: "score column to be '.' or a float",
            problem: format!("score '{}' is not a float", score_s),
            details: Some("score is column 6".to_string()),
            line_preview: preview_line(line, 320),
        })?)
    };

    let strand = match strand_s {
        "+" => Strand::Plus,
        "-" => Strand::Minus,
        "." | "?" => Strand::Unknown,
        _ => {
            return Err(ParseError::MalformedLine {
                line_no: Some(line_no),
                expected: "strand column to be one of: +, -, ., ?",
                problem: format!("invalid strand value '{}'", strand_s),
                details: Some("strand is column 7".to_string()),
                line_preview,
            })
        }
    };

    let phase = if phase_s == "." {
        None
    } else {
        let p: u8 = phase_s.parse().map_err(|_| ParseError::MalformedLine {
            line_no: Some(line_no),
            expected: "phase column to be '.' or 0/1/2",
            problem: format!("phase '{}' is not an integer", phase_s),
            details: Some("phase is column 8".to_string()),
            line_preview: preview_line(line, 320),
        })?;

        if p > 2 {
            return Err(ParseError::MalformedLine {
                line_no: Some(line_no),
                expected: "phase column to be '.' or 0/1/2",
                problem: format!("phase '{}' out of range", phase_s),
                details: Some("valid phase values are 0, 1, 2".to_string()),
                line_preview,
            });
        }
        Some(p)
    };

    let (dialect, attrs) = parse_attributes(attrs_s).map_err(|mut e| {
        // Attach line number + line preview if needed
        match &mut e {
            ParseError::MalformedLine { line_no: ln, line_preview: lp, .. } => {
                if ln.is_none() {
                    *ln = Some(line_no);
                }
                if lp.is_empty() {
                    *lp = preview_line(line, 320);
                }
            }
            ParseError::BadCoordinates { line_no: ln, line_preview: lp, .. } => {
                if ln.is_none() {
                    *ln = Some(line_no);
                }
                if lp.is_empty() {
                    *lp = preview_line(line, 320);
                }
            }
            _ => {}
        }
        e
    })?;

    Ok(AnnotationRecord {
        line_no,
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
        line_preview,
    })
}

/// Parse attributes field for either GFF3 or GTF.
///
/// Returns Result so it can fail *where the error occurs*.
pub fn parse_attributes(s: &str) -> Result<(Dialect, HashMap<String, String>), ParseError> {
    let raw = s.trim();

    // Dialect detection (cheap heuristics)
    let dialect = if raw.contains('=') {
        Dialect::Gff3
    } else if raw.contains('"') {
        Dialect::Gtf
    } else {
        Dialect::Unknown
    };

    // Strict GTF sanity: if quotes exist, they must be balanced.
    if dialect == Dialect::Gtf {
        let q = raw.chars().filter(|&c| c == '"').count();
        if q % 2 != 0 {
            return Err(ParseError::MalformedLine {
                line_no: None,
                expected: "GTF attributes with balanced double quotes",
                problem: "unbalanced double quote in attributes field".to_string(),
                details: Some(format!("found {} double-quote characters", q)),
                line_preview: preview_line(raw, 240),
            });
        }
    }

    let mut map = HashMap::new();

    match dialect {
        Dialect::Gff3 => {
            // key=value;key2=value2
            for part in raw.split(';') {
                let part = part.trim();
                if part.is_empty() {
                    continue;
                }
                if !part.contains('=') {
                    return Err(ParseError::MalformedLine {
                        line_no: None,
                        expected: "GFF3 attributes in the form key=value;key2=value2",
                        problem: "attribute segment missing '='".to_string(),
                        details: Some(format!("segment: '{}'", preview_line(part, 120))),
                        line_preview: preview_line(raw, 240),
                    });
                }
                let mut kv = part.splitn(2, '=');
                let k = kv.next().unwrap().trim();
                let v = kv.next().unwrap_or("").trim();
                if k.is_empty() {
                    return Err(ParseError::MalformedLine {
                        line_no: None,
                        expected: "GFF3 attribute keys to be non-empty",
                        problem: "empty attribute key".to_string(),
                        details: Some(format!("segment: '{}'", preview_line(part, 120))),
                        line_preview: preview_line(raw, 240),
                    });
                }
                let v = unquote(v);
                if !v.is_empty() {
                    map.insert(k.to_string(), v);
                }
            }
        }

        Dialect::Gtf => {
            // key "value"; key2 "value2";
            for part in raw.split(';') {
                let part = part.trim();
                if part.is_empty() {
                    continue;
                }
                let mut it = part.splitn(2, char::is_whitespace);
                let key = it.next().unwrap().trim();
                let rest = it.next().unwrap_or("").trim();

                if key.is_empty() {
                    return Err(ParseError::MalformedLine {
                        line_no: None,
                        expected: "GTF attributes in the form key \"value\";",
                        problem: "empty attribute key".to_string(),
                        details: Some(format!("segment: '{}'", preview_line(part, 120))),
                        line_preview: preview_line(raw, 240),
                    });
                }
                if rest.is_empty() {
                    return Err(ParseError::MalformedLine {
                        line_no: None,
                        expected: "GTF attributes in the form key \"value\"; (or key value;)",
                        problem: "attribute key without a value".to_string(),
                        details: Some(format!("key: '{}'", key)),
                        line_preview: preview_line(raw, 240),
                    });
                }

                let value = unquote(rest);
                if !value.is_empty() {
                    map.insert(key.to_string(), value);
                }
            }
        }

        Dialect::Unknown => {
            // best-effort
            for part in raw.split(';') {
                let part = part.trim();
                if part.is_empty() {
                    continue;
                }
                if part.contains('=') {
                    let mut kv = part.splitn(2, '=');
                    let k = kv.next().unwrap().trim();
                    let v = kv.next().unwrap_or("").trim();
                    if k.is_empty() {
                        return Err(ParseError::MalformedLine {
                            line_no: None,
                            expected: "attribute keys to be non-empty",
                            problem: "empty attribute key".to_string(),
                            details: Some(format!("segment: '{}'", preview_line(part, 120))),
                            line_preview: preview_line(raw, 240),
                        });
                    }
                    let v = unquote(v);
                    if !v.is_empty() {
                        map.insert(k.to_string(), v);
                    }
                } else {
                    let mut it = part.splitn(2, char::is_whitespace);
                    let key = it.next().unwrap().trim();
                    let rest = it.next().unwrap_or("").trim();
                    if key.is_empty() {
                        return Err(ParseError::MalformedLine {
                            line_no: None,
                            expected: "attribute keys to be non-empty",
                            problem: "empty attribute key".to_string(),
                            details: Some(format!("segment: '{}'", preview_line(part, 120))),
                            line_preview: preview_line(raw, 240),
                        });
                    }
                    if rest.is_empty() {
                        return Err(ParseError::MalformedLine {
                            line_no: None,
                            expected: "attribute segment to contain a value",
                            problem: "attribute key without a value".to_string(),
                            details: Some(format!("key: '{}'", key)),
                            line_preview: preview_line(raw, 240),
                        });
                    }
                    let value = unquote(rest);
                    if !value.is_empty() {
                        map.insert(key.to_string(), value);
                    }
                }
            }
        }
    }

    Ok((dialect, map))
}

// ----------------------------
// Helper functions
// ----------------------------

fn unquote(v: &str) -> String {
    let v = v.trim();
    let v = v.strip_prefix('"').unwrap_or(v);
    let v = v.strip_suffix('"').unwrap_or(v);
    v.to_string()
}

fn preview_line(line: &str, max: usize) -> String {
    let mut s = line.trim_end().to_string();
    if s.len() > max {
        s.truncate(max);
        s.push('…');
    }
    s
}

fn preview_attrs(attrs: &HashMap<String, String>, max: usize) -> String {
    let mut kv: Vec<_> = attrs.iter().collect();
    kv.sort_by(|a, b| a.0.cmp(b.0));

    let mut s = String::new();
    for (k, v) in kv {
        if !s.is_empty() {
            s.push_str("; ");
        }
        s.push_str(k);
        s.push('=');
        s.push_str(v);

        if s.len() > max {
            s.truncate(max);
            s.push('…');
            break;
        }
    }
    if s.is_empty() {
        "<no attributes>".to_string()
    } else {
        s
    }
}

// ----------------------------
// Tests
// ----------------------------

#[cfg(test)]
mod tests {
    use super::*;
    use std::io::Cursor;

    #[test]
    fn parse_gtf_line_ok() {
        let line = "chr1\tsrc\texon\t101\t150\t.\t+\t.\tgene_id \"G1\"; transcript_id \"T1\"; exon_number \"1\";";
        let rec = parse_record_line_with_lineno(line, 7).unwrap();

        assert_eq!(rec.dialect, Dialect::Gtf);
        assert_eq!(rec.line_no, 7);
        assert_eq!(rec.seqname, "chr1");
        assert_eq!(rec.feature_type, "exon");
        assert_eq!(rec.start0, 100);
        assert_eq!(rec.end0, 150);
        assert_eq!(rec.strand, Strand::Plus);

        assert_eq!(rec.attr("gene_id"), Some("G1"));
        assert_eq!(rec.attr("transcript_id"), Some("T1"));
        assert_eq!(rec.attr("exon_number"), Some("1"));
    }

    #[test]
    fn parse_gff3_line_ok() {
        let line = "chr2\tsrc\texon\t5\t20\t.\t-\t.\tID=ex1;Parent=tx1;gene_id=G9";
        let rec = parse_record_line_with_lineno(line, 1).unwrap();

        assert_eq!(rec.dialect, Dialect::Gff3);
        assert_eq!(rec.seqname, "chr2");
        assert_eq!(rec.feature_type, "exon");
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
        assert_eq!(recs[0].line_no, 2);
        assert_eq!(recs[1].line_no, 5);
        assert_eq!(recs[0].start0, 0);
        assert_eq!(recs[0].end0, 2);
        assert_eq!(recs[1].start0, 2);
        assert_eq!(recs[1].end0, 4);
    }

    #[test]
    fn parse_fails_on_wrong_column_count() {
        // 8 columns (missing attributes)
        let bad = "chr1\tsrc\texon\t101\t150\t.\t+\t.";
        let err = parse_record_line_with_lineno(bad, 42).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("wrong number of columns"));
        assert!(msg.contains("expected 9"));
        assert!(msg.contains("line 42"));
    }

    #[test]
    fn parse_fails_on_unbalanced_gtf_quotes() {
        let bad_attrs = "gene_id \"G1; transcript_id \"T1\";";
        let err = parse_attributes(bad_attrs).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("unbalanced double quote"));
    }

    #[test]
    fn parse_fails_on_invalid_strand() {
        let bad = "chr1\tsrc\texon\t101\t150\t.\tX\t.\tgene_id \"G1\"; transcript_id \"T1\";";
        let err = parse_record_line_with_lineno(bad, 9).unwrap_err();
        let msg = err.to_string();
        assert!(msg.contains("invalid strand"));
        assert!(msg.contains("line 9"));
    }
}

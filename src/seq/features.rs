//! Annotated features and the feature library (tab-separated name / sequence / type / colours).

use super::dna::{find_all, FindOpts};

pub type Rgb = [u8; 3];

/// The original's default colours when a feature has none: forward cyan, reverse green.
pub const DEFAULT_FWD: Rgb = [0, 255, 255];
pub const DEFAULT_REV: Rgb = [0, 255, 0];

/// Default (forward, reverse) colours per feature type, from the original's settings.
const TYPE_COLORS: &[(&str, &str, &str)] = &[
    ("V_segment", "cyan", "green"),
    ("repeat_unit", "#ff6666", "#ff6666"),
    ("promoter", "#346ee0", "#346ee0"),
    ("SNP", "#7e59d5", "#7e59d5"),
    ("misc_signal", "#008040", "#008040"),
    ("polyA_signal", "#ff3eee", "#ff3eee"),
    ("LTR", "#86f71d", "#d7336f"),
    ("primer_bind", "#14c0bd", "#4ec02b"),
    ("misc_recomb", "#8080ff", "#da6a3b"),
    ("repeat_region", "#0dfff7", "#0dfff7"),
    ("enhancer", "#5ac3fa", "#5ac3fa"),
    ("full_transcript", "#ff2074", "#ff2074"),
    ("coding_exon", "#ae1e5f", "#ae1e5f"),
    ("variation", "#ff9e39", "#ff9e39"),
    ("misc_binding", "#ff4809", "#ff4809"),
    ("rep_origin", "#999999", "#999999"),
    ("mRNA", "cyan", "cyan"),
    ("misc_feature", "#7eff74", "#7eff74"),
    ("gene", "#ff797d", "#ff102f"),
    ("RBS", "cyan", "cyan"),
    ("five_prime_UTR", "#d1d5ca", "#d1d5ca"),
    ("prim_transcript", "#ffe34b", "#ffe34b"),
    ("polyA_site", "cyan", "cyan"),
    ("exon", "#c83a90", "#c83a90"),
    ("CDS", "#e9d024", "#e9d024"),
    ("intron", "#eb606c", "#eb606c"),
    ("-35_signal", "cyan", "green"),
    ("UTR", "#526ba2", "#526ba2"),
    ("modified_base", "#92ff65", "#92ff65"),
    ("misc_structure", "#ff6462", "#ff6462"),
    ("deletion", "#1f5250", "#1f5250"),
    ("3'UTR", "#f76ca1", "#008040"),
    ("terminator", "#9d1b1c", "#9d1b1c"),
    ("substitution", "#ff6ff1", "#ff6ff1"),
    ("ncRNA_primary_transcript", "#ff29a8", "#5fd5a3"),
    ("three_prime_UTR", "#c5c0c6", "#c5c0c6"),
    ("5'UTR", "#ffd327", "#ffd327"),
    ("primer", "#00ffff", "#80ff00"),
];

/// Colour from `#rrggbb`, `#rgb`, `grayNN` or a common Tk colour name.
pub fn color_from_name(s: &str) -> Option<Rgb> {
    let s = s.trim().trim_matches('"').to_ascii_lowercase();
    if let Some(hex) = s.strip_prefix('#') {
        let digit = |i: usize, w: usize| u8::from_str_radix(&hex[i..i + w], 16).ok();
        return match hex.len() {
            6 => Some([digit(0, 2)?, digit(2, 2)?, digit(4, 2)?]),
            3 => Some([digit(0, 1)? * 17, digit(1, 1)? * 17, digit(2, 1)? * 17]),
            _ => None,
        };
    }
    if let Some(level) = s.strip_prefix("gray").or_else(|| s.strip_prefix("grey")) {
        let v = if level.is_empty() {
            190
        } else {
            (level.parse::<f32>().ok()? * 2.55).round() as u8
        };
        return Some([v, v, v]);
    }
    Some(match s.as_str() {
        "cyan" => [0, 255, 255],
        "green" => [0, 255, 0],
        "red" => [255, 0, 0],
        "blue" => [0, 0, 255],
        "yellow" => [255, 255, 0],
        "white" => [255, 255, 255],
        "black" => [0, 0, 0],
        "magenta" => [255, 0, 255],
        "orange" => [255, 165, 0],
        "purple" => [160, 32, 240],
        "pink" => [255, 192, 203],
        "brown" => [165, 42, 42],
        "gold" => [255, 215, 0],
        "navy" => [0, 0, 128],
        _ => return None,
    })
}

pub fn type_colors(kind: &str) -> (Rgb, Rgb) {
    TYPE_COLORS
        .iter()
        .find(|(k, _, _)| *k == kind)
        .map(|(_, f, r)| {
            (
                color_from_name(f).unwrap_or(DEFAULT_FWD),
                color_from_name(r).unwrap_or(DEFAULT_REV),
            )
        })
        .unwrap_or((DEFAULT_FWD, DEFAULT_REV))
}

#[derive(Clone, Debug, PartialEq)]
pub struct Feature {
    pub name: String,
    pub kind: String,
    /// 0-based half-open segments; a feature crossing the origin of a circle has two.
    pub spans: Vec<(usize, usize)>,
    pub reverse: bool,
    pub fwd_color: Rgb,
    pub rev_color: Rgb,
    pub qualifiers: Vec<(String, String)>,
}

impl Feature {
    pub fn new(
        name: impl Into<String>,
        kind: &str,
        start: usize,
        end: usize,
        reverse: bool,
    ) -> Feature {
        let (fwd_color, rev_color) = type_colors(kind);
        Feature {
            name: name.into(),
            kind: kind.to_string(),
            spans: vec![(start, end)],
            reverse,
            fwd_color,
            rev_color,
            qualifiers: Vec::new(),
        }
    }

    pub fn color(&self) -> Rgb {
        if self.reverse {
            self.rev_color
        } else {
            self.fwd_color
        }
    }

    pub fn start(&self) -> usize {
        self.spans.iter().map(|s| s.0).min().unwrap_or(0)
    }

    pub fn end(&self) -> usize {
        self.spans.iter().map(|s| s.1).max().unwrap_or(0)
    }

    pub fn len(&self) -> usize {
        self.spans.iter().map(|(s, e)| e - s).sum()
    }

    pub fn contains(&self, pos: usize) -> bool {
        self.spans.iter().any(|&(s, e)| pos >= s && pos < e)
    }

    /// Where the feature starts on its own strand (for labels and numbering), 1-based.
    pub fn location_text(&self) -> String {
        let parts: Vec<String> = self
            .spans
            .iter()
            .map(|&(s, e)| format!("{}..{}", s + 1, e))
            .collect();
        let joined = if parts.len() > 1 {
            format!("join({})", parts.join(","))
        } else {
            parts.join(",")
        };
        if self.reverse {
            format!("complement({joined})")
        } else {
            joined
        }
    }

    /// Keep the feature attached to its bases when `m` bases are inserted at `at`.
    pub fn after_insert(&mut self, at: usize, m: usize) {
        for (s, e) in &mut self.spans {
            if *s >= at {
                *s += m;
                *e += m;
            } else if *e > at {
                *e += m;
            }
        }
    }

    /// Shrink or move the feature after bases `[a, b)` were deleted. `false` if nothing is left.
    pub fn after_delete(&mut self, a: usize, b: usize) -> bool {
        let cut = b - a;
        let map = |p: usize| {
            if p >= b {
                p - cut
            } else if p > a {
                a
            } else {
                p
            }
        };
        self.spans = self
            .spans
            .iter()
            .map(|&(s, e)| (map(s), map(e)))
            .filter(|(s, e)| e > s)
            .collect();
        !self.spans.is_empty()
    }

    /// Move the feature when the sequence is rotated so base `at` becomes base 0.
    pub fn rotate(&mut self, at: usize, n: usize) {
        let mut out = Vec::new();
        for &(s, e) in &self.spans {
            if s >= at {
                out.push((s - at, e - at));
            } else if e <= at {
                out.push((s + n - at, e + n - at));
            } else {
                out.push((s + n - at, n));
                out.push((0, e - at));
            }
        }
        out.sort_unstable();
        self.spans = out;
    }

    /// Mirror the feature when the whole sequence is reverse-complemented.
    pub fn flip(&mut self, n: usize) {
        self.mirror(0, n);
    }

    /// Mirror the feature inside the range `[a, b)`, as when that stretch is reverse-complemented.
    /// Its direction flips (the original's `put_tags` with direction "reverse").
    pub fn mirror(&mut self, a: usize, b: usize) {
        for span in &mut self.spans {
            *span = (a + b - span.1, a + b - span.0);
        }
        self.spans.sort_unstable();
        self.reverse = !self.reverse;
    }

    /// The part of the feature inside `[a, b)` as a feature of its own; `None` if it misses the range.
    pub fn inside(&self, a: usize, b: usize) -> Option<Feature> {
        let spans: Vec<(usize, usize)> = self
            .spans
            .iter()
            .map(|&(s, e)| (s.max(a), e.min(b)))
            .filter(|(s, e)| e > s)
            .collect();
        (!spans.is_empty()).then(|| Feature {
            spans,
            ..self.clone()
        })
    }

    /// Drop the part inside `[a, b)`; `false` if nothing of the feature is left.
    pub fn remove_range(&mut self, a: usize, b: usize) -> bool {
        self.spans = self
            .spans
            .iter()
            .flat_map(|&(s, e)| [(s, e.min(a)), (s.max(b), e)])
            .filter(|(s, e)| e > s)
            .collect();
        !self.spans.is_empty()
    }

    /// The same feature `by` bases further along.
    pub fn shifted(mut self, by: usize) -> Feature {
        for (s, e) in &mut self.spans {
            *s += by;
            *e += by;
        }
        self
    }
}

// ---------------------------------------------------------------------------
// Library
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq)]
pub struct LibraryEntry {
    pub name: String,
    pub seq: Vec<u8>,
    pub kind: String,
    pub fwd_color: Rgb,
    pub rev_color: Rgb,
}

const BUILTIN_LIBRARY: &[(&str, &str, &str)] = &[
    ("T7 promoter", "TAATACGACTCACTATAGGG", "promoter"),
    ("T3 promoter", "AATTAACCCTCACTAAAGGG", "promoter"),
    ("SP6 promoter", "ATTTAGGTGACACTATAG", "promoter"),
    ("M13-fwd", "GTAAAACGACGGCCAGT", "primer_bind"),
    ("M13-rev", "CAGGAAACAGCTATGAC", "primer_bind"),
    ("lac operator", "AATTGTGAGCGGATAACAATT", "misc_binding"),
    (
        "T7 terminator",
        "CTAGCATAACCCCTTGGGGCCTCTAAACGGGTCTTGAGGGGTTTTTTG",
        "terminator",
    ),
    ("6xHis", "CAYCAYCAYCAYCAYCAY", "misc_feature"),
    ("FLAG", "GACTACAAGGACGACGATGACAAG", "misc_feature"),
    ("HA", "TACCCATACGATGTTCCAGATTACGCT", "misc_feature"),
    ("Myc", "GAACAAAAACTCATCTCAGAAGAGGATCTG", "misc_feature"),
    ("Strep-tag II", "TGGAGCCACCCGCAGTTCGAAAAA", "misc_feature"),
    (
        "V5",
        "GGTAAGCCTATCCCTAACCCTCTCCTCGGTCTCGATTCTACG",
        "misc_feature",
    ),
    ("TEV site", "GAAAACCTGTATTTTCAGGGC", "misc_feature"),
    ("Thrombin site", "CTGGTTCCGCGTGGATCC", "misc_feature"),
    ("SV40 NLS", "CCAAAAAAGAAGAGAAAGGTA", "misc_feature"),
    ("loxP", "ATAACTTCGTATAATGTATGCTATACGAAGTTAT", "misc_recomb"),
    ("FRT", "GAAGTTCCTATTCTCTAGAAAGTATAGGAACTTC", "misc_recomb"),
    ("attB1", "ACAAGTTTGTACAAAAAAGCAGGCT", "misc_recomb"),
    ("attB2", "ACCCAGCTTTCTTGTACAAAGTGGT", "misc_recomb"),
    ("Kozak", "GCCACCATGG", "misc_signal"),
];

pub fn builtin_library() -> Vec<LibraryEntry> {
    BUILTIN_LIBRARY
        .iter()
        .map(|&(name, seq, kind)| {
            let (fwd_color, rev_color) = type_colors(kind);
            LibraryEntry {
                name: name.into(),
                seq: seq.as_bytes().to_vec(),
                kind: kind.into(),
                fwd_color,
                rev_color,
            }
        })
        .collect()
}

/// Parse a library file in the original's tab-separated format.
/// Entries with gap markers (`#`, `+`), protein searches or fewer than four bases are skipped.
pub fn parse_library(text: &str) -> Vec<LibraryEntry> {
    text.lines()
        .filter(|l| !l.starts_with('#') && !l.trim().is_empty())
        .filter_map(|line| {
            let cols: Vec<&str> = line.split('\t').collect();
            if cols.len() < 2
                || cols.len() > 8
                || cols.get(6).is_some_and(|c| matches!(c.trim(), "1" | "2"))
            {
                return None;
            }
            let core = cols[1].split_once('<').map_or(cols[1], |(_, after)| after);
            let core = core.split_once('>').map_or(core, |(before, _)| before);
            let seq: Vec<u8> = core
                .bytes()
                .filter(|c| !c.is_ascii_whitespace())
                .map(|c| c.to_ascii_uppercase())
                .collect();
            let valid = seq.len() > 3 && seq.iter().all(|&c| super::dna::iupac_mask(c) != 0);
            if !valid {
                return None;
            }
            let kind = cols
                .get(2)
                .map(|k| k.trim())
                .filter(|k| !k.is_empty())
                .unwrap_or("misc_feature");
            let (default_f, default_r) = type_colors(kind);
            let color = |i: usize, fallback: Rgb| {
                cols.get(i)
                    .and_then(|c| color_from_name(c))
                    .unwrap_or(fallback)
            };
            Some(LibraryEntry {
                name: cols[0].trim().to_string(),
                seq,
                kind: kind.to_string(),
                fwd_color: color(3, default_f),
                rev_color: color(4, default_r),
            })
        })
        .collect()
}

/// Features found by matching library sequences on both strands (exact, IUPAC aware).
pub fn annotate(seq: &[u8], circular: bool, library: &[LibraryEntry]) -> Vec<Feature> {
    let n = seq.len();
    let opts = FindOpts {
        case_sensitive: false,
        both_strands: true,
        max_mismatch: 0,
        circular,
    };
    let mut found: Vec<Feature> = Vec::new();
    for entry in library {
        for hit in find_all(seq, &entry.seq, &opts) {
            let end = hit.start + hit.len;
            let spans = if end > n {
                vec![(0, end - n), (hit.start, n)]
            } else {
                vec![(hit.start, end)]
            };
            let f = Feature {
                name: entry.name.clone(),
                kind: entry.kind.clone(),
                spans,
                reverse: hit.reverse,
                fwd_color: entry.fwd_color,
                rev_color: entry.rev_color,
                qualifiers: Vec::new(),
            };
            if !found.iter().any(|g| g.name == f.name && g.spans == f.spans) {
                found.push(f);
            }
        }
    }
    found.sort_by_key(|f| (f.start(), f.name.clone()));
    found
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn colors_parse_hex_names_and_grays() {
        assert_eq!(color_from_name("#346ee0"), Some([0x34, 0x6e, 0xe0]));
        assert_eq!(color_from_name("#f0a"), Some([255, 0, 170]));
        assert_eq!(color_from_name("cyan"), Some([0, 255, 255]));
        assert_eq!(color_from_name("gray50"), Some([128, 128, 128]));
        assert_eq!(color_from_name("nonsense"), None);
        assert_eq!(type_colors("promoter").0, [0x34, 0x6e, 0xe0]);
        assert_eq!(type_colors("unknown_type"), (DEFAULT_FWD, DEFAULT_REV));
    }

    #[test]
    fn insertions_extend_or_shift_features() {
        let mut f = Feature::new("x", "gene", 10, 20, false);
        f.after_insert(5, 3);
        assert_eq!(f.spans, vec![(13, 23)]);
        f.after_insert(15, 2);
        assert_eq!(f.spans, vec![(13, 25)]);
        f.after_insert(25, 2);
        assert_eq!(f.spans, vec![(13, 25)]);
    }

    #[test]
    fn deletions_shrink_move_or_remove_features() {
        let mut f = Feature::new("x", "gene", 10, 20, false);
        assert!(f.after_delete(0, 4));
        assert_eq!(f.spans, vec![(6, 16)]);
        assert!(f.after_delete(8, 12));
        assert_eq!(f.spans, vec![(6, 12)]);
        assert!(!f.after_delete(0, 30));
    }

    #[test]
    fn rotating_splits_features_that_cross_the_new_origin() {
        let mut f = Feature::new("x", "gene", 8, 14, false);
        f.rotate(10, 20);
        assert_eq!(f.spans, vec![(0, 4), (18, 20)]);
        let mut g = Feature::new("y", "gene", 2, 5, true);
        g.rotate(10, 20);
        assert_eq!(g.spans, vec![(12, 15)]);
    }

    #[test]
    fn flipping_mirrors_coordinates_and_strand() {
        let mut f = Feature::new("x", "CDS", 2, 7, false);
        f.flip(20);
        assert_eq!((f.spans.clone(), f.reverse), (vec![(13, 18)], true));
        assert_eq!(f.location_text(), "complement(14..18)");
    }

    #[test]
    fn library_files_parse_like_the_original() {
        let text =
            "# comment\nT7\tgtgaattg<TAATACGACTCACTATAGGG\tprimer_bind\tcyan\tgreen\t\t0\t\n\
                    TrfA\tMNRTFDRKAYRQ\tCDS\tred\tblue\t\t1\t\nShort\tACG\t\nBad\tAC#GT\tgene\n\
                    Promoter\tTTGACAGCTTATCATCG\tpromoter\t#346ee0\t#346ee0\t\t0\t";
        let lib = parse_library(text);
        assert_eq!(
            lib.iter().map(|e| e.name.as_str()).collect::<Vec<_>>(),
            vec!["T7", "Promoter"]
        );
        assert_eq!(lib[0].seq, b"TAATACGACTCACTATAGGG");
        assert_eq!(lib[0].rev_color, [0, 255, 0]);
        assert_eq!(lib[1].fwd_color, [0x34, 0x6e, 0xe0]);
    }

    #[test]
    fn annotate_finds_features_on_both_strands_and_across_the_origin() {
        let lib = builtin_library();
        let t7 = b"TAATACGACTCACTATAGGG";
        let mut seq = b"CCCCCCCCCC".to_vec();
        seq.extend_from_slice(t7);
        seq.extend_from_slice(&super::super::dna::revcom(b"GACTACAAGGACGACGATGACAAG"));
        seq.extend_from_slice(b"CCCCCCCCCC");
        let found = annotate(&seq, false, &lib);
        let names: Vec<_> = found.iter().map(|f| (f.name.as_str(), f.reverse)).collect();
        assert_eq!(names, vec![("T7 promoter", false), ("FLAG", true)]);

        let mut wrapped = t7[10..].to_vec();
        wrapped.extend_from_slice(&[b'C'; 30]);
        wrapped.extend_from_slice(&t7[..10]);
        let found = annotate(&wrapped, true, &lib);
        assert_eq!(found.len(), 1);
        assert_eq!(found[0].spans, vec![(0, 10), (40, 50)]);
        assert!(annotate(&wrapped, false, &lib).is_empty());
    }
}

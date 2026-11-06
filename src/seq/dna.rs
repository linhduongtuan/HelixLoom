//! Sequence maths shared by every tool: complements, GC, melting temperature,
//! translation, ORF finding and pattern search.
//! Formulas follow the original program's procedures (revcom, translate, Tm, GC%, orf_list).

pub const STANDARD_CODE: &str = "FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG";

pub struct GeneticCode {
    pub id: u8,
    pub name: &'static str,
    /// 64 one-letter amino acids, codons ordered T,C,A,G in each position.
    pub table: &'static str,
}

pub const GENETIC_CODES: &[GeneticCode] = &[
    GeneticCode {
        id: 1,
        name: "Standard",
        table: STANDARD_CODE,
    },
    GeneticCode {
        id: 2,
        name: "Vertebrate mitochondrial",
        table: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
    },
    GeneticCode {
        id: 3,
        name: "Yeast mitochondrial",
        table: "FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    },
    GeneticCode {
        id: 4,
        name: "Mold / protozoan mitochondrial",
        table: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    },
    GeneticCode {
        id: 5,
        name: "Invertebrate mitochondrial",
        table: "FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
    },
    GeneticCode {
        id: 11,
        name: "Bacterial / plant plastid",
        table: STANDARD_CODE,
    },
    GeneticCode {
        id: 12,
        name: "Alternative yeast nuclear",
        table: "FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
    },
];

// ---------------------------------------------------------------------------
// IUPAC helpers
// ---------------------------------------------------------------------------

/// Bit mask of the bases an IUPAC letter stands for (A=1, C=2, G=4, T=8); 0 for non-bases.
pub fn iupac_mask(c: u8) -> u8 {
    match c.to_ascii_uppercase() {
        b'A' => 1,
        b'C' => 2,
        b'G' => 4,
        b'T' | b'U' => 8,
        b'R' => 5,
        b'Y' => 10,
        b'S' => 6,
        b'W' => 9,
        b'K' => 12,
        b'M' => 3,
        b'B' => 14,
        b'D' => 13,
        b'H' => 11,
        b'V' => 7,
        b'N' => 15,
        _ => 0,
    }
}

/// IUPAC letter for a set of bases.
pub fn iupac_from_set(bases: &[u8]) -> u8 {
    const BY_MASK: &[u8; 16] = b"-ACMGRSVTWYHKDBN";
    BY_MASK[bases.iter().fold(0u8, |m, &b| m | iupac_mask(b)) as usize]
}

/// Whether sequence base `s` is allowed by pattern letter `p`.
pub fn base_matches(p: u8, s: u8, case_sensitive: bool) -> bool {
    if case_sensitive && p.is_ascii_lowercase() != s.is_ascii_lowercase() {
        return false;
    }
    let (mp, ms) = (iupac_mask(p), iupac_mask(s));
    if mp == 0 || ms == 0 {
        return p.eq_ignore_ascii_case(&s);
    }
    ms & !mp == 0
}

/// The original's complement map (`com` and `revcom`, AppMain.tcl lines 14451 and 14468): the IUPAC letters
/// pair up, case is kept, and everything else, U included, stays as it is.
pub fn complement(c: u8) -> u8 {
    let up = match c.to_ascii_uppercase() {
        b'A' => b'T',
        b'T' => b'A',
        b'C' => b'G',
        b'G' => b'C',
        b'R' => b'Y',
        b'Y' => b'R',
        b'K' => b'M',
        b'M' => b'K',
        b'B' => b'V',
        b'V' => b'B',
        b'D' => b'H',
        b'H' => b'D',
        other => other,
    };
    if c.is_ascii_lowercase() {
        up.to_ascii_lowercase()
    } else {
        up
    }
}

pub fn revcom(seq: &[u8]) -> Vec<u8> {
    seq.iter().rev().map(|&c| complement(c)).collect()
}

// ---------------------------------------------------------------------------
// Composition and melting temperature
// ---------------------------------------------------------------------------

/// Percentage of bases that are not A/T/U or an ambiguity code containing only A/T.
pub fn gc_percent(seq: &[u8]) -> f64 {
    if seq.is_empty() {
        return 0.0;
    }
    let gc = seq
        .iter()
        .filter(|c| !b"ATUNBDHKMRVWY*".contains(&c.to_ascii_uppercase()))
        .count();
    100.0 * gc as f64 / seq.len() as f64
}

/// Extinction coefficient (L/mmol/cm), nearest-neighbour sum from the original's `extinction`.
/// `None` for an empty sequence or any base other than A, C, G, T (the original fails there too).
pub fn extinction(seq: &[u8]) -> Option<f64> {
    let single = |c: u8| match c.to_ascii_uppercase() {
        b'A' => Some(15.4),
        b'C' => Some(7.4),
        b'G' => Some(11.5),
        b'T' => Some(8.7),
        _ => None,
    };
    let pair = |a: u8, b: u8| match (a.to_ascii_uppercase(), b.to_ascii_uppercase()) {
        (b'A', b'A') => Some(13.7),
        (b'A', b'C') => Some(10.6),
        (b'A', b'G') => Some(12.5),
        (b'A', b'T') => Some(11.4),
        (b'C', b'A') => Some(10.6),
        (b'C', b'C') => Some(7.3),
        (b'C', b'G') => Some(9.0),
        (b'C', b'T') => Some(7.6),
        (b'G', b'A') => Some(12.6),
        (b'G', b'C') => Some(8.8),
        (b'G', b'G') => Some(10.8),
        (b'G', b'T') => Some(10.0),
        (b'T', b'A') => Some(11.7),
        (b'T', b'C') => Some(8.1),
        (b'T', b'G') => Some(9.5),
        (b'T', b'T') => Some(8.4),
        _ => None,
    };
    let mut e = single(*seq.first()?)?;
    for w in seq.windows(2) {
        e = (e + 2.0 * pair(w[0], w[1])?) - single(w[0])?;
    }
    Some(e)
}

/// Upper <-> lower for the letters the original's `swap` lists (the IUPAC codes); others stay.
pub fn swap_case(c: u8) -> u8 {
    if b"ABCDGHKMNRSTVWYabcdghkmnrstvwy".contains(&c) {
        c ^ 0x20
    } else {
        c
    }
}

/// (deltaH, deltaS) x10 for a dinucleotide, SantaLucia 1998 nearest-neighbour values.
pub(crate) fn nearest_neighbour(a: u8, b: u8) -> (i32, i32) {
    match (a, b) {
        (b'a', b'a') | (b't', b't') => (79, 222),
        (b'a', b't') => (72, 204),
        (b't', b'a') => (72, 213),
        (b'c', b'a') | (b't', b'g') => (85, 227),
        (b'g', b't') | (b'a', b'c') => (84, 224),
        (b'c', b't') | (b'a', b'g') => (78, 210),
        (b'g', b'a') | (b't', b'c') => (82, 222),
        (b'c', b'g') => (106, 272),
        (b'g', b'c') => (98, 244),
        _ => (80, 199),
    }
}

/// Melting temperature in degrees C, rounded to 0.1 and never below 0; `None` when the
/// sequence is empty or longer than 100. `dna_nm` and `salt_mm` default to 250 nM and 50 mM.
pub fn melting_temp(seq: &[u8], dna_nm: f64, salt_mm: f64) -> Option<f64> {
    melting_temp_raw(seq, dna_nm, salt_mm).map(|tm| (tm.max(0.0) * 10.0).round() / 10.0)
}

/// The same value before rounding and clamping (used by the primer filters).
pub fn melting_temp_raw(seq: &[u8], dna_nm: f64, salt_mm: f64) -> Option<f64> {
    let s: Vec<u8> = seq
        .iter()
        .filter(|c| !matches!(c, b'\n' | b'\r'))
        .map(|c| match c.to_ascii_lowercase() {
            b @ (b'a' | b'c' | b'g' | b't') => b,
            _ => b't',
        })
        .collect();
    if s.is_empty() || s.len() > 100 {
        return None;
    }
    let (mut h, mut e) = (0i32, 0i32);
    for end in [s[0], s[s.len() - 1]] {
        if end == b'a' || end == b't' {
            h -= 23;
            e -= 41;
        } else {
            h -= 1;
            e += 28;
        }
    }
    for w in s.windows(2) {
        let (dh, ds) = nearest_neighbour(w[0], w[1]);
        h += dh;
        e += ds;
    }
    let log_dna = 1.987 * (dna_nm / 1e9).ln();
    let log_salt = 12.0 * (salt_mm / 1000.0).log10();
    Some(h as f64 * 100.0 / (e as f64 * 0.1 - log_dna) - 273.15 + log_salt)
}

// ---------------------------------------------------------------------------
// Translation
// ---------------------------------------------------------------------------

fn codon_index(c: u8) -> Option<usize> {
    match c.to_ascii_uppercase() {
        b'T' | b'U' => Some(0),
        b'C' => Some(1),
        b'A' => Some(2),
        b'G' => Some(3),
        _ => None,
    }
}

/// One-letter translation; codons with ambiguous bases give `X`, a trailing partial codon is dropped.
pub fn translate(dna: &[u8], table: &str) -> Vec<u8> {
    let t = table.as_bytes();
    dna.as_chunks::<3>()
        .0
        .iter()
        .map(
            |c| match (codon_index(c[0]), codon_index(c[1]), codon_index(c[2])) {
                (Some(a), Some(b), Some(d)) => t[16 * a + 4 * b + d],
                _ => b'X',
            },
        )
        .collect()
}

pub fn three_letter(aa: u8) -> &'static str {
    match aa.to_ascii_uppercase() {
        b'F' => "Phe",
        b'L' => "Leu",
        b'S' => "Ser",
        b'Y' => "Tyr",
        b'C' => "Cys",
        b'*' => "End",
        b'W' => "Trp",
        b'P' => "Pro",
        b'H' => "His",
        b'Q' => "Gln",
        b'R' => "Arg",
        b'I' => "Ile",
        b'M' => "Met",
        b'T' => "Thr",
        b'N' => "Asn",
        b'K' => "Lys",
        b'V' => "Val",
        b'A' => "Ala",
        b'D' => "Asp",
        b'E' => "Glu",
        b'G' => "Gly",
        _ => "Xxx",
    }
}

fn is_stop(codon: &[u8]) -> bool {
    matches!(
        codon,
        [b'T', b'A', b'A'] | [b'T', b'A', b'G'] | [b'T', b'G', b'A']
    )
}

// ---------------------------------------------------------------------------
// ORFs
// ---------------------------------------------------------------------------

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum StartRule {
    /// ORF begins at ATG.
    Met,
    /// ORF begins right after any stop codon.
    AfterStop,
}

#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub enum Strand {
    Top,
    Bottom,
    Both,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Orf {
    /// Leftmost base on the top strand (0-based); wraps past the end on circular sequences.
    pub start: usize,
    pub length: usize,
    pub reverse: bool,
    /// Reading frame 0..=2 in the strand's own coordinates.
    pub frame: u8,
    /// No stop codon before the end of a linear sequence.
    pub open_end: bool,
}

impl Orf {
    /// Exclusive end on the top strand (may be smaller than `start` when wrapping).
    pub fn end(&self, n: usize) -> usize {
        (self.start + self.length) % n.max(1)
    }
}

/// ORFs longer than `min_bp` (the original tests `length - 1 > min_bp`).
/// With `nested` every start is reported; otherwise only the longest ORF per stop codon.
pub fn find_orfs(
    seq: &[u8],
    circular: bool,
    strand: Strand,
    min_bp: usize,
    rule: StartRule,
    nested: bool,
) -> Vec<Orf> {
    let n = seq.len();
    if n < 3 {
        return Vec::new();
    }
    let top: Vec<u8> = seq.to_ascii_uppercase();
    let mut out = Vec::new();
    if strand != Strand::Bottom {
        orfs_on_strand(&top, circular, min_bp, rule, nested, false, &mut out);
    }
    if strand != Strand::Top {
        orfs_on_strand(
            &revcom(&top),
            circular,
            min_bp,
            rule,
            nested,
            true,
            &mut out,
        );
    }
    out.sort_by_key(|o| (o.start, o.reverse));
    out
}

fn orfs_on_strand(
    text: &[u8],
    circular: bool,
    min_bp: usize,
    rule: StartRule,
    nested: bool,
    reverse: bool,
    out: &mut Vec<Orf>,
) {
    let n = text.len();
    let doubled: Vec<u8> = if circular {
        [text, text].concat()
    } else {
        text.to_vec()
    };
    let starts: Vec<usize> = match rule {
        StartRule::Met => (0..n)
            .filter(|&b| doubled.get(b..b + 3) == Some(b"ATG"))
            .collect(),
        StartRule::AfterStop => (0..n)
            .filter(|&p| doubled.get(p..p + 3).is_some_and(is_stop))
            .map(|p| p + 3)
            .filter(|&b| b < n)
            .collect(),
    };
    let mut used_stops = std::collections::HashSet::new();
    for b in starts {
        let stop = (b..doubled.len().saturating_sub(2))
            .step_by(3)
            .find(|&c| is_stop(&doubled[c..c + 3]));
        let (length, open_end) = match stop {
            Some(c) => (c + 3 - b, false),
            None if !circular => ((doubled.len() - b) / 3 * 3, true),
            None => continue,
        };
        if length < 3 || (length as isize - 1) <= min_bp as isize {
            continue;
        }
        if !nested && stop.is_some_and(|c| !used_stops.insert(c)) {
            continue;
        }
        let start = if reverse {
            n - 1 - ((b + length - 1) % n)
        } else {
            b
        };
        out.push(Orf {
            start,
            length,
            reverse,
            frame: (b % 3) as u8,
            open_end,
        });
    }
}

/// The original's selection summary, e.g. `M>85*/<85`: forward frame / reverse frame.
pub fn orf_summary(selection: &[u8]) -> String {
    let up = selection.to_ascii_uppercase();
    let n = up.len();
    let mut s = String::new();
    if up.starts_with(b"ATG") {
        s.push('M');
    }
    match (0..n.saturating_sub(2))
        .step_by(3)
        .find(|&i| is_stop(&up[i..i + 3]))
    {
        Some(i) if n - (i + 2) < 4 => s.push_str(&format!(">{}*", i / 3)),
        Some(_) => s.push('-'),
        None => s.push_str(&format!(">{}", n / 3)),
    }
    s.push('/');
    let reverse_stops: Vec<usize> = (0..n.saturating_sub(5))
        .filter(|&p| (n - p).is_multiple_of(3) && matches!(&up[p..p + 3], b"TCA" | b"CTA" | b"TTA"))
        .collect();
    if reverse_stops.iter().any(|&p| p >= 3) {
        s.push('-');
    } else if !reverse_stops.is_empty() {
        s.push_str(&format!("*{}<", (n / 3).saturating_sub(1)));
    } else {
        s.push_str(&format!("<{}", n / 3));
    }
    if up.ends_with(b"CAT") {
        s.push('M');
    }
    s
}

// ---------------------------------------------------------------------------
// Pattern search
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Hit {
    pub start: usize,
    pub len: usize,
    pub reverse: bool,
}

#[derive(Clone, Debug)]
pub struct FindOpts {
    pub case_sensitive: bool,
    pub both_strands: bool,
    pub max_mismatch: usize,
    pub circular: bool,
}

impl Default for FindOpts {
    fn default() -> Self {
        FindOpts {
            case_sensitive: false,
            both_strands: true,
            max_mismatch: 0,
            circular: false,
        }
    }
}

/// All places where `pattern` (IUPAC letters allowed) occurs, on both strands if asked.
pub fn find_all(seq: &[u8], pattern: &[u8], o: &FindOpts) -> Vec<Hit> {
    let pat: Vec<u8> = pattern
        .iter()
        .copied()
        .filter(|c| !c.is_ascii_whitespace())
        .collect();
    if pat.is_empty() || seq.is_empty() {
        return Vec::new();
    }
    let mut hits = scan(seq, &pat, o, false);
    let rc = revcom(&pat);
    if o.both_strands && !(rc.eq_ignore_ascii_case(&pat) && o.max_mismatch == 0) {
        hits.extend(scan(seq, &rc, o, true));
    }
    hits.sort_by_key(|h| (h.start, h.reverse));
    hits
}

fn scan(seq: &[u8], pat: &[u8], o: &FindOpts, reverse: bool) -> Vec<Hit> {
    let n = seq.len();
    let text: Vec<u8> = if o.circular && pat.len() > 1 {
        seq.iter()
            .chain(seq.iter().take(pat.len() - 1))
            .copied()
            .collect()
    } else {
        seq.to_vec()
    };
    if text.len() < pat.len() {
        return Vec::new();
    }
    (0..=text.len() - pat.len())
        .filter(|&i| i < n)
        .filter(|&i| {
            let bad = pat
                .iter()
                .zip(&text[i..])
                .filter(|(&p, &s)| !base_matches(p, s, o.case_sensitive))
                .take(o.max_mismatch + 1)
                .count();
            bad <= o.max_mismatch
        })
        .map(|start| Hit {
            start,
            len: pat.len(),
            reverse,
        })
        .collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tm(s: &str) -> Option<f64> {
        melting_temp(s.as_bytes(), 250.0, 50.0)
    }

    // Reference values come from running the original's Tm procedure.
    #[test]
    fn melting_temperature_matches_the_original() {
        assert_eq!(tm("ATGGCGTTTGAACGGCTCAA"), Some(59.3));
        assert_eq!(tm("GATCAAGGAA"), Some(21.6));
        assert_eq!(tm("GCGCGCGCGCGCGCGCGCGC"), Some(79.5));
        assert_eq!(tm("TAATACGACTCACTATAGGG"), Some(48.3));
        assert_eq!(tm("GGGGGGGGGGGGGGGGGGGGGG"), Some(81.1));
        assert_eq!(tm("AAAAAAAAAAAAAAAAAAAAAAAAA"), Some(44.7));
        assert_eq!(tm("CAGGAAACAGCTATGAC"), Some(47.3));
        assert_eq!(tm("GTAAAACGACGGCCAGT"), Some(53.0));
        assert_eq!(tm("ATGNNNCCC"), Some(18.5));
        assert_eq!(tm("ACGT"), Some(0.0));
        assert_eq!(tm("A"), Some(0.0));
        assert_eq!(
            melting_temp(b"ATGGCGTTTGAACGGCTCAA", 500.0, 100.0),
            Some(64.0)
        );
        assert_eq!(tm(&"ACGT".repeat(30)), None);
    }

    #[test]
    fn gc_percent_matches_the_original() {
        assert_eq!(gc_percent(b"ATGGCGTTTGAACGGCTCAA"), 50.0);
        assert!((gc_percent(b"CAGGAAACAGCTATGAC") - 47.05882352941177).abs() < 1e-9);
        assert!((gc_percent(b"ATGNNNCCC") - 400.0 / 9.0).abs() < 1e-9);
        assert_eq!(gc_percent(b""), 0.0);
    }

    // Reference values come from running the original's extinction procedure.
    #[test]
    fn extinction_matches_the_original() {
        let near = |s: &[u8], want: f64| {
            let got = extinction(s).expect("a valid DNA sequence");
            assert!((got - want).abs() < 1e-9, "{got} vs {want}");
        };
        near(b"ATGGCGTTTGAACGGCTCAA", 194.29999999999993);
        near(b"GATCAAGGAA", 111.29999999999998);
        near(b"A", 15.4);
        near(b"AC", 21.200000000000003);
        near(b"GCGCGC", 51.00000000000001);
        near(b"CAGGAAACAGCTATGAC", 173.59999999999997);
        near(b"acgt", 40.300000000000004);
        // the original stops with an error on these
        assert_eq!(extinction(b""), None);
        assert_eq!(extinction(b"ACGN"), None);
    }

    #[test]
    fn swap_case_only_touches_the_letters_the_original_lists() {
        let swapped: Vec<u8> = b"AcGtNuEx".iter().map(|&c| swap_case(c)).collect();
        assert_eq!(swapped, b"aCgTnuEx");
    }

    #[test]
    fn reverse_complement_handles_case_and_ambiguity() {
        assert_eq!(revcom(b"GAATTC"), b"GAATTC");
        assert_eq!(revcom(b"AAGCtt"), b"aaGCTT");
        assert_eq!(revcom(b"RYKMBVDHN"), b"NDHBVKMRY".to_vec());
    }

    #[test]
    fn iupac_sets_round_trip() {
        assert_eq!(iupac_from_set(b"AG"), b'R');
        assert_eq!(iupac_from_set(b"CT"), b'Y');
        assert_eq!(iupac_from_set(b"ACGT"), b'N');
        assert!(base_matches(b'R', b'A', false));
        assert!(!base_matches(b'R', b'C', false));
        assert!(!base_matches(b'A', b'N', false));
    }

    #[test]
    fn translation_uses_the_code_and_flags_ambiguous_codons() {
        assert_eq!(translate(b"ATGGCCTAAGG", STANDARD_CODE), b"MA*");
        assert_eq!(translate(b"ATGNCCTGA", STANDARD_CODE), b"MX*");
        assert_eq!(translate(b"TGA", GENETIC_CODES[1].table), b"W");
        assert_eq!(three_letter(b'*'), "End");
    }

    #[test]
    fn finds_orfs_on_both_strands() {
        let mut s = b"CCCATGAAATTTGGGCCCTAACCC".to_vec();
        let top = find_orfs(&s, false, Strand::Top, 10, StartRule::Met, false);
        assert_eq!(top.len(), 1);
        assert_eq!((top[0].start, top[0].length), (3, 18));
        s = revcom(&s);
        let bottom = find_orfs(&s, false, Strand::Both, 10, StartRule::Met, false);
        assert_eq!(bottom.len(), 1);
        assert!(bottom[0].reverse);
        assert_eq!((bottom[0].start, bottom[0].length), (3, 18));
    }

    #[test]
    fn orf_can_wrap_around_the_origin_of_a_circle() {
        // ATG at 24, then AAA CCC, and the TAA stop sits just past the origin (bases 3..6).
        let s = b"CCCTAAGGGGGGGGGGGGGGGGGGATGAAA".to_vec();
        let wrapped = find_orfs(&s, true, Strand::Top, 0, StartRule::Met, false);
        assert_eq!(wrapped.len(), 1);
        assert_eq!(
            (wrapped[0].start, wrapped[0].length, wrapped[0].end(s.len())),
            (24, 12, 6)
        );
        assert!(!wrapped[0].open_end);
        let linear = find_orfs(&s, false, Strand::Top, 0, StartRule::Met, false);
        assert_eq!(
            (linear[0].start, linear[0].length, linear[0].open_end),
            (24, 6, true)
        );
    }

    #[test]
    fn minimum_length_follows_the_original_rule() {
        let s = b"ATGAAATAA".to_vec();
        assert_eq!(
            find_orfs(&s, false, Strand::Top, 7, StartRule::Met, false).len(),
            1
        );
        assert!(find_orfs(&s, false, Strand::Top, 8, StartRule::Met, false).is_empty());
    }

    #[test]
    fn nested_starts_are_optional() {
        let s = b"ATGATGAAATAA".to_vec();
        assert_eq!(
            find_orfs(&s, false, Strand::Top, 0, StartRule::Met, true).len(),
            2
        );
        assert_eq!(
            find_orfs(&s, false, Strand::Top, 0, StartRule::Met, false).len(),
            1
        );
    }

    #[test]
    fn selection_summary_describes_both_frames() {
        assert_eq!(orf_summary(b"ATGAAATAA"), "M>2*/<3");
        assert_eq!(orf_summary(b"GGGAAACCC"), ">3/<3");
        assert_eq!(orf_summary(b"TTAGGGAAACAT"), ">4/*3<M");
    }

    #[test]
    fn finds_degenerate_patterns_on_both_strands() {
        let s = b"AAGAATTCAAGGATCCTT".to_vec();
        let hits = find_all(&s, b"GAATTC", &FindOpts::default());
        assert_eq!(
            hits,
            vec![Hit {
                start: 2,
                len: 6,
                reverse: false
            }]
        );
        let hits = find_all(&s, b"GGNNCC", &FindOpts::default());
        assert_eq!(
            hits,
            vec![Hit {
                start: 10,
                len: 6,
                reverse: false
            }]
        );
        let rc = find_all(b"GGATCC", b"GGTTCC", &FindOpts::default());
        assert_eq!(rc, vec![]);
        // One mismatch allowed: GGTTCC matches GGATCC directly and its reverse complement GGAACC too.
        let mm = find_all(
            b"GGATCC",
            b"GGTTCC",
            &FindOpts {
                max_mismatch: 1,
                ..Default::default()
            },
        );
        assert_eq!(
            mm.iter().map(|h| h.reverse).collect::<Vec<_>>(),
            vec![false, true]
        );
    }

    #[test]
    fn circular_search_finds_sites_across_the_origin() {
        let s = b"ATTCAAAAAAAAGA".to_vec();
        let linear = find_all(&s, b"GAATTC", &FindOpts::default());
        assert!(linear.is_empty());
        let circ = find_all(
            &s,
            b"GAATTC",
            &FindOpts {
                circular: true,
                ..Default::default()
            },
        );
        assert_eq!(
            circ,
            vec![Hit {
                start: 12,
                len: 6,
                reverse: false
            }]
        );
    }
}

//! PCR as the original's PCR dialog does it: where a primer anneals (`pcr_dialog_scan` @0x297522 with
//! `search_with_mm3` @0x2a43a7 and `pcr_dialog_precompute` @0x29863b) and the product two annealed primers
//! make (`pcr_dialog_do_pcr` @0x2944ee) with its warnings, its file comment and its two primer features.
//! The rules, quirks included, were checked against the original's own Tcl run on random templates and
//! primers (`pcr/golden.txt`, see MECHANICS.md section 8).

use super::dna::{melting_temp, revcom};
use super::features::{color_from_name, Feature};

/// `do_pcr` warns when a primer's first mismatch is nearer than this to its 3' end.
pub const PCR_WARN_BELOW: usize = 15;
/// On a circle, primers that overlap this little or less warn of linear amplification.
const OVERLAP_WARN_BASES: i64 = 10;
/// `pcr_dialog_precompute` repeats this many bases of a circular template (and of its reverse
/// complement) at the end, so a primer can anneal across the origin; a longer primer cannot.
const CIRCULAR_REPEAT: usize = 200;
/// The first four bases from a primer's 3' end must match exactly when no mismatches are allowed.
const EXACT_3_PRIME_WITHOUT_MISMATCHES: usize = 4;
/// The window name keeps this many characters of the template's name and of each primer's name.
const NAME_TEMPLATE_CHARS: usize = 16;
const NAME_PRIMER_CHARS: usize = 9;
/// `Tm` called without conditions: 250 nM DNA, 50 mM salt (its signature, @0xa72d3).
const TM_DNA_NM: f64 = 250.0;
const TM_SALT_MM: f64 = 50.0;
/// The colours of a `primer_bind` feature, forward and reverse (`feature_default_type_colors`, L67026).
const PRIMER_FORWARD_COLOR: &str = "#14c0bd";
const PRIMER_REVERSE_COLOR: &str = "#4ec02b";

pub const NOT_FACING: &str = "Primers don't point towards each other.";
const REVERSE_IN_FORWARD: &str =
    "Warning: Reverse primer amplifies part of forward primer. Primer dimers are very likely.";
const FORWARD_IN_REVERSE: &str =
    "Warning: Forwrd primer amplifies part of reverse primer. Primer dimers are very likely.";
const OVERLAP_SHORT: &str = "Warning: Reverse primer overlaps forward primer and extends less than 10bp from forward primer. Linear amplification is likely.";

/// How a primer is matched to the template: the dialog's three sliders (`pcr_dialog` @0x28656f,
/// L55413–55418).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct PcrScan {
    /// Shortest annealed stretch, counted from the primer's 3' end (`pcr_min_match`, 5 to 30).
    pub min_match: usize,
    /// Mismatches allowed inside it (`pcr_max_mm`, 0 to 10).
    pub max_mismatch: usize,
    /// Bases at the 3' end that must match exactly when mismatches are allowed (`pcr_min_3_match`, 0 to 4).
    pub min_3_match: usize,
}

impl Default for PcrScan {
    fn default() -> Self {
        PcrScan {
            min_match: 17,
            max_mismatch: 0,
            min_3_match: 4,
        }
    }
}

/// A primer as typed: its name and its bases, 5' to 3'.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct PcrPrimer {
    pub name: String,
    pub seq: Vec<u8>,
}

/// One place a primer anneals: a row of the dialog's table (`Direction`, `Location`, `Matchlength`, `MM`).
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
pub struct Site {
    /// `>>>`: the primer reads like the top strand (a forward primer). `<<<`: it pairs with the top
    /// strand (a reverse primer).
    pub forward: bool,
    /// 0-based. Forward: the top-strand position just after the annealed part (the 3' end). Reverse:
    /// where the annealed part starts.
    pub location: i64,
    /// Bases annealed, counted from the primer's 3' end (`Matchlength`).
    pub matched: i64,
    pub mismatches: usize,
}

// ---------------------------------------------------------------------------
// Finding the sites (`search_with_mm3`)
// ---------------------------------------------------------------------------

/// The template as `pcr_dialog_precompute` lists it: upper case, and for a circle followed by its
/// first bases again.
fn strand_list(text: &[u8], circular: bool) -> Vec<u8> {
    let mut list = text.to_ascii_uppercase();
    if circular {
        let again = list[..list.len().min(CIRCULAR_REPEAT)].to_vec();
        list.extend(again);
    }
    list
}

/// What the primer's base `c` may meet in the template; `search_with_mm3` expands the codes of a
/// degenerate primer this way.
fn letters(c: u8, degenerate: bool) -> Vec<u8> {
    if !degenerate {
        return vec![c];
    }
    match c {
        b'B' => b"BYSKCGT".to_vec(),
        b'D' => b"DRWKAGT".to_vec(),
        b'H' => b"HMWYACT".to_vec(),
        b'K' => b"KGT".to_vec(),
        b'M' => b"MAC".to_vec(),
        b'N' => b"NBDHVKMRYWSACGT".to_vec(),
        b'R' => b"RAG".to_vec(),
        b'S' => b"SGC".to_vec(),
        b'V' => b"VMSRACG".to_vec(),
        b'W' => b"WAT".to_vec(),
        b'Y' => b"YCT".to_vec(),
        other => vec![other],
    }
}

/// `search_with_mm3`: every start `i` in `list` where `s` (the primer's reverse complement, so that it
/// is read from the 3' end) anneals, with the annealed length `j` and the mismatches `mm`. A list
/// position outside the template reads as nothing, as `lindex` does; nothing equals only nothing.
fn search(
    s: &[u8],
    list: &[u8],
    template_len: usize,
    circular: bool,
    scan: &PcrScan,
) -> Vec<(i64, i64, usize)> {
    let degenerate = s.iter().any(|c| !b"ACGT".contains(c));
    let min_3 = if scan.max_mismatch > 0 {
        scan.min_3_match
    } else {
        EXACT_3_PRIME_WITHOUT_MISMATCHES
    };
    if s.len() < scan.min_match || s.len() < min_3 || template_len == 0 {
        return Vec::new();
    }
    let primer: Vec<Vec<u8>> = s.iter().map(|&c| letters(c, degenerate)).collect();
    let base = |i: i64| usize::try_from(i).ok().and_then(|i| list.get(i)).copied();
    let wanted = |j: i64| usize::try_from(j).ok().and_then(|j| primer.get(j));
    let mismatch = |k: usize, template: Option<u8>| match template {
        None => true,
        Some(t) => !primer[k].contains(&t),
    };
    // the base the walk back from the end of the stretch stops at: one the primer has as itself
    let same = |j: i64, i: i64| match (wanted(j), base(i + j)) {
        (None, None) => true,
        (Some(p), Some(t)) => p.len() == 1 && p[0] == t,
        _ => false,
    };
    // the first base in `from` that breaks the limit ends the walk; `j` counts the bases before it
    let walk = |i: i64, from: usize, mut j: i64| {
        let mut mm = 0;
        for k in from..primer.len() {
            if mismatch(k, base(i + k as i64)) {
                mm += 1;
            }
            if mm > scan.max_mismatch {
                mm -= 1;
                break;
            }
            j += 1;
        }
        (j, mm)
    };
    let needed = scan.min_match as i64 - 1;
    let mut found = Vec::new();
    if min_3 == 0 || degenerate {
        let upper = if circular {
            list.len() as i64 - CIRCULAR_REPEAT as i64
        } else {
            list.len() as i64 - primer.len() as i64
        };
        for i in 0..=upper {
            let (mut j, mm) = walk(i, 0, 0);
            if j >= needed {
                while !same(j, i) {
                    j -= 1;
                }
                found.push((i, j, mm));
            }
        }
    } else {
        for i in exact_starts(&s[..min_3], &list[..template_len], circular) {
            let (mut j, mm) = walk(i, min_3, min_3 as i64 - 1);
            if j >= needed {
                while !same(j, i) {
                    j -= 1;
                }
                found.push((i, j + 1, mm));
            }
        }
    }
    found
}

/// The starts, in order, where `head` occurs in `text` (across the origin of a circle): what the
/// original's precomputed dictionaries of the first one to four bases hold.
fn exact_starts(head: &[u8], text: &[u8], circular: bool) -> Vec<i64> {
    let n = text.len();
    (0..n)
        .filter(|&p| {
            head.iter().enumerate().all(|(d, &c)| {
                let q = p + d;
                let q = if circular {
                    Some(q % n)
                } else {
                    (q < n).then_some(q)
                };
                q.is_some_and(|q| text[q] == c)
            })
        })
        .map(|p| p as i64)
        .collect()
}

/// Every place `primer` anneals to the template, both directions, as the dialog's scan lists them.
pub fn find_sites(template: &[u8], circular: bool, primer: &[u8], scan: &PcrScan) -> Vec<Site> {
    let top = strand_list(template, circular);
    let bottom = strand_list(&revcom(&template.to_ascii_uppercase()), circular);
    let s = revcom(&primer.to_ascii_uppercase());
    let n = template.len() as i64;
    let mut sites = Vec::new();
    for (forward, list) in [(false, &top), (true, &bottom)] {
        for (i, j, mismatches) in search(&s, list, template.len(), circular, scan) {
            let location = if forward { n - i } else { i };
            sites.push(Site {
                forward,
                location,
                matched: j,
                mismatches,
            });
        }
    }
    sites
}

// ---------------------------------------------------------------------------
// The product (`pcr_dialog_do_pcr`)
// ---------------------------------------------------------------------------

/// What the original's PCR makes of two annealed primers: the new document's bases, name, comment and
/// the two primer features, and the warnings it gives first.
#[derive(Clone, Debug, PartialEq)]
pub struct Product {
    pub seq: Vec<u8>,
    /// Where the annealed parts lie on the template (0-based, end exclusive; the start is past the end
    /// when the part runs through the origin of a circle).
    pub fwd: (usize, usize),
    pub rev: (usize, usize),
    pub warnings: Vec<String>,
    /// The new document's name: `PCR`, the template, the two primers.
    pub name: String,
    /// The new document's file comment, one line after the other.
    pub comment: Vec<String>,
    /// The two `primer_bind` features, forward then reverse.
    pub features: Vec<Feature>,
    /// The stretch `[from, to)` of the template that the product copies (what `do_pcr` asks `get_tags2`
    /// to take the features from), when it does not run through the origin.
    pub template_span: Option<(usize, usize)>,
    /// Where that stretch starts in the product: the length of the forward primer's tail.
    pub template_offset: usize,
}

impl Product {
    pub fn length(&self) -> usize {
        self.seq.len()
    }
}

/// The bases from `a` up to `b`: a position before the start counts as the start and one past the end as
/// the end, an empty or backwards stretch is empty (`textarea_get` through `bp2ix`).
fn between(text: &[u8], a: i64, b: i64) -> Vec<u8> {
    let n = text.len() as i64;
    let (a, b) = (a.clamp(0, n), b.clamp(0, n));
    if b <= a {
        return Vec::new();
    }
    text[a as usize..b as usize].to_vec()
}

/// Bases that agree from the 3' end, comparing without regard to case: the length that `do_pcr` resets
/// the annealed stretch to ("the first mismatch between template and primer").
fn exact_from_3_prime(primer: &[u8], binding: &[u8]) -> i64 {
    let upper = |c: &u8| c.to_ascii_uppercase();
    primer
        .iter()
        .rev()
        .zip(binding.iter().rev())
        .take_while(|(p, b)| upper(p) == upper(b))
        .count() as i64
}

/// A melting temperature as the original prints it in a comment: one decimal, `---` over 100 bases, and
/// 0.0 for no bases at all.
fn tm_text(seq: &[u8]) -> String {
    if seq.is_empty() {
        return "0.0".to_string();
    }
    melting_temp(seq, TM_DNA_NM, TM_SALT_MM).map_or("---".to_string(), |t| format!("{t:.1}"))
}

fn tm_line(name: &str, seq: &[u8], exact: i64) -> String {
    let shown = String::from_utf8_lossy(seq);
    let tail = &seq[seq.len().saturating_sub(exact.max(0) as usize)..];
    if seq.len() as i64 == exact {
        format!("{name} {shown} Tm: {}°C.", tm_text(seq))
    } else {
        format!(
            "{name} {shown} Tm1: {}°C Tm2: {}°C.",
            tm_text(tail),
            tm_text(seq)
        )
    }
}

fn primer_feature(primer: &PcrPrimer, start: usize, end: usize, reverse: bool) -> Feature {
    let mut f = Feature::new(&primer.name, "primer_bind", start, end, reverse);
    if let (Some(fwd), Some(rev)) = (
        color_from_name(PRIMER_FORWARD_COLOR),
        color_from_name(PRIMER_REVERSE_COLOR),
    ) {
        (f.fwd_color, f.rev_color) = (fwd, rev);
    }
    let sequence = String::from_utf8_lossy(&primer.seq).into_owned();
    f.qualifiers.push((
        "PCR_conditions".into(),
        format!("primer sequence:{sequence}"),
    ));
    f
}

/// `pcr_dialog_do_pcr`: the product of the forward primer at `f` and the reverse primer at `r`.
/// `template_name` is the template's name without its extension. Errs when a linear template's
/// primers do not face each other.
pub fn make_product(
    template: &[u8],
    circular: bool,
    template_name: &str,
    fwd: &PcrPrimer,
    f: &Site,
    rev: &PcrPrimer,
    r: &Site,
) -> Result<Product, String> {
    let text_len = template.len() as i64;
    let len = text_len.max(1);
    let (fi, ri0) = (f.location, r.location);

    // reset each annealed stretch to the bases that agree from the 3' end
    let fwd_binding = if circular && fi - f.matched < 0 {
        let mut b = between(template, (fi - f.matched).rem_euclid(len), text_len);
        b.extend(between(template, fi - f.matched, fi));
        b
    } else {
        between(template, fi - f.matched, fi)
    };
    let fj = exact_from_3_prime(&fwd.seq, &fwd_binding.to_ascii_uppercase());
    let rev_binding = if circular && ri0 + r.matched > len {
        let mut b = between(template, ri0, text_len);
        b.extend(between(template, 0, (ri0 + r.matched).rem_euclid(len)));
        revcom(&b)
    } else {
        revcom(&between(template, ri0, ri0 + r.matched))
    };
    let rj = exact_from_3_prime(&rev.seq, &rev_binding.to_ascii_uppercase());

    let mut warnings = Vec::new();
    for (which, exact) in [("forward", fj), ("reverse", rj)] {
        if exact < PCR_WARN_BELOW as i64 {
            warnings.push(format!(
                "Warning: First {which} primer mismatch {exact} bases from 3' end."
            ));
        }
    }

    let (dna_first, dna_last);
    let mut ri = ri0;
    if circular {
        dna_first = (fi - fj).rem_euclid(len);
        if ri < fi {
            ri += len;
        }
        if (ri + rj) - (fi - fj) <= len {
            dna_last = (ri + rj).rem_euclid(len);
        } else {
            // round the world with primers that overlap: the overlap is copied (a negative end)
            dna_last = (ri + rj).rem_euclid(len) - len;
            if ri - (fi - fj) > len {
                warnings.push(REVERSE_IN_FORWARD.to_string());
            } else if ri + OVERLAP_WARN_BASES - (fi - fj) > len {
                warnings.push(OVERLAP_SHORT.to_string());
            }
            if (ri + rj) - fi > len {
                warnings.push(FORWARD_IN_REVERSE.to_string());
            } else if (ri + rj) - (fi - OVERLAP_WARN_BASES) > len {
                warnings.push(OVERLAP_SHORT.to_string());
            }
        }
    } else {
        if fi >= ri {
            return Err(NOT_FACING.to_string());
        }
        (dna_first, dna_last) = (fi - fj, ri + rj);
    }

    let mut between_primers = if dna_last > dna_first {
        between(template, dna_first, dna_last)
    } else {
        // through the origin; a negative end means more than once round
        let mut d = between(template, dna_first, text_len);
        let again = dna_last.rem_euclid(len);
        if dna_last < 0 && again < dna_first {
            d.extend(between(template, 0, text_len));
        }
        d.extend(between(template, 0, again));
        d
    };
    let fwd_tail = &fwd.seq[..fwd.seq.len().saturating_sub(fj.max(0) as usize)];
    let rev_tail = revcom(&rev.seq[..rev.seq.len().saturating_sub(rj.max(0) as usize)]);
    let mut seq = fwd_tail.to_vec();
    seq.append(&mut between_primers);
    seq.extend(rev_tail);

    let mut comment = vec!["PCR:".to_string(), format!("template: {template_name}")];
    comment.extend(warnings.iter().cloned());
    comment.push(tm_line(&fwd.name, &fwd.seq, fj));
    comment.push(tm_line(&rev.name, &rev.seq, rj));
    comment.push(format!("Product length: {}", seq.len()));

    let short = |name: &str, chars: usize| name.chars().take(chars).collect::<String>();
    let name = format!(
        "PCR {} {} {}",
        short(template_name, NAME_TEMPLATE_CHARS),
        short(&fwd.name, NAME_PRIMER_CHARS),
        short(&rev.name, NAME_PRIMER_CHARS)
    );
    let features = vec![
        primer_feature(fwd, 0, fwd.seq.len(), false),
        primer_feature(
            rev,
            seq.len().saturating_sub(rev.seq.len()),
            seq.len(),
            true,
        ),
    ];
    let start = |p: i64| p.rem_euclid(len) as usize;
    let end = |p: i64| {
        if p > len {
            (p - len) as usize
        } else {
            p.max(0) as usize
        }
    };
    let template_span = (0 <= dna_first && dna_first < dna_last && dna_last <= text_len)
        .then_some((dna_first as usize, dna_last as usize));
    Ok(Product {
        seq,
        fwd: (start(fi - f.matched), end(fi)),
        rev: (start(ri0), end(ri0 + r.matched)),
        warnings,
        name,
        comment,
        features,
        template_span,
        template_offset: fwd_tail.len(),
    })
}

/// Every product of a forward and a reverse primer: each forward site with each reverse site that
/// `make_product` accepts, shortest first.
pub fn pcr_products(
    template: &[u8],
    circular: bool,
    template_name: &str,
    fwd: &PcrPrimer,
    rev: &PcrPrimer,
    scan: &PcrScan,
) -> Vec<Product> {
    let forward: Vec<Site> = find_sites(template, circular, &fwd.seq, scan)
        .into_iter()
        .filter(|s| s.forward)
        .collect();
    let reverse: Vec<Site> = find_sites(template, circular, &rev.seq, scan)
        .into_iter()
        .filter(|s| !s.forward)
        .collect();
    let mut products: Vec<Product> = forward
        .iter()
        .flat_map(|f| reverse.iter().map(move |r| (f, r)))
        .filter_map(|(f, r)| make_product(template, circular, template_name, fwd, f, rev, r).ok())
        .collect();
    products.sort_by_key(Product::length);
    products
}

#[cfg(test)]
mod tests;

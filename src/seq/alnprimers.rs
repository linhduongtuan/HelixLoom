//! Primers and TaqMan probes picked from a multiple alignment of DNA: oligos that sit in the stretches the sequences share, with a count
//! of the sequences each one fits.
//!
//! The method is the usual one for assays that must work on a family of sequences (strains of a virus, copies of a gene), and is
//! this program's own arrangement of parts it already has:
//!
//! 1. Every column has a commonest base. A column is *usable* when at least `min_identity` of all the rows (a gap counts against)
//!    carry that base, and the gaps are not the majority. With `min_identity` 1 every sequence carries the base.
//! 2. The consensus of the usable columns, with `N` in every other column, is the template. Primer and probe candidates must be made
//!    of A, C, G and T (`primer3::design` and `probe::design_assays` refuse any other letter), so no oligo can lie on an unusable
//!    column; the amplicon between two primers may.
//! 3. Pairs (and probes) are picked from that template with primer3's rules (`seq::primer3`) and the probe rules of `seq::probe`;
//!    positions are alignment columns. Product sizes are asked for in columns: in a sequence with a gap inside the amplicon the product
//!    is shorter, and the true sizes are reported per design.
//! 4. Each oligo is then laid against every row: a row it fits exactly, and a row it fits *tolerably* (at most `max_mismatches`
//!    mismatches, a gap counting as one, and for a primer none in its last five bases, where a mismatch stops the polymerase
//!    most). A design covers the rows that every one of its oligos fits tolerably. Designs are ranked by the rows they cover, then
//!    by penalty.
//!
//! Not done: degenerate (IUPAC) oligos, which would cover more sequences with one order; mismatch-aware melting temperatures (the Tm is
//! that of the consensus oligo, as primer3 reports it); a check that the oligos are unique in a genome, which needs one.

use super::alnstats;
use super::primer3::{self, Oligo, Pair, Settings as Primers};
use super::probe::{self, Settings as ProbeSettings};

/// The most columns that are searched at once. primer3's rules look at every window of every size, and the search runs while the
/// window waits: a longer alignment is narrowed to a region first.
pub const MAX_COLUMNS: usize = 50_000;
/// How many bases at the 3' end of a primer must match.
const END_BASES: usize = 5;
/// Designs asked of the search for each one wanted, so that ranking by coverage has some to choose from, and the most asked.
const RANK_SPARE: usize = 4;
const MOST_ASKED: usize = 80;

#[derive(Clone, Debug, PartialEq)]
pub struct Settings {
    pub primers: Primers,
    /// Also pick a probe between the two primers of each design.
    pub with_probe: bool,
    pub probe: ProbeSettings,
    /// The share of the rows, 0 to 1, that must carry a column's commonest base for the column to be used.
    pub min_identity: f64,
    /// The mismatches an oligo may have in a row and still count as fitting it.
    pub max_mismatches: usize,
    /// Search only these columns (counted from 0, the end not included); none: all of them.
    pub window: Option<(usize, usize)>,
}

impl Default for Settings {
    fn default() -> Self {
        Settings {
            primers: Primers::default(),
            with_probe: false,
            probe: ProbeSettings::default(),
            min_identity: 1.0,
            max_mismatches: 1,
            window: None,
        }
    }
}

/// One oligo placed on the alignment, with the rows it fits.
#[derive(Clone, Debug, PartialEq)]
pub struct Site {
    /// As it is ordered, 5' to 3'; `start` and `end` are alignment columns (counted from 0, the end not included).
    pub oligo: Oligo,
    /// Rows it matches in every column.
    pub exact: usize,
    /// Rows it fits with the mismatches allowed.
    pub tolerated: usize,
    /// Where it lies in the reference row's own sequence, counted from 1, both ends included.
    pub in_reference: Option<(usize, usize)>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Design {
    pub left: Site,
    pub right: Site,
    pub probe: Option<Site>,
    /// The pair's penalty, with the probe's when there is one.
    pub penalty: f64,
    /// Rows that every oligo of the design matches perfectly, and rows that every oligo fits tolerably.
    pub exact: usize,
    pub tolerated: usize,
    /// The rows (counted from 0) that the design does not cover.
    pub uncovered: Vec<usize>,
    /// The fewest and the most bases from the start of the left primer to the end of the right one, over the rows covered.
    pub product: Option<(usize, usize)>,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Report {
    pub rows: usize,
    pub columns: usize,
    /// Columns that oligos may lie on, and the longest unbroken stretch of them.
    pub usable: usize,
    pub longest_run: usize,
    /// The commonest base of every column (`-` where gaps are most, `N` where no column has a base).
    pub consensus: Vec<u8>,
    pub designs: Vec<Design>,
}

/// A base as the comparison sees it: upper case, RNA's U as T.
fn norm(b: u8) -> u8 {
    match b.to_ascii_uppercase() {
        b'U' => b'T',
        c => c,
    }
}

fn is_gap(b: u8) -> bool {
    matches!(b, b'-' | b'.')
}

struct Facts {
    /// The commonest of A, C, G, T (the first of those that tie); 0 when no row has one.
    base: u8,
    /// The share of all the rows that carry it.
    share: f64,
    gaps_most: bool,
}

fn facts(rows: &[Vec<u8>], k: usize) -> Facts {
    let mut count = [0usize; 4];
    let mut gaps = 0;
    for row in rows {
        match norm(row[k]) {
            b'A' => count[0] += 1,
            b'C' => count[1] += 1,
            b'G' => count[2] += 1,
            b'T' => count[3] += 1,
            c if is_gap(c) => gaps += 1,
            _ => {}
        }
    }
    let (best, most) = count
        .iter()
        .enumerate()
        .fold((0, 0), |a, (i, &c)| if c > a.1 { (i, c) } else { a });
    Facts {
        base: if most == 0 { 0 } else { b"ACGT"[best] },
        share: most as f64 / rows.len() as f64,
        gaps_most: gaps * 2 > rows.len(),
    }
}

/// The consensus of the usable columns, `N` elsewhere.
fn template(facts: &[Facts], s: &Settings) -> Vec<u8> {
    let (from, to) = s.window.unwrap_or((0, facts.len()));
    facts
        .iter()
        .enumerate()
        .map(|(k, f)| {
            let usable = (from..to).contains(&k)
                && !f.gaps_most
                && f.base != 0
                && f.share + 1e-9 >= s.min_identity;
            if usable {
                f.base
            } else {
                b'N'
            }
        })
        .collect()
}

/// For each row, whether `oligo` matches it exactly and whether it fits it with the mismatches allowed. `primer`: the 3' end must
/// match as well.
fn fits(
    oligo: &Oligo,
    template: &[u8],
    rows: &[Vec<u8>],
    s: &Settings,
    primer: bool,
) -> Vec<(bool, bool)> {
    let (a, b) = (oligo.start, oligo.end);
    // a right primer is the reverse complement of its columns: its 3' end is the left end of them
    let end = if oligo.reverse {
        a..(a + END_BASES).min(b)
    } else {
        b.saturating_sub(END_BASES).max(a)..b
    };
    rows.iter()
        .map(|row| {
            let off = |k: usize| norm(row[k]) != template[k];
            let wrong = (a..b).filter(|&k| off(k)).count();
            let end_wrong = primer && end.clone().any(off);
            (wrong == 0, wrong <= s.max_mismatches && !end_wrong)
        })
        .collect()
}

/// Where columns `a..b` of `row` lie in the row's own sequence (counted from 1, both ends included); none when the row has no base there.
fn in_row(row: &[u8], a: usize, b: usize) -> Option<(usize, usize)> {
    let bases = |part: &[u8]| part.iter().filter(|&&c| !is_gap(c)).count();
    let before = bases(&row[..a]);
    let inside = bases(&row[a..b]);
    (inside > 0).then_some((before + 1, before + inside))
}

fn site(oligo: &Oligo, fit: &[(bool, bool)], reference: &[u8]) -> Site {
    Site {
        in_reference: in_row(reference, oligo.start, oligo.end),
        exact: fit.iter().filter(|f| f.0).count(),
        tolerated: fit.iter().filter(|f| f.1).count(),
        oligo: oligo.clone(),
    }
}

/// A pair, with its probe when it has one, laid against the rows.
fn lay(
    (pair, probe, penalty): (Pair, Option<Oligo>, f64),
    template: &[u8],
    rows: &[Vec<u8>],
    s: &Settings,
) -> Design {
    let fit = |o: &Oligo, primer: bool| fits(o, template, rows, s, primer);
    let (fl, fr) = (fit(&pair.left, true), fit(&pair.right, true));
    let fp = probe.as_ref().map(|p| fit(p, false));
    let every = |pick: fn(&(bool, bool)) -> bool| -> Vec<usize> {
        (0..rows.len())
            .filter(|&r| pick(&fl[r]) && pick(&fr[r]) && fp.as_ref().is_none_or(|f| pick(&f[r])))
            .collect()
    };
    let covered = every(|f| f.1);
    let sizes = covered.iter().map(|&r| {
        rows[r][pair.left.start..pair.right.end]
            .iter()
            .filter(|&&c| !is_gap(c))
            .count()
    });
    let product = sizes.clone().min().zip(sizes.max());
    let reference = &rows[0];
    Design {
        left: site(&pair.left, &fl, reference),
        right: site(&pair.right, &fr, reference),
        probe: probe
            .as_ref()
            .zip(fp.as_ref())
            .map(|(o, f)| site(o, f, reference)),
        penalty,
        exact: every(|f| f.0).len(),
        tolerated: covered.len(),
        uncovered: (0..rows.len()).filter(|r| !covered.contains(r)).collect(),
        product,
    }
}

/// The best designs first: those that cover more rows, then those with the lower penalty; no more than `want`.
fn rank(mut designs: Vec<Design>, want: usize) -> Vec<Design> {
    designs.sort_by(|a, b| {
        b.tolerated
            .cmp(&a.tolerated)
            .then(a.penalty.total_cmp(&b.penalty))
    });
    designs.truncate(want);
    designs
}

fn check(rows: &[(String, Vec<u8>)], s: &Settings) -> Result<(), String> {
    if rows.len() < 2 {
        return Err(
            "An alignment needs at least two sequences to have something in common.".into(),
        );
    }
    let columns = rows[0].1.len();
    if rows.iter().any(|r| r.1.len() != columns) {
        return Err("The rows are not all of one length, so this is not an alignment.".into());
    }
    let letters: Vec<Vec<u8>> = rows.iter().map(|r| r.1.clone()).collect();
    if alnstats::is_protein(&letters) {
        return Err("Primers and probes are for DNA: this alignment looks like protein.".into());
    }
    if !(s.min_identity > 0.0 && s.min_identity <= 1.0) {
        return Err("The share of sequences that must agree is more than 0 and at most 1.".into());
    }
    let (from, to) = s.window.unwrap_or((0, columns));
    if from >= to || to > columns {
        return Err(format!(
            "The region {}–{} is not inside the {columns} columns of the alignment.",
            from + 1,
            to
        ));
    }
    if to - from > MAX_COLUMNS {
        return Err(format!(
            "{} columns are too many to search at once: choose a region of at most {MAX_COLUMNS}.",
            to - from
        ));
    }
    Ok(())
}

/// The designs for an alignment of DNA (rows of one length, gaps as `-` or `.`), best first. The first row is the reference whose own
/// numbering the report gives positions in.
pub fn design(rows: &[(String, Vec<u8>)], s: &Settings) -> Result<Report, String> {
    check(rows, s)?;
    let letters: Vec<Vec<u8>> = rows.iter().map(|r| r.1.clone()).collect();
    let columns = letters[0].len();
    let all: Vec<Facts> = (0..columns).map(|k| facts(&letters, k)).collect();
    let template = template(&all, s);
    let want = s.primers.num_return;
    let asked = Primers {
        num_return: (want * RANK_SPARE).clamp(want, MOST_ASKED.max(want)),
        ..s.primers.clone()
    };
    let found: Vec<(Pair, Option<Oligo>, f64)> = if s.with_probe {
        probe::design_assays(&template, None, &asked, &s.probe)
            .into_iter()
            .map(|a| (a.pair, Some(a.probe), a.penalty))
            .collect()
    } else {
        primer3::design(&template, None, &asked)
            .into_iter()
            .map(|p| {
                let penalty = p.penalty;
                (p, None, penalty)
            })
            .collect()
    };
    let laid = found
        .into_iter()
        .map(|f| lay(f, &template, &letters, s))
        .collect();
    let usable = template.iter().filter(|&&c| c != b'N').count();
    let longest_run = template
        .split(|&c| c == b'N')
        .map(<[u8]>::len)
        .max()
        .unwrap_or(0);
    Ok(Report {
        rows: rows.len(),
        columns,
        usable,
        longest_run,
        consensus: all
            .iter()
            .map(|f| match (f.gaps_most, f.base) {
                (true, _) => b'-',
                (_, 0) => b'N',
                (_, base) => base,
            })
            .collect(),
        designs: rank(laid, want),
    })
}

/// The consensus without its gaps, and for each column the place in it (counted from 0) of that column's base, none for a column that
/// the consensus leaves out.
pub fn ungapped(consensus: &[u8]) -> (Vec<u8>, Vec<Option<usize>>) {
    let mut bases = Vec::with_capacity(consensus.len());
    let at = consensus
        .iter()
        .map(|&c| {
            (c != b'-').then(|| {
                bases.push(c);
                bases.len() - 1
            })
        })
        .collect();
    (bases, at)
}

/// The designs as a table of tab-separated text, one line to an oligo: for ordering, and for a spreadsheet. `names` are the rows' names, in
/// the order of the alignment.
pub fn table_text(report: &Report, names: &[String]) -> String {
    let mut out = String::from(
        "design\trole\tsequence_5to3\tcolumns\tstrand\ttm\tgc_percent\trows_exact\trows_fit\tposition_in_first_row\tproduct_min\tproduct_max\trows_covered\trows_not_covered\n",
    );
    for (i, d) in report.designs.iter().enumerate() {
        let product = d
            .product
            .map_or(("-".to_string(), "-".to_string()), |(a, b)| {
                (a.to_string(), b.to_string())
            });
        let missing: Vec<&str> = d
            .uncovered
            .iter()
            .filter_map(|&r| names.get(r).map(String::as_str))
            .collect();
        let sites = [
            ("left primer", Some(&d.left)),
            ("right primer", Some(&d.right)),
            ("probe", d.probe.as_ref()),
        ];
        for (role, site) in sites.into_iter().filter_map(|(r, s)| Some((r, s?))) {
            let o = &site.oligo;
            let place = site
                .in_reference
                .map_or("-".to_string(), |(a, b)| format!("{a}-{b}"));
            out.push_str(&format!(
                "{}\t{role}\t{}\t{}-{}\t{}\t{:.1}\t{:.1}\t{}\t{}\t{place}\t{}\t{}\t{}\t{}\n",
                i + 1,
                String::from_utf8_lossy(&o.seq),
                o.start + 1,
                o.end,
                if o.reverse { '-' } else { '+' },
                o.tm,
                o.gc,
                site.exact,
                site.tolerated,
                product.0,
                product.1,
                d.tolerated,
                missing.join(","),
            ));
        }
    }
    out
}

#[cfg(test)]
mod tests;

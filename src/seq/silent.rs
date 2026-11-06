//! Silent Sites (Enzymes menu): where a restriction site could be made without changing the protein.
//!
//! The rules are the original's `silent_sites` (AppMain.tcl L35791 @0x1a0b2c) with `coordinate_translation` (L35922
//! @0x1a24be), `degen_transl` (L64442 @0x2ea304) and `triplet_transl` (L64424 @0x2e9fae). Each site is translated in its three frames with every
//! ambiguity code allowed, which gives a set of residues for each codon; those sets are looked for in the translation of the
//! DNA; a place that already is a site is dropped; and for the others the fewest changed bases that make the site and
//! keep the residues are worked out. `tests` checks the reports against the original's own, run on random sequences.

use super::dna::{iupac_mask, translate, STANDARD_CODE};
use super::enzymes::Enzyme;
use super::sitesearch::{flat_pattern, header, pattern_covers, upper};

#[cfg(test)]
mod tests;

/// How many bases of the first codon come before the site, in the order the original looks.
const FRAMES: [usize; 3] = [2, 1, 0];
/// Most places looked at before the search stops. A long sequence with every enzyme would otherwise give millions of rows.
/// (HelixLoom's own limit; the original has none.)
pub const MAX_PLACES: usize = 100_000;
/// Most variants made of one place; each ambiguity code that does not fit grows them four-fold at most. (Own limit.)
const MAX_VARIANTS: usize = 1024;
const STOP_BIT: u32 = 26;
const ALL_RESIDUES: u32 = 21;

/// What `triplet_transl` makes of three site letters: the residues those codons can code for in the standard code.
/// All 21 (twenty residues and the stop) is the original's `.`, any character.
#[derive(Clone, Copy, Debug, PartialEq, Eq)]
enum Residues {
    Any,
    Of(u32),
}

fn residue_bit(aa: u8) -> u32 {
    match aa {
        b'A'..=b'Z' => 1 << (aa - b'A'),
        b'*' => 1 << STOP_BIT,
        _ => 0,
    }
}

impl Residues {
    fn of_triplet(triplet: &[u8]) -> Residues {
        const BASES: &[u8; 4] = b"TCAG";
        let mut set = 0;
        for (n, &aa) in STANDARD_CODE.as_bytes().iter().enumerate() {
            let codon = [BASES[n / 16], BASES[n / 4 % 4], BASES[n % 4]];
            let fits = triplet
                .iter()
                .zip(codon)
                .all(|(&letter, base)| iupac_mask(letter) & iupac_mask(base) != 0);
            if fits {
                set |= residue_bit(aa);
            }
        }
        if set.count_ones() < ALL_RESIDUES {
            Residues::Of(set)
        } else {
            Residues::Any
        }
    }

    fn allows(self, aa: u8) -> bool {
        match self {
            Residues::Any => true,
            Residues::Of(set) => set & residue_bit(aa) != 0,
        }
    }
}

/// `degen_transl`: the site with two `N` on each side, read in each of the three frames, a set of residues per codon.
fn degenerate_translation(site: &[u8]) -> [Vec<Residues>; 3] {
    let mut padded = b"NN".to_vec();
    padded.extend(
        site.iter()
            .map(u8::to_ascii_uppercase)
            .filter(|c| b"ABCDGHKMNRSTVWY".contains(c)),
    );
    padded.extend_from_slice(b"NN");
    std::array::from_fn(|start| {
        padded
            .get(start..)
            .unwrap_or_default()
            .as_chunks::<3>()
            .0
            .iter()
            .map(|triplet| Residues::of_triplet(triplet))
            .collect()
    })
}

/// Where the sets match the protein, like `regexp -all`: a match is not looked at again, the next starts after its end.
fn scan(protein: &[u8], sets: &[Residues]) -> Vec<usize> {
    let n = sets.len();
    let mut found = Vec::new();
    let mut at = 0;
    while n > 0 && at + n <= protein.len() {
        let fits = sets.iter().zip(&protein[at..]).all(|(s, &a)| s.allows(a));
        if fits {
            found.push(at);
            at += n;
        } else {
            at += 1;
        }
    }
    found
}

/// The bases an ambiguity code stands for, in the order the original's string map lists them (`S` is `GC`).
fn bases(c: u8) -> Option<&'static [u8]> {
    Some(match c {
        b'A' => b"A",
        b'B' => b"CGT",
        b'C' => b"C",
        b'D' => b"AGT",
        b'G' => b"G",
        b'H' => b"ACT",
        b'K' => b"GT",
        b'M' => b"AC",
        b'N' => b"ACGT",
        b'R' => b"AG",
        b'S' => b"GC",
        b'T' => b"T",
        b'V' => b"ACG",
        b'W' => b"AT",
        b'Y' => b"CT",
        _ => return None,
    })
}

fn is_plain(codon: &[u8]) -> bool {
    codon.iter().all(|c| b"ACGT".contains(c))
}

/// Keep the original's base wherever the site's letter allows it, take the site's letter where not. The original asks
/// `string match "*X*"` of the bases the letter stands for, so a `*` or `?` in the DNA matches anything and is kept.
fn substitute(original: &[u8], site: &[u8]) -> Vec<u8> {
    original
        .iter()
        .zip(site)
        .map(|(&o, &s)| {
            let allowed = match bases(s) {
                Some(set) => set.contains(&o),
                None => s == o,
            };
            if allowed || o == b'*' || o == b'?' {
                o
            } else {
                s
            }
        })
        .collect()
}

/// The codons of a residue, as the original lists them (L35946). Its stop list holds TGG, which codes for Trp; it is kept.
fn codons_for(aa: u8) -> &'static str {
    match aa {
        b'R' => "CGN CGV CGD CGB CGH CGW CGS CGR CGY CGK CGM CGA CGC CGG CGT AGR AGA AGG",
        b'S' => "TCN TCV TCD TCB TCH TCW TCS TCR TCY TCK TCM TCA TCC TCG TCT AGY AGC AGT",
        b'L' => "CTN CTV CTD CTB CTH CTW CTS CTR CTY CTK CTM CTA CTC CTG CTT TTR TTA TTG",
        b'G' => "GGN GGV GGD GGB GGH GGW GGS GGR GGY GGK GGM GGA GGC GGG GGT",
        b'V' => "GTN GTV GTD GTB GTH GTW GTS GTR GTY GTK GTM GTA GTC GTG GTT",
        b'P' => "CCN CCV CCD CCB CCH CCW CCS CCR CCY CCK CCM CCA CCC CCG CCT",
        b'T' => "ACN ACV ACD ACB ACH ACW ACS ACR ACY ACK ACM ACA ACC ACG ACT",
        b'A' => "GCN GCV GCD GCB GCH GCW GCS GCR GCY GCK GCM GCA GCC GCG GCT",
        b'I' => "ATH ATW ATY ATM ATA ATC ATT",
        b'F' => "TTY TTC TTC",
        b'Y' => "TAY TAC TAT",
        b'C' => "TGY TGC TGT",
        b'H' => "CAY CAC CAT",
        b'N' => "AAY AAC AAT",
        b'D' => "GAY GAC GAT",
        b'*' => "TAR TAA TAG TGG",
        b'Q' => "CAR CAA CAG",
        b'K' => "AAR AAA AAG",
        b'E' => "GAR GAA GAG",
        b'W' => "TGG",
        b'M' => "ATG",
        _ => "",
    }
}

fn residue_of(codon: &[u8], table: &str) -> Vec<u8> {
    translate(&codon.to_ascii_uppercase(), table)
}

/// The first codon of the residue that the site's letters allow.
fn first_codon(residue: u8, site: &[u8]) -> Option<Vec<u8>> {
    codons_for(residue)
        .split_whitespace()
        .find(|codon| {
            site.len() == codon.len()
                && site
                    .iter()
                    .zip(codon.bytes())
                    .all(|(&p, c)| pattern_covers(p, c))
        })
        .map(|codon| codon.as_bytes().to_vec())
}

/// One codon of the original's DNA and the same codon of the site, laid over each other: the bases to write.
fn fit_codon(original: &[u8], site: &[u8], table: &str) -> Vec<u8> {
    let simple = if is_plain(site) {
        site.to_vec()
    } else if site == b"NNN" {
        original.to_vec()
    } else {
        substitute(original, site)
    };
    let wanted = residue_of(original, table);
    if residue_of(&simple, table) == wanted {
        return simple;
    }
    // changing single bases does not keep the residue: take another codon for it that the site allows
    match wanted.first().and_then(|&aa| first_codon(aa, site)) {
        None => b"xxx".to_vec(),
        Some(codon) if is_plain(&codon) => codon,
        Some(codon) => substitute(original, &codon),
    }
}

/// Every plain sequence an ambiguous one can be; the first letter changes fastest, as in the original.
fn expand_all(fitted: &[u8]) -> Vec<Vec<u8>> {
    if is_plain(fitted) {
        return vec![fitted.to_vec()];
    }
    let mut all: Vec<Vec<u8>> = vec![Vec::new()];
    for &c in fitted {
        let choices = bases(c).map_or_else(|| vec![c], <[u8]>::to_vec);
        let mut next = Vec::new();
        for &base in &choices {
            for stem in &all {
                next.push([stem.as_slice(), &[base]].concat());
            }
        }
        next.truncate(MAX_VARIANTS);
        all = next;
    }
    all
}

/// `coordinate_translation` (L35922): the sequences that carry the site and keep the residues of `original`, which is
/// whole codons; `site` is the site padded with `N` to the same length.
pub fn coordinate(original: &[u8], site: &[u8], table: &str) -> Vec<Vec<u8>> {
    let original = original.to_ascii_uppercase();
    let (old, new) = (original.as_chunks::<3>().0, site.as_chunks::<3>().0);
    let fitted: Vec<u8> = old
        .iter()
        .zip(new)
        .flat_map(|(o, s)| fit_codon(o, s, table))
        .collect();
    expand_all(&fitted)
}

/// What to search: the sequence, the stretch of it (`from..to`, 0-based), the genetic code, and the enzymes with the 1-based
/// starts of the sites each one already has in the whole sequence, **in order** (so those are not offered again).
pub struct Query<'a> {
    pub seq: &'a [u8],
    pub region: (usize, usize),
    pub table: &'a str,
    pub enzymes: &'a [&'a Enzyme],
    pub existing: &'a [Vec<usize>],
}

/// One place where an enzyme's site could be made.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Hit {
    pub enzyme: String,
    pub frame: usize,
    /// The first base of the site, counting from 1.
    pub site_start: usize,
    /// The first base of the first codon the site touches, counting from 1: what the list shows.
    pub position: usize,
    /// The codons the site touches, `start..end` (0-based), for selecting them.
    pub codons_at: (usize, usize),
    pub codons: Vec<u8>,
    pub variants: Vec<Vec<u8>>,
    /// The site written under the codons, in threes.
    pub layout: String,
    /// The residues of the codons, three spaces apart.
    pub residues: String,
}

#[derive(Clone, Debug, Default, PartialEq, Eq)]
pub struct Search {
    pub hits: Vec<Hit>,
    /// The search stopped at `MAX_PLACES`.
    pub truncated: bool,
}

struct Place {
    enzyme: usize,
    frame: usize,
    start: usize,
    first: usize,
    len: usize,
}

fn places(q: &Query, protein: &[u8]) -> (Vec<Place>, bool) {
    let mut found = Vec::new();
    let mut truncated = false;
    'enzymes: for (k, e) in q.enzymes.iter().enumerate() {
        let sets = degenerate_translation(&flat_pattern(e));
        for (sets, frame) in sets.iter().zip(FRAMES) {
            for first in scan(protein, sets) {
                let start = q.region.0 + 1 + 3 * first + frame;
                if q.existing
                    .get(k)
                    .is_some_and(|sites| sites.binary_search(&start).is_ok())
                {
                    continue;
                }
                if found.len() >= MAX_PLACES {
                    truncated = true;
                    break 'enzymes;
                }
                found.push(Place {
                    enzyme: k,
                    frame,
                    start,
                    first,
                    len: sets.len(),
                });
            }
        }
    }
    // the original sorts by the start of the site, equal starts stay in the order they were found
    found.sort_by_key(|p| p.start);
    (found, truncated)
}

/// `seq_len_mod` (L58532): 1 to `len`, zero is `len`.
fn wrap(x: usize, len: usize) -> usize {
    match x % len {
        0 => len,
        r => r,
    }
}

fn layout(frame: usize, flat: &[u8]) -> String {
    let mut text = " ".repeat(frame).into_bytes();
    text.extend_from_slice(flat);
    text.extend_from_slice(b"   ");
    in_threes(&text)
}

fn hit_of(q: &Query, p: &Place) -> Hit {
    let e = q.enzymes[p.enzyme];
    let flat = flat_pattern(e);
    let begin = q.region.0 + 3 * p.first;
    let codons = q.seq[begin..begin + 3 * p.len].to_vec();
    let mut site = vec![b'N'; p.frame];
    site.extend_from_slice(&flat);
    site.resize(site.len().div_ceil(3) * 3, b'N');
    let residues: Vec<String> = translate(&codons, q.table)
        .iter()
        .map(|&a| char::from(a).to_string())
        .collect();
    Hit {
        enzyme: e.name.clone(),
        frame: p.frame,
        site_start: p.start,
        position: wrap(p.start - p.frame, q.seq.len()),
        codons_at: (begin, begin + codons.len()),
        variants: coordinate(&codons, &site, q.table),
        layout: layout(p.frame, &flat),
        residues: residues.join("   "),
        codons,
    }
}

/// The places in the stretch where one of the enzymes' sites could be made without changing the protein, in order of
/// position.
pub fn search(q: &Query) -> Search {
    let (from, to) = q.region;
    if from > to || to > q.seq.len() {
        return Search::default();
    }
    let protein = translate(&q.seq[from..to], q.table);
    let (found, truncated) = places(q, &protein);
    Search {
        hits: found.iter().map(|p| hit_of(q, p)).collect(),
        truncated,
    }
}

/// Each base of a variant with whether it differs from the DNA it replaces (upper case, as the original shows it).
pub fn marks(variant: &[u8], codons: &[u8]) -> Vec<(char, bool)> {
    variant
        .iter()
        .enumerate()
        .map(|(i, &c)| {
            let same = codons.get(i).is_some_and(|o| o.eq_ignore_ascii_case(&c));
            (upper(c), !same)
        })
        .collect()
}

fn spaced(marks: &[(char, bool)]) -> String {
    let mut text = String::new();
    for (i, &(c, changed)) in marks.iter().enumerate() {
        if i > 0 && i % 3 == 0 {
            text.push(' ');
        }
        if changed {
            text.push_str(&format!("<{c}>"));
        } else {
            text.push(c);
        }
    }
    text
}

/// Whole codons with a space between them (`regexp -inline -all {(?:...)}`, joined): a trailing part of a codon is dropped.
pub fn in_threes(bases: &[u8]) -> String {
    let threes: Vec<String> = bases
        .as_chunks::<3>()
        .0
        .iter()
        .map(|c| String::from_utf8_lossy(c).into_owned())
        .collect();
    threes.join(" ")
}

/// What the original's window shows for one hit, a changed base between `<` and `>`.
pub fn render_hit(h: &Hit) -> String {
    let mut text = format!("\t{}\t{}\n", h.position, in_threes(&h.codons));
    for v in &h.variants {
        text.push_str(&format!("\t\t{}\n", spaced(&marks(v, &h.codons))));
    }
    text.push_str(&format!("\t{}\t{}\n", h.enzyme, h.layout));
    text.push_str(&format!("\t\t {}\n\n", h.residues));
    text
}

/// The whole text of the original's Silent Sites window.
pub fn report(date: &str, name: &str, from: usize, to: usize, hits: &[Hit]) -> String {
    let mut text = header(date, name, from, to, "Silent Sites");
    for h in hits {
        text.push_str(&render_hit(h));
    }
    if hits.is_empty() {
        text.push_str("No sites Found");
    }
    text
}

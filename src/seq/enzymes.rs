//! Restriction enzymes: site finding, Dam/Dcm blocking and digestion.
//! Rules follow findenzymes / enz_cut_offset / digest_make_fraglist of the original.

use super::dna::{base_matches, revcom};

/// The original's default enzyme set (116 commercial enzymes): `name=pattern`.
/// `^` marks the top-strand cut; `(a/b)` gives top/bottom cut offsets past the site.
const DEFAULT_ENZYMES: &str = "\
AatII=GACGT^C AbsI=CC^TCGAGG Acc65I=G^GTACC AciI=C^CGC AclI=AA^CGTT AfeI=AGC^GCT AflII=C^TTAAG \
AgeI=A^CCGGT AluI=AG^CT ApaI=GGGCC^C ApaLI=G^TGCAC ApoI=R^AATTY AscI=GG^CGCGCC AseI=AT^TAAT \
AvaI=C^YCGRG AvrII=C^CTAGG BamHI=G^GATCC BbsI=GAAGAC(2/6) BclI=T^GATCA BfaI=C^TAG BglII=A^GATCT \
BmgBI=CAC^GTC BsaAI=YAC^GTR BsaBI=GATNN^NNATC BsaI=GGTCTC(1/5) BsaWI=W^CCGGW BsiWI=C^GTACG \
BsmBI=CGTCTC(1/5) BspEI=T^CCGGA BspHI=T^CATGA BsrBI=CCG^CTC BsrGI=T^GTACA BssHII=G^CGCGC \
BssSI=C^ACGAG BstBI=TT^CGAA BstUI=CG^CG BstYI=R^GATCY BstZ17I=GTA^TAC Cac8I=GCN^NGC ClaI=AT^CGAT \
DpnI=GA^TC DraI=TTT^AAA EagI=C^GGCCG EcoRI=G^AATTC EcoRV=GAT^ATC FseI=GGCCGG^CC FspI=TGC^GCA \
HaeII=RGCGC^Y HaeIII=GG^CC HhaI=GCG^C HincII=GTY^RAC HindIII=A^AGCTT HpaI=GTT^AAC HpaII=C^CGG \
Hpy8I=GTN^NAC HpyCH4IV=A^CGT HpyCH4V=TG^CA KasI=G^GCGCC KpnI=GGTAC^C MauBI=CG^CGCGCG MboI=^GATC \
MfeI=C^AATTG MluCI=^AATT MluI=A^CGCGT MreI=CG^CCGGCG MscI=TGG^CCA MseI=T^TAA NaeI=GCC^GGC \
NarI=GG^CGCC NcoI=C^CATGG NdeI=CA^TATG NgoMIV=G^CCGGC NheI=G^CTAGC NlaIII=CATG^ NotI=GC^GGCCGC \
NruI=TCG^CGA NsiI=ATGCA^T NspI=RCATG^Y PacI=TTAAT^TAA PaqCI=CACCTGC(4/8) PciI=A^CATGT PmeI=GTTT^AAAC \
PmlI=CAC^GTG PshAI=GACNN^NNGTC PsiI=TTA^TAA PstI=CTGCA^G PvuI=CGAT^CG PvuII=CAG^CTG RsaI=GT^AC \
SacI=GAGCT^C SacII=CCGC^GG SalI=G^TCGAC SapI=GCTCTTC(1/4) SbfI=CCTGCA^GG ScaI=AGT^ACT \
SfiI=GGCCNNNN^NGGCC SfoI=GGC^GCC SgfI=GCGAT^CGC SgrAI=CR^CCGGYG SgrDI=CG^TCGACG SmaI=CCC^GGG \
SnaBI=TAC^GTA SpeI=A^CTAGT SphI=GCATG^C SrfI=GCCC^GGGC SspI=AAT^ATT StuI=AGG^CCT SwaI=ATTT^AAAT \
TaqI=T^CGA TatI=W^GTACW XbaI=T^CTAGA XcmI=CCANNNNN^NNNNTGG XhoI=C^TCGAG XmaI=C^CCGGG \
XmnI=GAANN^NNTTC ZraI=GAC^GTC";

/// Enzymes whose sites are blocked when they overlap a Dam (GATC) site: (name, shift1, shift2).
const DAM_BLOCKED: &[(&str, i32, i32)] = &[
    ("BclI", 1, 1),
    ("BsaBI", 0, 6),
    ("BspEI", -2, 4),
    ("BspHI", -2, 4),
    ("ClaI", -1, 3),
    ("MboI", 0, 0),
    ("NruI", -2, 4),
    ("TaqI", -2, 2),
    ("XbaI", -2, 4),
];

/// Same for Dcm (CCWGG) sites.
const DCM_BLOCKED: &[(&str, i32, i32)] = &[
    ("Acc65I", -3, 4),
    ("ApaI", -3, 4),
    ("BsaI", -3, 4),
    ("MscI", -2, 3),
    ("StuI", -2, 3),
];

/// Enzymes that only cut methylated DNA.
const NEEDS_METHYLATION: &[&str] = &["DpnI", "MalI"];

/// Bases copied from the start onto the end when a circular sequence is scanned (findenzymes).
const CIRCULAR_SCAN_EXTRA: usize = 63;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Enzyme {
    pub name: String,
    /// Recognition sequence, upper case IUPAC, without cut marks.
    pub site: Vec<u8>,
    /// Top / bottom strand cut offsets from the start of the site.
    pub top: i32,
    pub bottom: i32,
}

impl Enzyme {
    pub fn parse(name: &str, pattern: &str) -> Option<Enzyme> {
        let pattern = pattern.to_ascii_uppercase();
        let (site, top, bottom) = if let Some(open) = pattern.find('(') {
            let (a, b) = pattern[open + 1..].trim_end_matches(')').split_once('/')?;
            let len = open as i32;
            (
                pattern[..open].to_string(),
                len + a.trim().parse::<i32>().ok()?,
                len + b.trim().parse::<i32>().ok()?,
            )
        } else if let Some(caret) = pattern.find('^') {
            let flat: String = pattern.chars().filter(|&c| c != '^').collect();
            let len = flat.len() as i32;
            (flat, caret as i32, len - caret as i32)
        } else {
            (pattern.clone(), 0, 0)
        };
        let valid = !site.is_empty() && site.bytes().all(|c| super::dna::iupac_mask(c) != 0);
        valid.then(|| Enzyme {
            name: name.to_string(),
            site: site.into_bytes(),
            top,
            bottom,
        })
    }

    pub fn is_palindromic(&self) -> bool {
        revcom(&self.site) == self.site
    }

    /// Positive: 5' overhang, negative: 3' overhang, zero: blunt.
    pub fn overhang(&self) -> i32 {
        self.bottom - self.top
    }

    /// Number of fully specified bases (A/C/G/T) in the recognition sequence.
    pub fn cutter_size(&self) -> usize {
        self.site.iter().filter(|c| b"ACGT".contains(c)).count()
    }

    pub fn is_degenerate(&self) -> bool {
        self.site.iter().any(|c| !b"ACGT".contains(c))
    }

    pub fn site_text(&self) -> String {
        String::from_utf8_lossy(&self.site).into_owned()
    }

    /// The bases `insert_enzyme_site` writes: each ambiguity code becomes one fixed base
    /// (N W R D H M V -> A, S Y B -> C, K -> G); the reverse strand when asked.
    pub fn site_to_insert(&self, reverse: bool) -> Vec<u8> {
        let fixed: Vec<u8> = self
            .site
            .iter()
            .map(|&c| match c {
                b'N' | b'W' | b'R' | b'D' | b'H' | b'M' | b'V' => b'A',
                b'S' | b'Y' | b'B' => b'C',
                b'K' => b'G',
                other => other,
            })
            .collect();
        if reverse {
            revcom(&fixed)
        } else {
            fixed
        }
    }

    /// The pattern with its cut mark: `G^AATTC`, or `GGTCTC(1/5)` for cuts away from the site.
    pub fn pattern_text(&self) -> String {
        let len = self.site.len() as i32;
        let site = self.site_text();
        if self.top >= 0 && self.top <= len && self.bottom == len - self.top {
            let at = self.top as usize;
            format!("{}^{}", &site[..at], &site[at..])
        } else {
            format!("{site}({}/{})", self.top - len, self.bottom - len)
        }
    }
}

pub fn default_enzymes() -> Vec<Enzyme> {
    DEFAULT_ENZYMES
        .split_whitespace()
        .filter_map(|tok| tok.split_once('='))
        .filter_map(|(n, p)| Enzyme::parse(n, p))
        .collect()
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Site {
    /// Start of the recognition sequence (0-based).
    pub start: usize,
    pub reverse: bool,
    /// Cut positions as boundaries (a cut at `k` falls before base `k`), top and bottom strand.
    pub top_cut: usize,
    pub bottom_cut: usize,
    /// False for the one site per enzyme that the original lists although its cut lies a base before the start of a linear
    /// sequence (see `Cuts::Edge`): it is counted, but it has no cut, so no digest or map may use `top_cut` and `bottom_cut`
    /// (both 0 then).
    pub cuts_inside: bool,
}

/// Where the cuts of a site fall in a sequence of `n` bases.
enum Cuts {
    /// The two cuts, as boundaries.
    Inside(usize, usize),
    /// A reverse site of an enzyme that cuts away from its site, in a linear sequence, whose top-strand cut would fall one
    /// base before the start. The original's `findenzymes` takes out the reverse sites with `site < reach - length` (L43080-43094,
    /// the site numbered from 1), one more than it should, so this one stays in its list; its digest rejects it (L17128-17134).
    Edge,
}

/// The cuts of a site, or none when a linear sequence drops it for cutting beyond its ends.
fn cut_positions(
    e: &Enzyme,
    start: usize,
    reverse: bool,
    n: usize,
    circular: bool,
) -> Option<Cuts> {
    let len = e.site.len() as i32;
    let (t, b) = if reverse {
        (len - e.bottom, len - e.top)
    } else {
        (e.top, e.bottom)
    };
    let (t, b) = (start as i32 + t, start as i32 + b);
    if circular {
        let m = n as i32;
        Some(Cuts::Inside(
            t.rem_euclid(m) as usize,
            b.rem_euclid(m) as usize,
        ))
    } else if t >= 0 && b >= 0 && t <= n as i32 && b <= n as i32 {
        Some(Cuts::Inside(t as usize, b as usize))
    } else {
        // numbered from 1, a reverse site stays unless its number is below `reach - length`
        let reach = e.top.max(e.bottom);
        (reverse && start as i32 + 1 >= reach - len).then_some(Cuts::Edge)
    }
}

/// Recognition sequences longer than this are not scanned base by base: the original finds them with a regular-expression
/// search (`longnamelist` in the enzyme file, `findenzymes` L42960 @0x1f09a2), which goes on after the end of a site it has
/// found and so never reports one that overlaps it. The 21 enzymes of the default file that this applies to are exactly the
/// ones with more than six bases (the others go through the hash scan, L43001, which sees every start).
const HASH_SCAN_LIMIT: usize = 6;

/// The starts, from 0, where `pat` matches `text`, one strand's worth. Beyond `HASH_SCAN_LIMIT` the next start looked at is
/// the one after the end of the last match.
fn match_starts(pat: &[u8], text: &[u8], limit: usize) -> Vec<usize> {
    let skip = if pat.len() > HASH_SCAN_LIMIT {
        pat.len()
    } else {
        1
    };
    let mut starts = Vec::new();
    let mut start = 0;
    while start < limit {
        let hit = pat
            .iter()
            .zip(&text[start..])
            .all(|(&p, &s)| base_matches(p, s, false));
        if hit {
            starts.push(start);
        }
        start += if hit { skip } else { 1 };
    }
    starts
}

/// All sites of `e` in `seq`; linear sequences drop sites whose cuts fall outside the ends.
pub fn find_sites(e: &Enzyme, seq: &[u8], circular: bool) -> Vec<Site> {
    let n = seq.len();
    let k = e.site.len();
    if n < 4 || n < k && !circular {
        return Vec::new();
    }
    let text: Vec<u8> = if circular {
        seq.iter()
            .chain(seq.iter().take(CIRCULAR_SCAN_EXTRA))
            .copied()
            .collect()
    } else {
        seq.to_vec()
    };
    let rc = revcom(&e.site);
    let orientations: Vec<(&[u8], bool)> = if e.is_palindromic() {
        vec![(&e.site, false)]
    } else {
        vec![(&e.site, false), (&rc, true)]
    };
    let limit = n.min(text.len().saturating_sub(k - 1));
    let mut found: Vec<(usize, bool)> = Vec::new();
    for &(pat, reverse) in &orientations {
        let starts = match_starts(pat, &text, limit);
        found.extend(starts.into_iter().map(|start| (start, reverse)));
    }
    // by start; at one start the forward strand comes first, as the original sorts its two lists into one
    found.sort_by_key(|&(start, _)| start);
    found
        .into_iter()
        .filter_map(|(start, reverse)| {
            let (top_cut, bottom_cut, cuts_inside) =
                match cut_positions(e, start, reverse, n, circular)? {
                    Cuts::Inside(top, bottom) => (top, bottom, true),
                    Cuts::Edge => (0, 0, false),
                };
            Some(Site {
                start,
                reverse,
                top_cut,
                bottom_cut,
                cuts_inside,
            })
        })
        .collect()
}

fn methylase_starts(seq: &[u8], motif: &[u8], circular: bool) -> Vec<usize> {
    let n = seq.len();
    let text: Vec<u8> = if circular {
        seq.iter()
            .chain(seq.iter().take(motif.len() - 1))
            .copied()
            .collect()
    } else {
        seq.to_vec()
    };
    (0..n)
        .filter(|&i| i + motif.len() <= text.len())
        .filter(|&i| {
            motif
                .iter()
                .zip(&text[i..])
                .all(|(&p, &s)| base_matches(p, s, false))
        })
        .collect()
}

/// Sites for every enzyme (index aligned with `enzymes`). With `methylated`, sites overlapping
/// Dam/Dcm methylation are removed and DpnI-type enzymes are allowed to cut.
pub fn all_sites(
    enzymes: &[Enzyme],
    seq: &[u8],
    circular: bool,
    methylated: bool,
) -> Vec<Vec<Site>> {
    let dam = methylase_starts(seq, b"GATC", circular);
    let dcm = methylase_starts(seq, b"CCWGG", circular);
    let n = seq.len() as i32;
    let blocked = |name: &str, start: usize| -> bool {
        let hit = |table: &[(&str, i32, i32)], starts: &[usize]| {
            table
                .iter()
                .filter(|(nm, _, _)| *nm == name)
                .any(|&(_, a, b)| {
                    [a, b].iter().any(|&d| {
                        let p = start as i32 + d;
                        let p = if circular { p.rem_euclid(n) } else { p };
                        p >= 0 && starts.contains(&(p as usize))
                    })
                })
        };
        methylated && (hit(DAM_BLOCKED, &dam) || hit(DCM_BLOCKED, &dcm))
    };
    enzymes
        .iter()
        .map(|e| {
            if NEEDS_METHYLATION.contains(&e.name.as_str()) && !methylated {
                return Vec::new();
            }
            let mut sites = find_sites(e, seq, circular);
            sites.retain(|s| !blocked(&e.name, s.start));
            sites
        })
        .collect()
}

// ---------------------------------------------------------------------------
// Digestion
// ---------------------------------------------------------------------------

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Cut {
    pub pos: usize,
    pub label: String,
    /// Recognition-site start (0-based) the cut came from.
    pub site: usize,
}

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Fragment {
    pub length: usize,
    pub start: usize,
    pub end: usize,
    pub left: String,
    pub right: String,
    pub mass_pct: u32,
}

/// Cuts (top strand) for the chosen enzymes, merged where two enzymes cut at the same place.
pub fn collect_cuts(chosen: &[(&Enzyme, &[Site])]) -> Vec<Cut> {
    let mut cuts: Vec<Cut> = chosen
        .iter()
        .flat_map(|(e, sites)| {
            sites.iter().filter(|s| s.cuts_inside).map(|s| Cut {
                pos: s.top_cut,
                label: e.name.clone(),
                site: s.start,
            })
        })
        .collect();
    cuts.sort_by_key(|c| (c.pos, c.label.clone()));
    let mut merged: Vec<Cut> = Vec::new();
    for c in cuts {
        match merged.last_mut() {
            Some(last) if last.pos == c.pos => {
                if !last.label.split(" + ").any(|l| l == c.label) {
                    last.label = format!("{} + {}", last.label, c.label);
                }
            }
            _ => merged.push(c),
        }
    }
    merged
}

/// Fragments of a complete digest, largest first. Sizes are measured between top-strand cuts.
pub fn digest(n: usize, circular: bool, cuts: &[Cut]) -> Vec<Fragment> {
    if n == 0 {
        return Vec::new();
    }
    let make =
        |start: usize, end: usize, length: usize, left: &str, right: &str| -> Option<Fragment> {
            (length > 0).then(|| Fragment {
                length,
                start,
                end,
                left: left.into(),
                right: right.into(),
                mass_pct: ((100.0 * length as f64 / n as f64).round() as u32).max(1),
            })
        };
    let uncut = || make(0, n, n, "uncut", "").into_iter().collect::<Vec<_>>();
    let mut frags: Vec<Fragment> = if circular && cuts.is_empty() {
        uncut()
    } else if circular {
        cuts.iter()
            .enumerate()
            .filter_map(|(i, c)| {
                let next = &cuts[(i + 1) % cuts.len()];
                let length = if cuts.len() == 1 {
                    n
                } else {
                    (next.pos + n - c.pos) % n
                };
                make(
                    c.pos,
                    next.pos,
                    if length == 0 { n } else { length },
                    &c.label,
                    &next.label,
                )
            })
            .collect()
    } else {
        let mut bounds: Vec<(usize, &str)> = vec![(0, "Start")];
        bounds.extend(cuts.iter().map(|c| (c.pos, c.label.as_str())));
        bounds.push((n, "End"));
        let pieces: Vec<Fragment> = bounds
            .windows(2)
            .filter_map(|w| make(w[0].0, w[1].0, w[1].0 - w[0].0, w[0].1, w[1].1))
            .collect();
        if pieces.is_empty() {
            uncut()
        } else {
            pieces
        }
    };
    frags.sort_by_key(|f| std::cmp::Reverse(f.length));
    frags
}

#[cfg(test)]
mod tests {
    use super::*;

    fn enz(name: &str) -> Enzyme {
        default_enzymes()
            .into_iter()
            .find(|e| e.name == name)
            .unwrap()
    }

    #[test]
    fn inserted_sites_replace_ambiguity_codes_with_one_fixed_base() {
        assert_eq!(enz("EcoRI").site_to_insert(false), b"GAATTC");
        assert_eq!(enz("BsaI").site_to_insert(false), b"GGTCTC");
        assert_eq!(enz("BsaI").site_to_insert(true), b"GAGACC");
        // GTYRAC: Y -> C, R -> A
        assert_eq!(enz("HincII").site_to_insert(false), b"GTCAAC");
        // SfiI GGCCNNNNNGGCC: N -> A
        assert_eq!(enz("SfiI").site_to_insert(false), b"GGCCAAAAAGGCC");
    }

    #[test]
    fn the_default_set_has_all_116_enzymes() {
        let all = default_enzymes();
        assert_eq!(all.len(), 116);
        let eco = enz("EcoRI");
        assert_eq!(
            (eco.site_text().as_str(), eco.top, eco.bottom),
            ("GAATTC", 1, 5)
        );
        assert_eq!(eco.overhang(), 4);
        let bsai = enz("BsaI");
        assert_eq!(
            (bsai.site_text().as_str(), bsai.top, bsai.bottom),
            ("GGTCTC", 7, 11)
        );
        assert!(!bsai.is_palindromic());
        assert_eq!(enz("SfiI").site_text(), "GGCCNNNNNGGCC");
        assert_eq!(enz("PstI").overhang(), -4);
        assert_eq!(enz("EcoRV").overhang(), 0);
        assert_eq!(enz("MboI").top, 0);
        assert_eq!(enz("NlaIII").top, 4);
    }

    #[test]
    fn patterns_print_with_their_cut_marks() {
        let marks = |name: &str| enz(name).pattern_text();
        assert_eq!(marks("EcoRI"), "G^AATTC");
        assert_eq!(marks("MboI"), "^GATC");
        assert_eq!(marks("NlaIII"), "CATG^");
        assert_eq!(marks("BsaI"), "GGTCTC(1/5)");
        assert_eq!(marks("PaqCI"), "CACCTGC(4/8)");
        assert_eq!(marks("SfiI"), "GGCCNNNN^NGGCC");
        for e in default_enzymes() {
            let again = Enzyme::parse(&e.name, &e.pattern_text()).unwrap();
            assert_eq!(again, e, "{} did not round-trip", e.name);
        }
    }

    #[test]
    fn finds_sites_on_both_strands_and_cut_positions() {
        let mut seq = vec![b'A'; 60];
        seq[3..9].copy_from_slice(b"GAATTC");
        seq[12..18].copy_from_slice(b"GGTCTC");
        seq[40..46].copy_from_slice(b"GAGACC");
        let eco = find_sites(&enz("EcoRI"), &seq, false);
        assert_eq!(
            eco,
            vec![Site {
                start: 3,
                reverse: false,
                top_cut: 4,
                bottom_cut: 8,
                cuts_inside: true
            }]
        );
        let bsai = find_sites(&enz("BsaI"), &seq, false);
        assert_eq!(bsai.len(), 2);
        assert!(!bsai[0].reverse && bsai[1].reverse);
        // forward GGTCTC at 12 cuts 1 base past the site on top (19) and 5 on the bottom (23)
        assert_eq!((bsai[0].top_cut, bsai[0].bottom_cut), (19, 23));
        // reverse site GAGACC at 40 cuts upstream of it: top 40+6-11, bottom 40+6-7
        assert_eq!((bsai[1].top_cut, bsai[1].bottom_cut), (35, 39));
    }

    /// A linear sequence of As with GAGACC (the reverse of BsaI's GGTCTC) put in at `at`.
    fn reverse_bsai_at(at: usize) -> Vec<u8> {
        let mut seq = vec![b'A'; 40];
        seq[at..at + 6].copy_from_slice(b"GAGACC");
        seq
    }

    #[test]
    fn a_reverse_site_that_cuts_before_the_start_is_dropped_but_the_one_the_original_keeps_has_no_cut(
    ) {
        // BsaI cuts 5 bases upstream of a reverse site on the top strand: the cut fits from the sixth base on. The original
        // takes out the sites numbered from 1 below 5, so the site at 0-based 4, whose cut would be a base outside, stays
        let sites = |at: usize| find_sites(&enz("BsaI"), &reverse_bsai_at(at), false);
        assert!(sites(3).is_empty(), "three bases before it: out");
        let edge = sites(4);
        assert_eq!(edge.len(), 1, "four bases before it: listed");
        assert!(!edge[0].cuts_inside && (edge[0].top_cut, edge[0].bottom_cut) == (0, 0));
        let inside = sites(5);
        assert_eq!(inside.len(), 1);
        assert!(inside[0].cuts_inside);
        assert_eq!((inside[0].top_cut, inside[0].bottom_cut), (0, 4));
    }

    #[test]
    fn the_edge_site_is_counted_but_makes_no_cut_in_a_digest() {
        let bsai = enz("BsaI");
        let seq = reverse_bsai_at(4);
        let sites = find_sites(&bsai, &seq, false);

        let cuts = collect_cuts(&[(&bsai, &sites)]);

        assert_eq!(sites.len(), 1, "it is a site of the sequence");
        assert!(cuts.is_empty(), "but it cuts nowhere");
        assert_eq!(digest(seq.len(), false, &cuts).len(), 1);
    }

    #[test]
    fn the_edge_site_is_found_on_the_other_type_iis_enzymes_at_their_own_distance() {
        // the start (from 0) where each enzyme's reverse site is the edge one: reach - length - 1
        let cases: [(&str, &[u8], usize); 4] = [
            ("BbsI", b"GTCTTC", 5),
            ("BsmBI", b"GAGACG", 4),
            ("PaqCI", b"GCAGGTG", 7),
            ("SapI", b"GAAGAGC", 3),
        ];
        for (name, site, at) in cases {
            let mut seq = vec![b'A'; 50];
            seq[at..at + site.len()].copy_from_slice(site);
            let found = find_sites(&enz(name), &seq, false);
            assert!(found.len() == 1 && !found[0].cuts_inside, "{name}");
            let mut one_less = vec![b'A'; 50];
            one_less[at - 1..at - 1 + site.len()].copy_from_slice(site);
            assert!(
                find_sites(&enz(name), &one_less, false).is_empty(),
                "{name}"
            );
        }
    }

    #[test]
    fn a_circular_sequence_has_no_edge_site_and_a_forward_site_near_the_end_is_dropped_as_before() {
        let circular = find_sites(&enz("BsaI"), &reverse_bsai_at(4), true);
        assert!(circular.iter().all(|s| s.cuts_inside));
        // forward GGTCTC cuts 11 bases from its start on the bottom strand: it needs them
        let mut seq = vec![b'A'; 40];
        seq[40 - 11..40 - 5].copy_from_slice(b"GGTCTC");
        assert_eq!(find_sites(&enz("BsaI"), &seq, false).len(), 1);
        seq = vec![b'A'; 40];
        seq[40 - 10..40 - 4].copy_from_slice(b"GGTCTC");
        assert!(find_sites(&enz("BsaI"), &seq, false).is_empty());
    }

    #[test]
    fn finds_degenerate_sites() {
        let seq = b"TTGAATTTTCAAAGGGCCCAAA";
        assert_eq!(find_sites(&enz("ApoI"), seq, false).len(), 1);
        assert_eq!(find_sites(&enz("SfiI"), b"GGCCAAAATGGCC", false).len(), 1);
    }

    fn starts(name: &str, seq: &[u8]) -> Vec<usize> {
        let sites = find_sites(&enz(name), seq, false);
        sites.iter().map(|s| s.start).collect()
    }

    #[test]
    fn an_enzyme_with_more_than_six_bases_misses_a_site_that_overlaps_the_one_before() {
        // MauBI CGCGCGCG occurs at 4 and, overlapping it by six bases, at 6: the original's regular-expression search goes
        // on after the first and finds only that one
        let seq = b"AAAACGCGCGCGCGAAAA";
        assert_eq!(starts("MauBI", seq), [4]);
        // a site that merely follows the first is found
        assert_eq!(starts("MauBI", b"AAAACGCGCGCGCGCGCGCGAAAA"), [4, 12]);
        // NotI GCGGCCGC overlaps itself by two bases
        assert_eq!(starts("NotI", b"AAAAGCGGCCGCGGCCGCAAAA"), [4]);
    }

    #[test]
    fn an_enzyme_with_six_bases_or_fewer_finds_every_start() {
        // BstUI CGCG at 2, 4 and 6 in CGCGCGCG... (the hash scan sees each start)
        assert_eq!(starts("BstUI", b"AACGCGCGCGAA"), [2, 4, 6]);
        // AatII GACGTC cannot overlap itself; one site
        assert_eq!(starts("AatII", b"AAGACGTCAA"), [2]);
    }

    #[test]
    fn the_two_strands_of_a_long_enzyme_are_scanned_apart() {
        // SapI GCTCTTC and its reverse GAAGAGC overlap each other by two bases here: both are found (the cuts are well
        // inside the ends), in the order of their starts
        let seq = b"AAAAAAAAGAAGAGCTCTTCAAAAAAAA";
        let sites = find_sites(&enz("SapI"), seq, false);
        let found: Vec<(usize, bool)> = sites.iter().map(|s| (s.start, s.reverse)).collect();
        assert_eq!(found, [(8, true), (13, false)]);
    }

    #[test]
    fn circular_sequences_find_sites_across_the_origin() {
        let seq = b"ATTCAAAAAAAAAAAAGA";
        assert!(find_sites(&enz("EcoRI"), seq, false).is_empty());
        let sites = find_sites(&enz("EcoRI"), seq, true);
        assert_eq!(sites.len(), 1);
        assert_eq!((sites[0].start, sites[0].top_cut), (16, 17));
    }

    #[test]
    fn linear_sequences_drop_sites_that_cut_past_the_end() {
        // GGTCTC at the very end would cut beyond the sequence on the top strand
        assert!(find_sites(&enz("BsaI"), b"AAAAAAAAAGGTCTC", false).is_empty());
        assert_eq!(
            find_sites(&enz("BsaI"), b"AAAAAAAAAGGTCTCAAAAAAAAA", false).len(),
            1
        );
    }

    #[test]
    fn dam_methylation_blocks_overlapping_sites() {
        let seq = b"AAAAAATCTAGATCAAAAAAAAAA"; // XbaI TCTAGA followed by TC -> GATC overlaps
        let xba = [enz("XbaI")];
        assert_eq!(all_sites(&xba, seq, false, false)[0].len(), 1);
        assert_eq!(all_sites(&xba, seq, false, true)[0].len(), 0);
        let clear = b"AAAAAATCTAGAAAAAAAAAAAAA";
        assert_eq!(all_sites(&xba, clear, false, true)[0].len(), 1);
        let dpn = [enz("DpnI")];
        assert_eq!(all_sites(&dpn, b"AAAAGATCAAAA", false, false)[0].len(), 0);
        assert_eq!(all_sites(&dpn, b"AAAAGATCAAAA", false, true)[0].len(), 1);
    }

    #[test]
    fn dcm_methylation_blocks_stui() {
        let blocked = b"AAAACCAGGCCTAAAAAAAA"; // CCAGG + CCT: StuI AGGCCT overlaps CCWGG
        let stu = [enz("StuI")];
        assert_eq!(all_sites(&stu, blocked, false, false)[0].len(), 1);
        assert_eq!(all_sites(&stu, blocked, false, true)[0].len(), 0);
    }

    fn digest_of(seq: &[u8], circular: bool, names: &[&str]) -> Vec<Fragment> {
        let enzymes: Vec<Enzyme> = names.iter().map(|n| enz(n)).collect();
        let sites = all_sites(&enzymes, seq, circular, false);
        let chosen: Vec<(&Enzyme, &[Site])> = enzymes
            .iter()
            .zip(sites.iter().map(|s| s.as_slice()))
            .collect();
        digest(seq.len(), circular, &collect_cuts(&chosen))
    }

    #[test]
    fn circular_digest_gives_fragments_that_sum_to_the_whole() {
        let mut seq = vec![b'A'; 200];
        seq[20..26].copy_from_slice(b"GAATTC");
        seq[120..126].copy_from_slice(b"GAATTC");
        let f = digest_of(&seq, true, &["EcoRI"]);
        assert_eq!(
            f.iter().map(|x| x.length).collect::<Vec<_>>(),
            vec![100, 100]
        );
        assert_eq!(f.iter().map(|x| x.length).sum::<usize>(), 200);
        assert!(f.iter().all(|x| x.mass_pct == 50));
    }

    #[test]
    fn linear_digest_adds_start_and_end() {
        let mut seq = vec![b'A'; 300];
        seq[100..106].copy_from_slice(b"GGATCC");
        let f = digest_of(&seq, false, &["BamHI"]);
        // BamHI cuts after the first G: 101 bases to the left, 199 to the right
        assert_eq!(
            f.iter().map(|x| x.length).collect::<Vec<_>>(),
            vec![199, 101]
        );
        assert_eq!(f[0].left, "BamHI");
        assert_eq!(f[0].right, "End");
        assert_eq!(f[1].left, "Start");
    }

    #[test]
    fn a_single_cut_in_a_circle_gives_one_full_length_fragment() {
        let mut seq = vec![b'A'; 120];
        seq[50..56].copy_from_slice(b"GAATTC");
        let f = digest_of(&seq, true, &["EcoRI"]);
        assert_eq!(f.len(), 1);
        assert_eq!(f[0].length, 120);
    }

    #[test]
    fn uncut_sequences_report_the_whole_length() {
        let seq = vec![b'A'; 90];
        assert_eq!(digest_of(&seq, true, &["EcoRI"])[0].left, "uncut");
        assert_eq!(digest_of(&seq, false, &["EcoRI"])[0].length, 90);
    }

    #[test]
    fn enzymes_cutting_the_same_place_are_merged() {
        // SmaI CCC^GGG and XmaI C^CCGGG share the site CCCGGG but cut at different offsets
        let mut seq = vec![b'T'; 100];
        seq[40..46].copy_from_slice(b"CCCGGG");
        let f = digest_of(&seq, false, &["SmaI", "XmaI"]);
        assert_eq!(f.len(), 3);
        // isoschizomers with identical cuts: HpaII C^CGG and MspI-like are merged by label
        let seq2 = {
            let mut s = vec![b'T'; 50];
            s[10..14].copy_from_slice(b"CCGG");
            s
        };
        let cuts = collect_cuts(&[
            (&enz("HpaII"), &find_sites(&enz("HpaII"), &seq2, false)),
            (&enz("HpaII"), &find_sites(&enz("HpaII"), &seq2, false)),
        ]);
        assert_eq!(cuts.len(), 1);
        assert_eq!(cuts[0].label, "HpaII");
    }
}

//! The tests of typing, and the made-up schemes and genomes that the tests of the parts use: no real scheme is in the repository (PubMLST's
//! data may not be bundled), so a scheme is made here of random alleles that differ from each other by a few bases.

use super::call::{type_contigs, Call, Params, Status};
use super::scheme::Scheme;
use super::search::Contigs;
use crate::seq::dna::revcom;
use crate::seq::minimap::sim::{genome, Rng};

pub const LENGTH: usize = 450;

/// The sequence of allele `id` (1 to ...) of locus `l`: a base sequence of the locus with `id - 1` bases changed, a different set for each id.
pub fn allele(l: usize, id: usize) -> Vec<u8> {
    allele_of(1000, l, id)
}

/// `allele`, for another scheme: the locus's base sequence comes from `base`.
pub fn allele_of(base: u32, l: usize, id: usize) -> Vec<u8> {
    let mut s = genome(LENGTH, base + l as u32);
    for j in 0..id - 1 {
        let at = (j * 61 + id * 17 + 29) % (LENGTH - 60) + 30;
        s[at] = if s[at] == b'A' { b'C' } else { b'A' };
    }
    s
}

/// A scheme of `loci` loci (`L0`, `L1` ...) with `alleles` alleles each, and five sequence types: type `s` has the allele `(s + l) % alleles + 1`
/// at locus `l`.
pub fn scheme(name: &str, loci: usize, alleles: usize) -> Scheme {
    scheme_from(name, 1000, loci, alleles)
}

/// `scheme`, with the loci's base sequences from `base` on.
pub fn scheme_from(name: &str, base: u32, loci: usize, alleles: usize) -> Scheme {
    let mut table = String::from("ST");
    (0..loci).for_each(|l| table += &format!("\tL{l}"));
    table += "\tclonal_complex\n";
    for s in 1..=5 {
        table += &s.to_string();
        (0..loci).for_each(|l| table += &format!("\t{}", (s + l) % alleles + 1));
        table += "\tCC\n";
    }
    let fastas = |locus: &str| {
        let l: usize = locus.strip_prefix('L')?.parse().ok()?;
        Some(
            (1..=alleles)
                .map(|id| {
                    format!(
                        ">{locus}_{id}\n{}\n",
                        String::from_utf8(allele_of(base, l, id)).unwrap()
                    )
                })
                .collect::<String>(),
        )
    };
    Scheme::from_texts(name, &table, &fastas).unwrap()
}

/// What to put in a genome: the sequence, whether on the reverse strand, and the contig.
pub struct Put {
    pub seq: Vec<u8>,
    pub rev: bool,
    pub contig: usize,
}

pub fn put(seq: Vec<u8>, rev: bool, contig: usize) -> Put {
    Put { seq, rev, contig }
}

/// A genome of three contigs of random bases (about 30 kb each, from `seed`) with the sequences put in them at places that are spread out.
pub fn contigs(seed: u32, puts: &[Put]) -> Vec<(String, Vec<u8>)> {
    let mut c: Vec<Vec<u8>> = (0..3).map(|i| genome(30_000, seed * 10 + i)).collect();
    let mut r = Rng::new(seed);
    let mut used: Vec<(usize, usize, usize)> = Vec::new();
    for p in puts {
        let seq = if p.rev { revcom(&p.seq) } else { p.seq.clone() };
        // a place that does not touch the others
        let at = loop {
            let at = 1000 + r.below(c[p.contig].len() - 2000 - seq.len());
            if used
                .iter()
                .all(|&(k, a, b)| k != p.contig || at + seq.len() + 500 < a || at > b + 500)
            {
                break at;
            }
        };
        used.push((p.contig, at, at + seq.len()));
        c[p.contig].splice(at..at + seq.len(), seq);
    }
    c.into_iter()
        .enumerate()
        .map(|(i, s)| (format!("contig{}", i + 1), s))
        .collect()
}

/// Put allele `ids[l]` of every locus `l` in the genome, on alternate strands.
pub fn with_alleles(ids: &[usize], seed: u32) -> Vec<(String, Vec<u8>)> {
    let puts: Vec<Put> = ids
        .iter()
        .enumerate()
        .map(|(l, &id)| put(allele(l, id), l % 2 == 1, l % 3))
        .collect();
    contigs(seed, &puts)
}

/// The scheme written as a folder of the layout that `Scheme::load` reads.
pub fn write_scheme(dir: &std::path::Path, s: &Scheme) {
    let folder = dir.join(&s.name);
    std::fs::create_dir_all(&folder).unwrap();
    let mut table = String::from("ST");
    s.loci
        .iter()
        .for_each(|l| table += &format!("\t{}", l.name));
    table.push('\n');
    for (st, alleles) in &s.profiles {
        table += &format!("{st}\t{}\n", alleles.join("\t"));
    }
    std::fs::write(folder.join(format!("{}.txt", s.name)), table).unwrap();
    for l in &s.loci {
        let fasta: String = l
            .alleles
            .iter()
            .map(|a| {
                format!(
                    ">{}_{}\n{}\n",
                    l.name,
                    a.id,
                    String::from_utf8_lossy(&a.seq)
                )
            })
            .collect();
        std::fs::write(folder.join(format!("{}.tfa", l.name)), fasta).unwrap();
    }
}

/// A genome as a FASTA file's text.
pub fn fasta(contigs: &[(String, Vec<u8>)]) -> String {
    contigs
        .iter()
        .map(|(n, s)| format!(">{n} a contig\n{}\n", String::from_utf8_lossy(s)))
        .collect()
}

fn typed(s: &Scheme, c: Vec<(String, Vec<u8>)>) -> super::call::Typing {
    type_contigs(s, &Contigs::new(c).unwrap(), &Params::default())
}

fn texts(t: &super::call::Typing) -> Vec<String> {
    t.calls.iter().map(|(_, c)| c.text()).collect()
}

#[test]
fn a_genome_with_the_alleles_of_a_sequence_type_is_typed_as_it_with_all_the_points() {
    let s = scheme("synth", 7, 6);
    // type 2: allele (2 + l) % 6 + 1 at locus l
    let ids: Vec<usize> = (0..7).map(|l| (2 + l) % 6 + 1).collect();
    let t = typed(&s, with_alleles(&ids, 1));
    assert_eq!(
        texts(&t),
        ids.iter().map(usize::to_string).collect::<Vec<_>>()
    );
    assert_eq!(
        (t.sequence_type.as_deref(), t.status, t.scheme.as_str()),
        (Some("2"), Status::Perfect, "synth")
    );
    // seven exact loci and a known type: 90 and 10
    assert!((t.score - 100.0).abs() < 1e-9, "{}", t.score);
}

#[test]
fn exact_alleles_that_are_in_no_profile_are_a_novel_combination_with_no_type_and_ninety_points() {
    let s = scheme("synth", 7, 6);
    // every locus has its allele 1: no type has that
    let t = typed(&s, with_alleles(&[1; 7], 2));
    assert_eq!(texts(&t), ["1"; 7]);
    assert_eq!((t.sequence_type, t.status), (None, Status::Novel));
    assert!((t.score - 90.0).abs() < 1e-9);
}

#[test]
fn an_allele_that_is_not_known_but_close_is_novel_and_named_for_the_closest() {
    let s = scheme("synth", 7, 6);
    let ids: Vec<usize> = (0..7).map(|l| (2 + l) % 6 + 1).collect();
    let mut puts: Vec<Put> = ids
        .iter()
        .enumerate()
        .map(|(l, &id)| put(allele(l, id), l % 2 == 1, l % 3))
        .collect();
    // locus 2: its allele with three more bases changed (99.3 % alike)
    let mut changed = allele(2, ids[2]);
    for at in [100, 200, 300] {
        changed[at] = if changed[at] == b'G' { b'T' } else { b'G' };
    }
    puts[2] = put(changed.clone(), false, 2);
    let t = typed(&s, contigs(3, &puts));
    assert_eq!(t.calls[2].1.text(), format!("~{}", ids[2]));
    assert_eq!(t.sequence_type, None);
    // six exact and a novel: (6 × 90 + 63) / 7
    assert!(
        (t.score - (6.0 * 90.0 + 63.0) / 7.0).abs() < 1e-9,
        "{}",
        t.score
    );
    assert_eq!(t.status, Status::Ok);
    // the sequence that was found is the changed allele, as it was written
    let Call::Novel(found) = &t.calls[2].1 else {
        panic!("a novel call")
    };
    assert_eq!(found[0].1, changed);
}

#[test]
fn an_allele_cut_short_by_the_end_of_a_contig_is_partial_and_one_that_is_not_there_is_missing() {
    let s = scheme("synth", 7, 6);
    let ids: Vec<usize> = (0..7).map(|l| (2 + l) % 6 + 1).collect();
    let mut puts: Vec<Put> = ids
        .iter()
        .enumerate()
        .skip(2)
        .map(|(l, &id)| put(allele(l, id), l % 2 == 1, l % 3))
        .collect();
    // locus 0: only the first 300 of its 450 bases, at the very end of a contig
    let mut c = contigs(4, &puts);
    let shorter = c[0].1.len() - 1;
    c[0].1.truncate(shorter);
    let part = allele(0, ids[0])[..300].to_vec();
    let end = c[0].1.len();
    c[0].1.splice(end - 300..end, part);
    puts.clear();
    let t = typed(&s, c);
    assert_eq!(t.calls[0].1.text(), format!("{}?", ids[0]));
    // locus 1 is nowhere
    assert_eq!(t.calls[1].1, Call::Missing);
    assert_eq!(t.sequence_type, None);
    assert!(
        (t.score - (5.0 * 90.0 + 18.0) / 7.0).abs() < 1e-9,
        "{}",
        t.score
    );
    assert_eq!(t.status, Status::Missing);
}

#[test]
fn two_different_alleles_at_one_locus_are_both_called_in_order_and_there_is_no_type() {
    let s = scheme("synth", 7, 6);
    let ids: Vec<usize> = (0..7).map(|l| (2 + l) % 6 + 1).collect();
    let mut puts: Vec<Put> = ids
        .iter()
        .enumerate()
        .map(|(l, &id)| put(allele(l, id), l % 2 == 1, l % 3))
        .collect();
    // locus 3 has allele 6 and allele 2 on two contigs
    puts.push(put(allele(3, 6), false, 1));
    puts.push(put(allele(3, 2), true, 2));
    let t = typed(&s, contigs(5, &puts));
    let mut both = vec![ids[3], 6, 2];
    both.sort_unstable();
    both.dedup();
    assert_eq!(
        t.calls[3].1.text(),
        both.iter()
            .map(usize::to_string)
            .collect::<Vec<_>>()
            .join(",")
    );
    assert_eq!((t.sequence_type, t.status), (None, Status::Mixed));
    // the same allele twice, on both strands, is one allele
    let again = typed(
        &s,
        contigs(
            6,
            &[put(allele(0, 4), false, 0), put(allele(0, 4), true, 1)],
        ),
    );
    assert_eq!(again.calls[0].1.text(), "4");
}

#[test]
fn an_allele_that_is_too_unlike_is_missing_until_the_limit_is_lowered() {
    let s = scheme("synth", 1, 6);
    let mut changed = allele(0, 3);
    // forty-four bases of 450 (10 %) changed, spread out
    // (the first and the last bases are left as they are: the end of an alignment that does not match is left out of it)
    for at in (15..LENGTH - 15).step_by(10) {
        changed[at] = if changed[at] == b'A' { b'C' } else { b'A' };
    }
    let c = contigs(7, &[put(changed, false, 0)]);
    let strict = type_contigs(&s, &Contigs::new(c.clone()).unwrap(), &Params::default());
    assert_eq!(strict.calls[0].1, Call::Missing);
    let lenient = type_contigs(
        &s,
        &Contigs::new(c).unwrap(),
        &Params {
            min_identity: 85.0,
            ..Params::default()
        },
    );
    // novel, and named for whichever of the alleles it is closest to
    assert!(
        lenient.calls[0].1.text().starts_with('~'),
        "{}",
        lenient.calls[0].1.text()
    );
}

#[test]
fn a_single_unknown_base_makes_a_novel_allele_and_an_allele_split_over_two_contigs_is_partial() {
    let s = scheme("synth", 1, 6);
    let mut with_n = allele(0, 4);
    with_n[200] = b'N';
    let t = typed(&s, contigs(8, &[put(with_n, false, 0)]));
    assert_eq!(t.calls[0].1.text(), "~4");
    // the first 200 bases at the end of one contig, the rest at the start of another
    let a = allele(0, 5);
    let mut c = contigs(9, &[]);
    let n = c[0].1.len();
    c[0].1.splice(n - 200..n, a[..200].to_vec());
    c[1].1.splice(0..250, a[200..].to_vec());
    let t = typed(&s, c);
    assert_eq!(t.calls[0].1.text(), "5?");
}

#[test]
fn an_allele_with_a_deletion_is_found_exactly_and_the_search_is_quick_enough_for_a_real_scheme() {
    // a scheme with alleles of two lengths: allele 2 lacks three bases
    let s = scheme("synth", 7, 60);
    let mut short = allele(0, 2);
    short.drain(200..203);
    let hit = super::search::find(
        &Contigs::new(contigs(10, &[put(short, false, 0)])).unwrap(),
        &s.loci[0],
    );
    // the deleted allele is not among the alleles: it is found as the closest, with one gap
    assert!(hit.iter().all(|f| f.exact.is_none()));
    // seven loci of sixty alleles against 90 kb
    let ids: Vec<usize> = (0..7).map(|l| (3 + l) % 60 + 1).collect();
    let started = std::time::Instant::now();
    let t = typed(&s, with_alleles(&ids, 11));
    assert!(started.elapsed().as_secs() < 20, "{:?}", started.elapsed());
    assert_eq!(t.sequence_type.as_deref(), Some("3"));
}

/// A contig of 20 kb of random bases with sequences put in at the places given (the reverse strand when the flag is set).
fn contig_with(name: &str, seed: u32, puts: &[(usize, &[u8], bool)]) -> (String, Vec<u8>) {
    let mut c = genome(20_000, seed);
    for &(at, seq, rev) in puts {
        let s = if rev { revcom(seq) } else { seq.to_vec() };
        c.splice(at..at + s.len(), s);
    }
    (name.to_string(), c)
}

fn changed(seq: &[u8], at: &[usize]) -> Vec<u8> {
    let mut s = seq.to_vec();
    for &i in at {
        s[i] = if s[i] == b'A' { b'C' } else { b'A' };
    }
    s
}

#[test]
fn places_that_overlap_in_coordinates_but_not_in_strand_or_contig_stay_apart() {
    let s = scheme("synth", 1, 6);
    // allele 1 on the reverse strand low on a contig, allele 2 forward high on it: two places
    let same = vec![contig_with(
        "c1",
        21,
        &[(2000, &allele(0, 1), true), (12_000, &allele(0, 2), false)],
    )];
    assert_eq!(texts(&typed(&s, same)), ["1,2"]);
    // allele 1 high on one contig, allele 2 low on another: two places
    let apart = vec![
        contig_with("c1", 22, &[(15_000, &allele(0, 1), false)]),
        contig_with("c2", 23, &[(2000, &allele(0, 2), false)]),
    ];
    assert_eq!(texts(&typed(&s, apart)), ["1,2"]);
}

#[test]
fn an_allele_that_holds_another_is_found_as_the_longer_one() {
    // allele 2 is the middle of allele 1: a genome with allele 1 has allele 2 in it too, and the call is the whole of it
    let long = allele(0, 1);
    let short = long[100..350].to_vec();
    let fastas = |l: &str| {
        (l == "L0").then(|| {
            format!(
                ">L0_1\n{}\n>L0_2\n{}\n",
                String::from_utf8_lossy(&long),
                String::from_utf8_lossy(&short)
            )
        })
    };
    let s = Scheme::from_texts("nested", "ST\tL0\n1\t1\n2\t2\n", &fastas).unwrap();
    let t = typed(&s, vec![contig_with("c1", 24, &[(5000, &long, false)])]);
    assert_eq!(texts(&t), ["1"]);
    assert_eq!(t.sequence_type.as_deref(), Some("1"));
}

#[test]
fn the_closest_allele_is_named_among_more_than_the_few_that_are_aligned() {
    // twelve alleles; the genome has allele 12 with one base changed, and allele 12 is the closest
    let s = scheme("synth", 1, 12);
    let novel = changed(&allele(0, 12), &[5]);
    let t = typed(&s, vec![contig_with("c1", 25, &[(8000, &novel, false)])]);
    assert_eq!(texts(&t), ["~12"]);
}

#[test]
fn of_the_alleles_aligned_the_one_with_most_bases_the_same_is_the_closest() {
    // allele 1 differs from the genome's allele at 4 bases far apart, allele 2 at 10 bases together: allele 2 shares more of its
    // minimizers with it (so it is the first to be aligned), but allele 1 is the closer
    let g = genome(LENGTH, 5000);
    let p = changed(&g, &[80, 160, 240, 320]);
    let q = changed(&g, &[150, 151, 152, 153, 154, 155, 156, 157, 158, 159]);
    let fastas = |l: &str| {
        (l == "L0").then(|| {
            format!(
                ">L0_1\n{}\n>L0_2\n{}\n",
                String::from_utf8_lossy(&p),
                String::from_utf8_lossy(&q)
            )
        })
    };
    let s = Scheme::from_texts("close", "ST\tL0\n1\t1\n2\t2\n", &fastas).unwrap();
    let t = typed(&s, vec![contig_with("c1", 26, &[(6000, &g, false)])]);
    assert_eq!(texts(&t), ["~1"]);
}

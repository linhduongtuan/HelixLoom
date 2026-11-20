use super::*;
use crate::seq::dna::revcom;
use std::collections::BTreeMap;

struct Rng(u64);

impl Rng {
    fn next(&mut self) -> u64 {
        self.0 = self
            .0
            .wrapping_mul(6_364_136_223_846_793_005)
            .wrapping_add(1_442_695_040_888_963_407);
        self.0 >> 33
    }

    fn below(&mut self, n: usize) -> usize {
        self.next() as usize % n
    }

    fn bases(&mut self, n: usize) -> Vec<u8> {
        (0..n).map(|_| b"ACGT"[self.below(4)]).collect()
    }
}

/// The number that a k-mer of A, C, G, T has as the counter writes it (two bits a base, the first base highest).
fn number(kmer: &[u8]) -> u64 {
    kmer.iter().fold(0, |n, &b| {
        (n << 2) | b"ACGT".iter().position(|&c| c == b).expect("a base") as u64
    })
}

/// The counts by the plainest way there is, on strings: each window of letters that are all bases, upper case, U as T, the smaller of
/// it and its reverse complement as text (which orders as the numbers do).
fn reference(seqs: &[&[u8]], k: usize) -> BTreeMap<u64, u32> {
    let mut out = BTreeMap::new();
    for seq in seqs {
        let seq: Vec<u8> = seq
            .iter()
            .map(|b| b.to_ascii_uppercase())
            .map(|b| if b == b'U' { b'T' } else { b })
            .collect();
        for w in seq.windows(k) {
            if w.iter().all(|b| b"ACGT".contains(b)) {
                let rc = revcom(w);
                *out.entry(number(w.min(&rc[..]))).or_insert(0) += 1;
            }
        }
    }
    out
}

#[test]
fn the_counts_are_those_of_the_plain_way_for_every_length_with_other_letters_and_lower_case() {
    let mut rng = Rng(5);
    for k in [1usize, 2, 3, 5, 11, 21, 27, 31, 32] {
        let mut seqs: Vec<Vec<u8>> = Vec::new();
        for _ in 0..6 {
            let length = 40 + rng.below(60);
            seqs.push(rng.bases(length));
        }
        // letters that break windows, lower case, U, and a sequence shorter than k
        seqs[0][10] = b'N';
        seqs[1][20] = b'R';
        seqs[2] = seqs[2].to_ascii_lowercase();
        seqs[3][5] = b'U';
        seqs.push(b"AC".to_vec());
        let refs: Vec<&[u8]> = seqs.iter().map(Vec::as_slice).collect();
        let c = Counts::count(&refs, k).expect("counted");
        let want = reference(&refs, k);
        assert_eq!(c.distinct(), want.len(), "k {k}");
        assert_eq!(
            c.total(),
            want.values().map(|&n| u64::from(n)).sum::<u64>(),
            "k {k}"
        );
        for (kmer, n) in &want {
            assert_eq!(c.get(*kmer), *n, "k {k}, k-mer {kmer}");
        }
    }
}

#[test]
fn a_hand_worked_case_a_sequence_and_its_reverse_complement_are_one_and_a_letter_that_is_no_base_breaks_windows(
) {
    // ACG and CGT are reverse complements: four windows of one k-mer; GTA and TAC: two of another
    let c = Counts::count(&[b"ACGTACGT"], 3).expect("counted");
    assert_eq!(
        (c.get(number(b"ACG")), c.get(number(b"GTA")), c.distinct()),
        (4, 2, 2)
    );
    // the other strand is the same k-mers
    let both = Counts::count(&[b"ACGTTGCAAGCT", &revcom(b"ACGTTGCAAGCT")], 4).expect("counted");
    let one = Counts::count(&[b"ACGTTGCAAGCT"], 4).expect("counted");
    assert_eq!(both.distinct(), one.distinct());
    assert_eq!(both.total(), 2 * one.total());
    // an N takes out the windows that hold it, and every window has its place
    let w = windows(b"ACGNACG", 3);
    assert_eq!(w.len(), 5);
    assert_eq!(
        w.iter().map(Option::is_some).collect::<Vec<_>>(),
        [true, false, false, false, true]
    );
    assert_eq!(Counts::count(&[b"ACGNACG"], 3).expect("counted").total(), 2);
    assert!(windows(b"AC", 3).is_empty());
    assert!(windows(b"ACGT", 0).is_empty() && windows(b"ACGT", 33).is_empty());
}

#[test]
fn the_length_of_the_k_mer_must_fit_in_sixty_four_bits() {
    for k in [0usize, 33, 100] {
        assert!(
            check_k(k).is_err() && Counts::count(&[b"ACGT"], k).is_err(),
            "{k}"
        );
    }
    assert!(check_k(1).is_ok() && check_k(32).is_ok());
}

#[test]
fn the_spectrum_the_gc_matrix_and_the_comparison_are_what_they_are_worked_out_to_be() {
    // k = 2 on AAAAAAAC: AA six times, AC once
    let a = Counts::count(&[b"AAAAAAAC"], 2).expect("counted");
    let s = spectrum(&a, 1, 10);
    assert_eq!((s[1], s[6], s.iter().sum::<u64>()), (1, 1, 2));
    // the last bin holds those that occur as often or more; fewer than `low` are left out
    let cut = spectrum(&a, 1, 4);
    assert_eq!((cut[1], cut[4], cut.len()), (1, 1, 5));
    assert_eq!(spectrum(&a, 2, 10).iter().sum::<u64>(), 1);
    // AA has no G or C, AC has one
    let g = gc_matrix(&a, 10);
    assert_eq!((g[0][6], g[1][1], g.len()), (1, 1, 3));
    assert_eq!(g.iter().flatten().sum::<u64>(), 2);
    // against ACACAC: AC three times, CA twice
    let b = Counts::count(&[b"ACACAC"], 2).expect("counted");
    let m = comparison(&a, &b, 10, 10).expect("compared");
    assert_eq!((m[6][0], m[1][3], m[0][2]), (1, 1, 1));
    assert_eq!(m.iter().flatten().sum::<u64>(), 3);
    // the last row and column hold all the larger; another k is refused
    let capped = comparison(&a, &b, 3, 2).expect("compared");
    assert_eq!((capped[3][0], capped[1][2], capped[0][2]), (1, 1, 1));
    let other = Counts::count(&[b"ACGT"], 3).expect("counted");
    assert!(comparison(&a, &other, 5, 5)
        .expect_err("another k")
        .contains("different k"));
}

#[test]
fn a_gc_count_is_the_same_for_a_k_mer_and_its_reverse_complement() {
    let mut rng = Rng(9);
    for _ in 0..50 {
        let kmer = rng.bases(21);
        let c = Counts::count(&[&kmer], 21).expect("counted");
        let r = Counts::count(&[&revcom(&kmer)], 21).expect("counted");
        assert_eq!(gc_matrix(&c, 3), gc_matrix(&r, 3));
    }
}

#[test]
fn a_section_gives_the_counts_along_a_sequence_with_its_statistics() {
    let a = Counts::count(&[b"AAAAAAAC"], 2).expect("counted");
    let s = section(&a, b"AAC");
    assert_eq!(s.coverage, vec![Some(6), Some(1)]);
    assert_eq!((s.median, s.mean, s.invalid, s.non_zero), (3.5, 3.5, 0, 2));
    assert!((s.gc - 100.0 / 3.0).abs() < 1e-9);
    // a letter that is no base: the windows with it have no number
    let n = section(&a, b"AANC");
    assert_eq!(n.coverage, vec![Some(6), None, None]);
    assert_eq!((n.invalid, n.non_zero), (2, 1));
    // a sequence that has nothing in common with the reads
    let none = section(&a, b"GGGG");
    assert_eq!(
        (none.coverage.clone(), none.median, none.non_zero),
        (vec![Some(0); 3], 0.0, 0)
    );
    // shorter than k: nothing to say
    assert!(section(&a, b"A").coverage.is_empty());
}

/// The Poisson probabilities for a mean of `lambda`, from 0 to `top`.
fn poisson(lambda: f64, top: usize) -> Vec<f64> {
    let mut p = vec![(-lambda).exp()];
    for m in 1..=top {
        p.push(p[m - 1] * lambda / m as f64);
    }
    p
}

#[test]
fn a_spectrum_of_errors_and_one_peak_gives_the_valley_the_peak_and_the_size() {
    // 1,000 distinct k-mers each read about 20 times, and errors that halve at every step from 1: the valley is where the two add up to least
    let genome = poisson(20.0, 60);
    let mut spec = vec![0u64; 61];
    for m in 1..=60 {
        spec[m] = (1000.0 * genome[m]).round() as u64 + (50_000.0 * 0.5f64.powi(m as i32)) as u64;
    }
    let s = summarize(&spec);
    let (valley, peak) = (s.valley.expect("a valley"), s.peak.expect("a peak"));
    assert!((6..=14).contains(&valley), "{valley}");
    assert!((19..=21).contains(&peak), "{peak}");
    let size = s.genome_size.expect("a size") as f64;
    assert!((size - 1000.0).abs() < 80.0, "{size}");
    assert!(s.error_share > 0.5 && s.half_peak.is_none(), "{s:?}");
}

#[test]
fn a_second_peak_at_half_the_first_is_found_and_a_spectrum_with_no_peak_has_no_size() {
    // most of the genome with one copy of each k-mer (peak at 40), some where the two copies differ (each k-mer at half that, 20)
    let (full, half) = (poisson(40.0, 100), poisson(20.0, 100));
    let mut spec = vec![0u64; 101];
    for m in 1..=100 {
        spec[m] = (900.0 * full[m] + 100.0 * half[m]).round() as u64
            + (50_000.0 * 0.5f64.powi(m as i32)) as u64;
    }
    let s = summarize(&spec);
    assert!((39..=41).contains(&s.peak.expect("a peak")), "{s:?}");
    assert!(s.half_peak.is_some_and(|h| (19..=21).contains(&h)), "{s:?}");
    // only errors, falling all the way: no valley, no peak, no size
    let falling: Vec<u64> = (0..30)
        .map(|m| {
            if m == 0 {
                0
            } else {
                100_000 / (m as u64 * m as u64)
            }
        })
        .collect();
    let s = summarize(&falling);
    assert_eq!((s.valley, s.peak, s.genome_size), (None, None, None));
    // nothing at all, and spectra too short to have a shape
    for tiny in [vec![], vec![0], vec![0, 5], vec![0, 5, 3], vec![0, 5, 3, 1]] {
        let s = summarize(&tiny);
        assert_eq!(s.genome_size, None, "{tiny:?}");
    }
}

#[test]
fn reads_of_a_simulated_genome_give_a_peak_at_the_coverage_of_the_k_mers_and_the_size_of_the_genome(
) {
    let mut rng = Rng(77);
    let genome = rng.bases(20_000);
    // 30x of 100 base reads from either strand, one base in 200 wrong: a k-mer of 21 is in a read with the coverage 30 * 80 / 100 = 24
    let reads: Vec<Vec<u8>> = (0..6000)
        .map(|_| {
            let at = rng.below(genome.len() - 100);
            let mut read = genome[at..at + 100].to_vec();
            for b in read.iter_mut() {
                if rng.below(200) == 0 {
                    *b = b"ACGT"[rng.below(4)];
                }
            }
            if rng.below(2) == 1 {
                read = revcom(&read);
            }
            read
        })
        .collect();
    let refs: Vec<&[u8]> = reads.iter().map(Vec::as_slice).collect();
    let c = Counts::count(&refs, 21).expect("counted");
    let s = summarize(&spectrum(&c, 1, 200));
    let peak = s.peak.expect("a peak");
    assert!((21..=27).contains(&peak), "{peak}");
    let size = s.genome_size.expect("a size") as f64;
    assert!((size - 20_000.0).abs() / 20_000.0 < 0.10, "{size}");
    assert!(s.error_share > 0.1, "{}", s.error_share);
}

/// A folder of its own for a test, removed when the test ends.
struct Folder(std::path::PathBuf);

impl Folder {
    fn new(name: &str) -> Folder {
        let path =
            std::env::temp_dir().join(format!("helixloom-kat-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&path).expect("a folder");
        Folder(path)
    }

    fn write(&self, name: &str, text: &str) -> String {
        let path = self.0.join(name);
        std::fs::write(&path, text).expect("a file");
        path.to_string_lossy().into_owned()
    }

    fn read(&self, name: &str) -> String {
        std::fs::read_to_string(self.0.join(name)).unwrap_or_else(|e| format!("({e})"))
    }

    fn run(&self, args: &[&str]) -> i32 {
        cli::run(&args.iter().map(|a| a.to_string()).collect::<Vec<_>>())
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).ok();
    }
}

#[test]
fn the_commands_write_the_spectrum_the_matrices_and_the_sections() {
    let f = Folder::new("commands");
    let a = f.write("a.fa", ">r1\nAAAAAAAC\n");
    let b = f.write("b.fq", "@r2\nACACAC\n+\nIIIIII\n");
    let out = f.0.join("out.txt").to_string_lossy().into_owned();
    // hist: the header says what was counted, then ### and the bins up to the last that is used
    assert_eq!(f.run(&["hist", "-k", "2", "-o", &out, &a]), 0);
    let hist = f.read("out.txt");
    assert!(
        hist.contains("# k: 2\n# distinct k-mers: 2\n# total k-mers: 7\n"),
        "{hist}"
    );
    let bins: Vec<&str> = hist.split("###\n").nth(1).expect("bins").lines().collect();
    assert_eq!((bins[0], bins[5]), ("1\t1", "6\t1"));
    assert_eq!(bins.len(), 6);
    // gcp: a row for each number of G and C
    assert_eq!(f.run(&["gcp", "-k", "2", "-o", &out, &a]), 0);
    let gcp = f.read("out.txt");
    let rows: Vec<&str> = gcp.split("###\n").nth(1).expect("rows").lines().collect();
    assert_eq!(rows[0], "gc\t1\t2\t3\t4\t5\t6");
    assert_eq!(rows[1], "0\t0\t0\t0\t0\t0\t1");
    assert_eq!(rows[2], "1\t1\t0\t0\t0\t0\t0");
    // comp: the k-mers of both, and of one only, in the header
    assert_eq!(f.run(&["comp", "-k", "2", "-o", &out, &a, &b]), 0);
    let comp = f.read("out.txt");
    assert!(
        comp.contains(
            "# distinct k-mers in both: 1\n# only in the first: 1\n# only in the second: 1\n"
        ),
        "{comp}"
    );
    // sect: statistics for each sequence, and the numbers along it
    let cov = f.0.join("cov.txt").to_string_lossy().into_owned();
    let seqs = f.write("s.fa", ">s1\nAAC\n>s2\nGGGG\n");
    assert_eq!(
        f.run(&[
            "sect",
            "-k",
            "2",
            "-c",
            &a,
            "--coverage",
            &cov,
            "-o",
            &out,
            &seqs
        ]),
        0
    );
    let stats = f.read("out.txt");
    let lines: Vec<&str> = stats.lines().collect();
    assert!(
        lines[0].starts_with("seq_name\tmedian\tmean\tgc%"),
        "{stats}"
    );
    assert!(
        lines[1].starts_with("s1\t3.5\t3.50\t33.33\t3\t2\t0\t0.00\t2\t100.00"),
        "{stats}"
    );
    assert!(
        lines[2].starts_with("s2\t0.0\t0.00\t100.00\t4\t3\t0\t0.00\t0\t0.00"),
        "{stats}"
    );
    assert_eq!(f.read("cov.txt"), ">s1\n6 1\n>s2\n0 0 0\n");
    // several files joined with commas are one input
    let joined = format!("{a},{b}");
    assert_eq!(f.run(&["hist", "-k", "2", "-o", &out, &joined]), 0);
    assert!(f.read("out.txt").contains("# total k-mers: 12\n"));
}

#[test]
fn a_command_that_cannot_be_run_says_so_and_stops() {
    let f = Folder::new("errors");
    let a = f.write("a.fa", ">r\nACGTACGT\n");
    let junk = f.write("junk.txt", "not a sequence file\n");
    for args in [
        vec!["hist"],
        vec!["nothing", a.as_str()],
        vec!["hist", "-k", "0", a.as_str()],
        vec!["hist", "-k", "33", a.as_str()],
        vec!["hist", "-k", "many", a.as_str()],
        vec!["hist", "--high", "0", a.as_str()],
        vec!["hist", "--unknown", a.as_str()],
        vec!["hist", junk.as_str()],
        vec!["hist", "missing.fa"],
        vec!["comp", a.as_str()],
        vec!["sect", a.as_str()],
        vec!["sect", "-c", a.as_str()],
    ] {
        assert_eq!(f.run(&args), 1, "{args:?}");
    }
    assert_eq!(f.run(&["--help"]), 0);
    assert_eq!(f.run(&[]), 0);
}

#[test]
fn a_k_mer_with_a_g_has_it_counted_among_the_g_and_c() {
    // AGGC: AG (one G), GG (CC as the canonical: two), GC (two)
    let c = Counts::count(&[b"AGGC"], 2).expect("counted");
    let g = gc_matrix(&c, 5);
    assert_eq!((g[0][1], g[1][1], g[2][1]), (0, 1, 2));
}

#[test]
fn the_valley_goes_to_the_end_of_a_flat_bottom_and_starts_the_count_from_it() {
    // errors fall away to a long flat bottom of 20 (bins 4 to 9), then the peak at 13
    let bins: [(usize, u64); 17] = [
        (1, 1000),
        (2, 400),
        (3, 100),
        (4, 20),
        (5, 20),
        (6, 20),
        (7, 20),
        (8, 20),
        (9, 20),
        (10, 60),
        (11, 200),
        (12, 500),
        (13, 800),
        (14, 500),
        (15, 200),
        (16, 60),
        (17, 10),
    ];
    let mut spec = vec![0u64; 18];
    for (m, n) in bins {
        spec[m] = n;
    }
    let s = summarize(&spec);
    assert!(s.valley.is_some_and(|v| (7..=9).contains(&v)), "{s:?}");
    assert_eq!(s.peak, Some(13));
}

#[test]
fn a_spectrum_that_rises_from_the_start_has_its_valley_at_two_and_its_size_is_worked_out_exactly() {
    // no errors to speak of: the counts rise from 1; the k-mers at 1 are not the genome's
    let spec = [0u64, 40, 40, 200, 500, 300, 100, 20];
    let s = summarize(&spec);
    assert_eq!((s.valley, s.peak), (Some(2), Some(4)));
    // (2*40 + 3*200 + 4*500 + 5*300 + 6*100 + 7*20) / 4
    assert_eq!(s.genome_size, Some(1230));
    assert!(
        (s.error_share - 40.0 / 1200.0).abs() < 1e-9,
        "{}",
        s.error_share
    );
}

#[test]
fn a_half_peak_must_be_a_tenth_as_high_as_the_peak_to_be_one() {
    let (full, half) = (poisson(40.0, 100), poisson(20.0, 100));
    let spectrum_with = |at_half: f64| -> Vec<u64> {
        let mut spec = vec![0u64; 101];
        for m in 1..=100 {
            spec[m] = (900.0 * full[m] + at_half * half[m]).round() as u64
                + (50_000.0 * 0.5f64.powi(m as i32)) as u64;
        }
        spec
    };
    // 100 k-mers at half: a bump of about 16 % of the peak; 30: about 5 %
    assert!(summarize(&spectrum_with(100.0)).half_peak.is_some());
    assert_eq!(summarize(&spectrum_with(30.0)).half_peak, None);
}

use super::*;
use crate::pipeline::asm::tests::fastq_reads;
use crate::seq::debruijn::sim::random;

fn quality_tail(seq: &[u8], good: usize) -> FastqRead {
    let mut qual = vec![40u8; good];
    qual.resize(seq.len(), 5);
    FastqRead {
        name: "r".into(),
        seq: seq.to_vec(),
        qual,
    }
}

fn kept(v: Value) -> Vec<FastqRead> {
    match v {
        Value::Reads(r) => r,
        other => panic!("{other:?}"),
    }
}

#[test]
fn the_adapter_is_cut_from_the_end_of_a_read() {
    let body = random(60, 1);
    let read = [&body[..], b"AGATCGGAAGAGC"].concat();
    let (out, note) = trim(&fastq_reads(&[read]), "AGATCGGAAGAGC", 20, 20, Cut::Off).unwrap();
    let r = kept(out);
    assert_eq!(r[0].seq, body);
    assert!(
        note.contains("1 of 1 reads kept") && note.contains("1 had the adapter cut"),
        "{note}"
    );
}

#[test]
fn with_none_for_the_adapter_it_stays() {
    let body = random(60, 2);
    let read = [&body[..], b"AGATCGGAAGAGC"].concat();
    let r = kept(
        trim(
            &fastq_reads(std::slice::from_ref(&read)),
            "none",
            20,
            20,
            Cut::Off,
        )
        .unwrap()
        .0,
    );
    assert_eq!(r[0].seq, read);
}

#[test]
fn a_tail_of_poor_quality_is_cut_off_with_cut_tail_and_right_and_kept_with_off() {
    let seq = random(100, 3);
    let read = quality_tail(&seq, 70);
    let tail = kept(
        trim(std::slice::from_ref(&read), "none", 20, 20, Cut::Tail)
            .unwrap()
            .0,
    );
    assert!(
        tail[0].seq.len() >= 66 && tail[0].seq.len() <= 72,
        "{}",
        tail[0].seq.len()
    );
    assert_eq!(tail[0].seq, seq[..tail[0].seq.len()]);
    let right = kept(
        trim(std::slice::from_ref(&read), "none", 20, 20, Cut::Right)
            .unwrap()
            .0,
    );
    assert!(
        right[0].seq.len() < 100 && right[0].seq.len() >= 60,
        "{}",
        right[0].seq.len()
    );
    let off = kept(trim(&[read], "none", 20, 20, Cut::Off).unwrap().0);
    assert_eq!(off[0].seq.len(), 100);
}

#[test]
fn the_quality_asked_for_decides_where_the_cut_is() {
    let seq = random(100, 4);
    let mut qual = vec![40u8; 100];
    qual[60..].iter_mut().for_each(|q| *q = 25);
    let read = FastqRead {
        name: "r".into(),
        seq,
        qual,
    };
    let lax = kept(
        trim(std::slice::from_ref(&read), "none", 20, 20, Cut::Tail)
            .unwrap()
            .0,
    );
    let strict = kept(trim(&[read], "none", 30, 20, Cut::Tail).unwrap().0);
    assert_eq!(lax[0].seq.len(), 100);
    assert!(strict[0].seq.len() < 70, "{}", strict[0].seq.len());
}

#[test]
fn a_read_that_is_too_short_is_dropped_and_when_none_is_left_the_step_fails_and_says_why() {
    let reads = fastq_reads(&[random(30, 5), random(80, 6)]);
    let (out, note) = trim(&reads, "none", 20, 50, Cut::Off).unwrap();
    assert_eq!(kept(out).len(), 1);
    assert!(
        note.contains("1 of 2 reads kept") && note.contains("1 too short"),
        "{note}"
    );
    let e = trim(&reads, "none", 20, 500, Cut::Off).unwrap_err();
    assert!(
        e.contains("no read is left of 2 after trimming") && e.contains("2 too short"),
        "{e}"
    );
}

#[test]
fn the_long_reads_are_selected_as_the_paper_does() {
    let mut reads = fastq_reads(&[
        random(1500, 7),
        random(3000, 8),
        random(3000, 9),
        random(3000, 10),
    ]);
    reads[2].qual = vec![10; 3000];
    let (out, note) = filter(&reads, 2000, 6000, 100).unwrap();
    let r = kept(out);
    assert_eq!(
        r.iter().map(|x| x.name.as_str()).collect::<Vec<_>>(),
        ["r1", "r3"]
    );
    assert!(
        note.contains("2 of 4 reads kept")
            && note.contains("10500 to 6000 bases")
            && note.contains("1 too short")
            && note.contains("1 dropped for quality"),
        "{note}"
    );
}

#[test]
fn when_no_read_is_long_enough_the_step_fails() {
    let e = filter(&fastq_reads(&[random(100, 11)]), 2000, 500_000_000, 100).unwrap_err();
    assert!(e.contains("no read of 1 is 2000 bases or longer"), "{e}");
}

#[test]
fn a_window_whose_mean_is_the_quality_asked_for_is_not_cut() {
    // quality 30 all through and 30 asked for: nothing is cut; 31 is asked for, and the whole read goes
    let seq = random(100, 12);
    let read = FastqRead {
        name: "r".into(),
        seq,
        qual: vec![30; 100],
    };
    let same = kept(
        trim(std::slice::from_ref(&read), "none", 30, 20, Cut::Tail)
            .unwrap()
            .0,
    );
    assert_eq!(same[0].seq.len(), 100);
    let over = trim(std::slice::from_ref(&read), "none", 31, 20, Cut::Tail);
    assert!(over.is_err(), "{over:?}");
}

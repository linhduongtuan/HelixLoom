use super::*;

/// 10 bases, then the exons and the intron of a gene, then 9 bases: the mRNA is ATGAAACCCGGGTTT + AAATAG, the protein MKPGFK.
const FLANK_LEFT: &str = "CCCCCCCCCC";
const EXON_1: &str = "ATGAAACCCGGGTTT";
const INTRON: &str = "GTCCCCCCAG";
const EXON_2: &str = "AAATAG";
const FLANK_RIGHT: &str = "CCCCCCCCC";

fn plus_strand_chromosome() -> String {
    format!("{FLANK_LEFT}{EXON_1}{INTRON}{EXON_2}{FLANK_RIGHT}")
}

fn genome_of(parts: &[(&str, String)]) -> Genome {
    let text: String = parts
        .iter()
        .map(|(name, bases)| format!(">{name} some description\n{bases}\n"))
        .collect();
    Genome::from_fasta(&text)
}

/// The gene of the constants on the plus strand as GFF3: the exons are 11-25 and 36-41, the CDS the same.
const GFF3_PLUS: &str = "##gff-version 3
chr1\ttest\tgene\t11\t41\t.\t+\t.\tID=gene1;Name=alpha
chr1\ttest\tmRNA\t11\t41\t.\t+\t.\tID=rna1;Parent=gene1
chr1\ttest\texon\t11\t25\t.\t+\t.\tParent=rna1
chr1\ttest\texon\t36\t41\t.\t+\t.\tParent=rna1
chr1\ttest\tCDS\t11\t25\t.\t+\t0\tParent=rna1
chr1\ttest\tCDS\t36\t41\t.\t+\t0\tParent=rna1
";

#[test]
fn a_gene_of_gff3_is_a_transcript_with_its_exons_its_cds_and_the_gene_it_belongs_to() {
    let annotation = parse(GFF3_PLUS).expect("parsed");
    assert_eq!(annotation.transcripts.len(), 1);
    let t = &annotation.transcripts[0];
    assert_eq!(
        (t.id.as_str(), t.kind.as_str(), t.seqname.as_str(), t.strand),
        ("rna1", "mRNA", "chr1", Strand::Plus)
    );
    assert_eq!((t.start, t.end), (11, 41));
    assert_eq!(t.exons, vec![(11, 25), (36, 41)]);
    assert_eq!(
        t.cds,
        vec![
            Cds {
                start: 11,
                end: 25,
                phase: 0
            },
            Cds {
                start: 36,
                end: 41,
                phase: 0
            }
        ]
    );
    assert_eq!(t.gene_id.as_deref(), Some("gene1"));
    assert_eq!(t.gene_name.as_deref(), Some("alpha"));
    assert!(t.is_coding());
}

#[test]
fn the_spliced_transcript_the_cds_and_the_protein_are_cut_from_the_genome_as_the_annotation_says() {
    let genome = genome_of(&[("chr1", plus_strand_chromosome())]);
    let t = &parse(GFF3_PLUS).expect("parsed").transcripts[0];
    assert_eq!(spliced(t, &genome), Ok(b"ATGAAACCCGGGTTTAAATAG".to_vec()));
    assert_eq!(
        cds_sequence(t, &genome),
        Ok(Some(b"ATGAAACCCGGGTTTAAATAG".to_vec()))
    );
    // the stop codon at the end is left out of the protein, and a stop in the middle is shown with the symbol asked for
    assert_eq!(protein(t, &genome, b'.'), Ok(Some(b"MKPGFK".to_vec())));
}

#[test]
fn a_gene_on_the_minus_strand_gives_the_same_transcript_and_protein() {
    // the same gene written on the other strand: the chromosome is the reverse complement of the first, and its coordinates are mirrored
    let plus = plus_strand_chromosome();
    let length = plus.len();
    let chromosome = String::from_utf8(crate::seq::dna::revcom(plus.as_bytes())).expect("text");
    let mirrored = |start: usize, end: usize| (length + 1 - end, length + 1 - start);
    let (e1, e2) = (mirrored(11, 25), mirrored(36, 41));
    let gff = format!(
        "##gff-version 3
chr1\t.\tmRNA\t{}\t{}\t.\t-\t.\tID=rna1
chr1\t.\texon\t{}\t{}\t.\t-\t.\tParent=rna1
chr1\t.\texon\t{}\t{}\t.\t-\t.\tParent=rna1
chr1\t.\tCDS\t{}\t{}\t.\t-\t0\tParent=rna1
chr1\t.\tCDS\t{}\t{}\t.\t-\t0\tParent=rna1
",
        e2.0, e1.1, e1.0, e1.1, e2.0, e2.1, e1.0, e1.1, e2.0, e2.1
    );
    let t = &parse(&gff).expect("parsed").transcripts[0];
    assert_eq!(t.strand, Strand::Minus);
    // the exons are listed from the left, whichever the strand
    assert_eq!(t.exons, vec![e2, e1]);
    let genome = genome_of(&[("chr1", chromosome)]);
    assert_eq!(spliced(t, &genome), Ok(b"ATGAAACCCGGGTTTAAATAG".to_vec()));
    assert_eq!(protein(t, &genome, b'*'), Ok(Some(b"MKPGFK".to_vec())));
}

#[test]
fn a_cds_that_begins_in_the_middle_of_a_codon_loses_the_bases_its_phase_says() {
    // 1 base before the ATG belongs to a codon that is not in this CDS: phase 1 skips it
    let genome = genome_of(&[("c", "TATGAAACCCTAACCC".to_string())]);
    let gff = "##gff-version 3\nc\t.\tmRNA\t1\t13\t.\t+\t.\tID=t\nc\t.\texon\t1\t13\t.\t+\t.\tParent=t\nc\t.\tCDS\t1\t13\t.\t+\t1\tParent=t\n";
    let t = &parse(gff).expect("parsed").transcripts[0];
    assert_eq!(cds_sequence(t, &genome), Ok(Some(b"ATGAAACCCTAA".to_vec())));
    assert_eq!(protein(t, &genome, b'.'), Ok(Some(b"MKP".to_vec())));
    // and phase 2 skips two
    let gff2 = gff.replace("\t+\t1\tParent", "\t+\t2\tParent");
    let t2 = &parse(&gff2).expect("parsed").transcripts[0];
    assert_eq!(cds_sequence(t2, &genome), Ok(Some(b"TGAAACCCTAA".to_vec())));
}

#[test]
fn a_stop_in_the_middle_is_shown_and_a_transcript_without_cds_has_no_protein() {
    let genome = genome_of(&[("c", "ATGTAAAAATAG".to_string())]);
    let gff = "##gff-version 3\nc\t.\tmRNA\t1\t12\t.\t+\t.\tID=t\nc\t.\texon\t1\t12\t.\t+\t.\tParent=t\nc\t.\tCDS\t1\t12\t.\t+\t0\tParent=t\n\
               c\t.\tncRNA\t1\t12\t.\t+\t.\tID=n\nc\t.\texon\t1\t12\t.\t+\t.\tParent=n\n";
    let annotation = parse(gff).expect("parsed");
    let (t, n) = (&annotation.transcripts[0], &annotation.transcripts[1]);
    assert_eq!(protein(t, &genome, b'.'), Ok(Some(b"M.K".to_vec())));
    assert_eq!(protein(t, &genome, b'*'), Ok(Some(b"M*K".to_vec())));
    assert!(!n.is_coding());
    assert_eq!(protein(n, &genome, b'.'), Ok(None));
    assert_eq!(cds_sequence(n, &genome), Ok(None));
    assert_eq!(spliced(n, &genome), Ok(b"ATGTAAAAATAG".to_vec()));
}

const GTF_PLUS: &str = "chr1\ttest\texon\t11\t25\t.\t+\t.\tgene_id \"gene1\"; transcript_id \"rna1\"; gene_name \"alpha\";
chr1\ttest\texon\t36\t41\t.\t+\t.\tgene_id \"gene1\"; transcript_id \"rna1\"; gene_name \"alpha\";
chr1\ttest\tCDS\t11\t25\t.\t+\t0\tgene_id \"gene1\"; transcript_id \"rna1\";
chr1\ttest\tCDS\t36\t41\t.\t+\t0\tgene_id \"gene1\"; transcript_id \"rna1\";
chr1\ttest\tstart_codon\t11\t13\t.\t+\t0\tgene_id \"gene1\"; transcript_id \"rna1\";
chr1\ttest\tstop_codon\t39\t41\t.\t+\t0\tgene_id \"gene1\"; transcript_id \"rna1\";
";

#[test]
fn gtf_gives_the_same_transcript_and_start_and_stop_codon_lines_are_ignored_when_there_is_cds() {
    let t = &parse(GTF_PLUS).expect("parsed").transcripts[0];
    assert_eq!(
        (t.id.as_str(), t.seqname.as_str(), t.strand),
        ("rna1", "chr1", Strand::Plus)
    );
    assert_eq!(t.exons, vec![(11, 25), (36, 41)]);
    assert_eq!(
        t.cds.len(),
        2,
        "the codons are not CDS segments of their own"
    );
    assert_eq!(t.gene_id.as_deref(), Some("gene1"));
    assert_eq!(t.gene_name.as_deref(), Some("alpha"));
    let genome = genome_of(&[("chr1", plus_strand_chromosome())]);
    assert_eq!(protein(t, &genome, b'.'), Ok(Some(b"MKPGFK".to_vec())));
}

#[test]
fn utr_lines_are_ignored_when_there_are_exons_and_exons_are_made_from_cds_and_utr_when_there_are_none(
) {
    // exons given: the UTR lines add nothing
    let with_exons = format!(
        "{GTF_PLUS}chr1\ttest\tUTR\t1\t10\t.\t+\t.\tgene_id \"gene1\"; transcript_id \"rna1\";\n"
    );
    assert_eq!(
        parse(&with_exons).expect("parsed").transcripts[0].exons,
        vec![(11, 25), (36, 41)]
    );
    // no exons: the CDS and the UTRs make them, and segments that touch are one exon
    let gtf = "c\t.\t5UTR\t1\t10\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";
c\t.\tCDS\t11\t25\t.\t+\t0\tgene_id \"g\"; transcript_id \"t\";
c\t.\tCDS\t36\t41\t.\t+\t0\tgene_id \"g\"; transcript_id \"t\";
c\t.\t3UTR\t42\t50\t.\t+\t.\tgene_id \"g\"; transcript_id \"t\";
";
    let t = &parse(gtf).expect("parsed").transcripts[0];
    assert_eq!(t.exons, vec![(1, 25), (36, 50)]);
    // only CDS: the CDS are the exons
    let only = "c\t.\tCDS\t11\t25\t.\t+\t0\tgene_id \"g\"; transcript_id \"t\";\n";
    assert_eq!(
        parse(only).expect("parsed").transcripts[0].exons,
        vec![(11, 25)]
    );
}

#[test]
fn what_is_written_is_read_back_the_same_as_gff3_and_as_gtf() {
    let annotation = parse(GFF3_PLUS).expect("parsed");
    let again = parse(&write_gff3(&annotation, false)).expect("parsed");
    assert_eq!(again.transcripts, annotation.transcripts);
    let gtf = write_gtf(&annotation, false);
    assert!(
        gtf.contains("transcript_id \"rna1\"") && gtf.contains("gene_id \"gene1\""),
        "{gtf}"
    );
    let back = parse(&gtf).expect("parsed");
    let (a, b) = (&annotation.transcripts[0], &back.transcripts[0]);
    assert_eq!(
        (&a.id, &a.exons, &a.cds, a.strand),
        (&b.id, &b.exons, &b.cds, b.strand)
    );
    assert_eq!((&a.gene_id, &a.gene_name), (&b.gene_id, &b.gene_name));
}

#[test]
fn names_with_the_characters_gff3_escapes_come_back_whole() {
    let gff = "##gff-version 3\nc\t.\tmRNA\t1\t9\t.\t+\t.\tID=a%3Bb%3Dc;Name=x%2Cy\nc\t.\texon\t1\t9\t.\t+\t.\tParent=a%3Bb%3Dc\n";
    let t = &parse(gff).expect("parsed").transcripts[0];
    assert_eq!(t.id, "a;b=c");
    let written = write_gff3(&parse(gff).expect("parsed"), false);
    assert!(written.contains("ID=a%3Bb%3Dc"), "{written}");
    assert_eq!(parse(&written).expect("parsed").transcripts[0].id, "a;b=c");
}

#[test]
fn a_cds_without_a_parent_is_a_transcript_of_its_own_as_bacterial_annotations_write_them_and_features_in_several_lines_share_an_id(
) {
    let gff = "##gff-version 3
p\tProdigal\tCDS\t101\t400\t.\t+\t0\tID=cds1;product=thing
p\tProdigal\tCDS\t501\t700\t.\t-\t0\tID=cds2
p\t.\tmRNA\t1000\t1100\t.\t+\t.\tID=m
p\t.\tCDS\t1000\t1040\t.\t+\t0\tID=shared;Parent=m
p\t.\tCDS\t1061\t1100\t.\t+\t1\tID=shared;Parent=m
";
    let a = parse(gff).expect("parsed");
    assert_eq!(a.transcripts.len(), 3);
    let first = &a.transcripts[0];
    assert_eq!(
        (first.id.as_str(), first.exons.clone(), first.cds.len()),
        ("cds1", vec![(101, 400)], 1)
    );
    assert_eq!(a.transcripts[1].strand, Strand::Minus);
    assert_eq!(a.transcripts[2].cds.len(), 2);
}

#[test]
fn an_exon_that_has_two_parents_belongs_to_both() {
    let gff = "##gff-version 3
c\t.\tmRNA\t1\t30\t.\t+\t.\tID=a
c\t.\tmRNA\t1\t30\t.\t+\t.\tID=b
c\t.\texon\t1\t10\t.\t+\t.\tParent=a,b
c\t.\texon\t21\t30\t.\t+\t.\tParent=a
";
    let t = parse(gff).expect("parsed").transcripts;
    assert_eq!((t[0].exons.len(), t[1].exons.len()), (2, 1));
}

#[test]
fn a_file_that_cannot_be_read_is_refused_with_its_line() {
    for (text, why) in [
        ("", "no annotation"),
        ("##gff-version 3\n", "no annotation"),
        ("c\t.\texon\t1\t9\t.\t+\n", "Line 1"),
        // eight columns: the attributes are missing
        ("c\t.\texon\t1\t9\t.\t+\t.\n", "Line 1"),
        (
            "##gff-version 3\nc\t.\texon\tx\t9\t.\t+\t.\tParent=a\n",
            "Line 2",
        ),
        (
            "##gff-version 3\nc\t.\texon\t9\t1\t.\t+\t.\tParent=a\n",
            "start",
        ),
        ("c\t.\texon\t1\t9\t.\t+\t.\tnothing here\n", "GFF3 or GTF"),
    ] {
        let said = parse(text).expect_err(text);
        assert!(said.contains(why), "{text:?}: {said}");
    }
}

#[test]
fn a_transcript_cannot_be_cut_from_a_genome_that_lacks_its_sequence_or_is_too_short() {
    let t = &parse(GFF3_PLUS).expect("parsed").transcripts[0];
    let other = genome_of(&[("chr2", plus_strand_chromosome())]);
    assert!(spliced(t, &other).expect_err("no chr1").contains("chr1"));
    let short = genome_of(&[("chr1", "ACGT".to_string())]);
    assert!(spliced(t, &short)
        .expect_err("too short")
        .contains("beyond"));
}

fn many() -> Annotation {
    let gff = "##gff-version 3
c\t.\tmRNA\t11\t41\t.\t+\t.\tID=coding
c\t.\texon\t11\t25\t.\t+\t.\tParent=coding
c\t.\texon\t36\t41\t.\t+\t.\tParent=coding
c\t.\tCDS\t11\t25\t.\t+\t0\tParent=coding
c\t.\tCDS\t36\t41\t.\t+\t0\tParent=coding
c\t.\tncRNA\t11\t20\t.\t+\t.\tID=noncoding
c\t.\texon\t11\t20\t.\t+\t.\tParent=noncoding
c\t.\tmRNA\t60\t90\t.\t-\t.\tID=other
c\t.\texon\t60\t70\t.\t-\t.\tParent=other
c\t.\texon\t80\t90\t.\t-\t.\tParent=other
";
    parse(gff).expect("parsed")
}

fn kept(a: &Annotation, f: &Filters, genome: Option<&Genome>) -> Vec<String> {
    a.transcripts
        .iter()
        .filter(|t| f.keeps(t, genome).expect("judged"))
        .map(|t| t.id.clone())
        .collect()
}

#[test]
fn the_filters_keep_what_they_say() {
    let a = many();
    let f = |f: Filters| kept(&a, &f, None);
    assert_eq!(f(Filters::default()), ["coding", "noncoding", "other"]);
    assert_eq!(
        f(Filters {
            coding_only: true,
            ..Filters::default()
        }),
        ["coding"]
    );
    assert_eq!(
        f(Filters {
            noncoding_only: true,
            ..Filters::default()
        }),
        ["noncoding", "other"]
    );
    assert_eq!(
        f(Filters {
            no_single_exon: true,
            ..Filters::default()
        }),
        ["coding", "other"]
    );
    // introns: coding has one of 10, other has one of 9
    assert_eq!(
        f(Filters {
            max_intron: Some(9),
            ..Filters::default()
        }),
        ["noncoding", "other"]
    );
    assert_eq!(
        f(Filters {
            min_intron: Some(10),
            ..Filters::default()
        }),
        ["coding", "noncoding"]
    );
    // a range: overlapping 1-20 keeps those that touch it, and contained keeps those inside it
    let overlapping = Range {
        seqname: "c".into(),
        start: 1,
        end: 20,
        strand: None,
    };
    assert_eq!(
        f(Filters {
            range: Some((overlapping.clone(), false)),
            ..Filters::default()
        }),
        ["coding", "noncoding"]
    );
    assert_eq!(
        f(Filters {
            range: Some((overlapping, true)),
            ..Filters::default()
        }),
        ["noncoding"]
    );
    let on_minus = Range {
        seqname: "c".into(),
        start: 1,
        end: 100,
        strand: Some(Strand::Minus),
    };
    assert_eq!(
        f(Filters {
            range: Some((on_minus, false)),
            ..Filters::default()
        }),
        ["other"]
    );
    // lists of ids
    let ids: std::collections::HashSet<String> = ["coding".to_string()].into();
    assert_eq!(
        f(Filters {
            only_ids: Some(ids.clone()),
            ..Filters::default()
        }),
        ["coding"]
    );
    assert_eq!(
        f(Filters {
            not_ids: Some(ids),
            ..Filters::default()
        }),
        ["noncoding", "other"]
    );
}

#[test]
fn the_filters_that_read_the_genome_judge_splice_sites_and_stops() {
    // GT...AG introns are the usual ones; an intron that begins CC does not follow the rule
    let good = plus_strand_chromosome();
    let bad = good.replace("GTCCCCCCAG", "CCCCCCCCAG");
    let a = parse(GFF3_PLUS).expect("parsed");
    let canonical = Filters {
        canonical_splice: true,
        ..Filters::default()
    };
    assert_eq!(
        kept(&a, &canonical, Some(&genome_of(&[("chr1", good.clone())]))),
        ["rna1"]
    );
    assert!(kept(&a, &canonical, Some(&genome_of(&[("chr1", bad)]))).is_empty());
    // a stop in the middle of the CDS
    let in_frame = good.replace("AAACCC", "TAACCC");
    let no_stop = Filters {
        no_inframe_stop: true,
        ..Filters::default()
    };
    assert_eq!(
        kept(&a, &no_stop, Some(&genome_of(&[("chr1", good.clone())]))),
        ["rna1"]
    );
    assert!(kept(&a, &no_stop, Some(&genome_of(&[("chr1", in_frame)]))).is_empty());
    // a complete CDS: starts with a start codon, ends with a stop, no stop inside
    let complete = Filters {
        complete_cds: true,
        ..Filters::default()
    };
    assert_eq!(
        kept(&a, &complete, Some(&genome_of(&[("chr1", good.clone())]))),
        ["rna1"]
    );
    let no_start = good.replace("ATGAAA", "CTGAAA");
    assert!(kept(&a, &complete, Some(&genome_of(&[("chr1", no_start)]))).is_empty());
    // these need the genome, and say so
    assert!(canonical
        .keeps(&a.transcripts[0], None)
        .expect_err("needs a genome")
        .contains("genome"));
}

#[test]
fn a_range_is_read_from_its_text() {
    let r = Range::parse("chr1:100-2000").expect("a range");
    assert_eq!(
        (r.seqname.as_str(), r.start, r.end, r.strand),
        ("chr1", 100, 2000, None)
    );
    let r = Range::parse("-chr1:5-6").expect("a range");
    assert_eq!((r.strand, r.start), (Some(Strand::Minus), 5));
    assert!(Range::parse("chr1").is_err() && Range::parse("chr1:9-3").is_err());
}

#[test]
fn the_order_of_the_output_is_by_sequence_as_met_and_then_by_position() {
    let gff = "##gff-version 3
b\t.\tmRNA\t50\t60\t.\t+\t.\tID=b2
b\t.\texon\t50\t60\t.\t+\t.\tParent=b2
a\t.\tmRNA\t80\t90\t.\t+\t.\tID=a2
a\t.\texon\t80\t90\t.\t+\t.\tParent=a2
a\t.\tmRNA\t10\t20\t.\t+\t.\tID=a1
a\t.\texon\t10\t20\t.\t+\t.\tParent=a1
";
    let mut a = parse(gff).expect("parsed");
    a.sort_by_location();
    let ids: Vec<&str> = a.transcripts.iter().map(|t| t.id.as_str()).collect();
    assert_eq!(ids, ["b2", "a1", "a2"]);
}

#[test]
fn the_fasta_of_a_genome_is_found_by_the_first_word_of_its_header() {
    let g = genome_of(&[("chr1", "ACGT".to_string())]);
    assert_eq!(g.sequence("chr1"), Some(&b"ACGT"[..]));
    assert_eq!(g.sequence("chr2"), None);
}

/// A folder of its own for a test, with the genome and the annotation in it, removed when the test ends.
struct Folder(std::path::PathBuf);

impl Folder {
    fn new(name: &str) -> Folder {
        let path =
            std::env::temp_dir().join(format!("helixloom-gff-{}-{name}", std::process::id()));
        std::fs::create_dir_all(&path).expect("a folder");
        std::fs::write(
            path.join("genome.fa"),
            format!(">chr1 a description\n{}\n", plus_strand_chromosome()),
        )
        .expect("a file");
        std::fs::write(path.join("a.gff3"), GFF3_PLUS).expect("a file");
        Folder(path)
    }

    fn p(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }

    fn text(&self, name: &str) -> String {
        std::fs::read_to_string(self.0.join(name)).unwrap_or_else(|e| format!("({e})"))
    }

    fn run(&self, args: &[&str]) -> i32 {
        let args: Vec<String> = args.iter().map(|a| a.to_string()).collect();
        cli::run(&args)
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).ok();
    }
}

#[test]
fn the_command_writes_the_sequences_that_were_asked_for_and_not_the_annotation() {
    let f = Folder::new("sequences");
    let code = f.run(&[
        &f.p("a.gff3"),
        "-g",
        &f.p("genome.fa"),
        "-w",
        &f.p("w.fa"),
        "-x",
        &f.p("x.fa"),
        "-y",
        &f.p("y.fa"),
    ]);
    assert_eq!(code, 0);
    assert_eq!(f.text("w.fa"), ">rna1\nATGAAACCCGGGTTTAAATAG\n");
    assert_eq!(f.text("x.fa"), ">rna1\nATGAAACCCGGGTTTAAATAG\n");
    assert_eq!(f.text("y.fa"), ">rna1\nMKPGFK\n");
    // with the coordinates in the names
    let code = f.run(&[
        &f.p("a.gff3"),
        "-g",
        &f.p("genome.fa"),
        "-W",
        "-w",
        &f.p("w2.fa"),
        "-x",
        &f.p("x2.fa"),
    ]);
    assert_eq!(code, 0);
    assert_eq!(
        f.text("w2.fa").lines().next(),
        Some(">rna1 loc:chr1(+)11-41 exons:11-25,36-41")
    );
    assert_eq!(
        f.text("x2.fa").lines().next(),
        Some(">rna1 loc:chr1(+)11-41 segs:1-15,16-21")
    );
    // a stop inside is a * when -S says so
    let stop = Folder::new("stop");
    std::fs::write(stop.0.join("genome.fa"), ">chr1\nATGTAAAAATAG\n").expect("a file");
    std::fs::write(stop.0.join("a.gff3"), "##gff-version 3\nchr1\t.\tmRNA\t1\t12\t.\t+\t.\tID=t\nchr1\t.\tCDS\t1\t12\t.\t+\t0\tParent=t\n").expect("a file");
    stop.run(&[
        &stop.p("a.gff3"),
        "-g",
        &stop.p("genome.fa"),
        "-y",
        &stop.p("dot.fa"),
    ]);
    stop.run(&[
        &stop.p("a.gff3"),
        "-g",
        &stop.p("genome.fa"),
        "-S",
        "-y",
        &stop.p("star.fa"),
    ]);
    assert_eq!(
        (stop.text("dot.fa"), stop.text("star.fa")),
        (">t\nM.K\n".to_string(), ">t\nM*K\n".to_string())
    );
}

#[test]
fn the_command_writes_gff3_or_gtf_of_what_passes_the_filters() {
    let f = Folder::new("write");
    let mixed = format!("{GFF3_PLUS}chr1\t.\tncRNA\t50\t60\t.\t+\t.\tID=nc\nchr1\t.\texon\t50\t60\t.\t+\t.\tParent=nc\n");
    std::fs::write(f.0.join("m.gff3"), mixed).expect("a file");
    assert_eq!(f.run(&[&f.p("m.gff3"), "-C", "-o", &f.p("coding.gff3")]), 0);
    let written = f.text("coding.gff3");
    assert!(
        written.contains("ID=rna1") && !written.contains("ID=nc"),
        "{written}"
    );
    assert_eq!(
        f.run(&[&f.p("m.gff3"), "--nc", "-T", "-o", &f.p("nc.gtf")]),
        0
    );
    let gtf = f.text("nc.gtf");
    assert!(
        gtf.contains("transcript_id \"nc\"") && !gtf.contains("rna1"),
        "{gtf}"
    );
    // the ids of a file, kept and left out
    std::fs::write(f.0.join("ids.txt"), "nc\n").expect("a file");
    assert_eq!(
        f.run(&[
            &f.p("m.gff3"),
            "--ids",
            &f.p("ids.txt"),
            "-o",
            &f.p("only.gff3")
        ]),
        0
    );
    assert!(f.text("only.gff3").contains("ID=nc") && !f.text("only.gff3").contains("rna1"));
    assert_eq!(
        f.run(&[
            &f.p("m.gff3"),
            "--nids",
            &f.p("ids.txt"),
            "-o",
            &f.p("not.gff3")
        ]),
        0
    );
    assert!(f.text("not.gff3").contains("rna1") && !f.text("not.gff3").contains("ID=nc"));
    // a range
    assert_eq!(
        f.run(&[&f.p("m.gff3"), "-r", "chr1:45-100", "-o", &f.p("r.gff3")]),
        0
    );
    assert!(f.text("r.gff3").contains("ID=nc") && !f.text("r.gff3").contains("rna1"));
}

#[test]
fn the_table_has_a_line_for_each_transcript_with_the_fields_asked_for() {
    let f = Folder::new("table");
    let code = f.run(&[
        &f.p("a.gff3"),
        "--table",
        "@id,@chr,@start,@end,@strand,@numexons,@exons,@cds,@covlen,@cdslen,@geneid,@genename",
        "-o",
        &f.p("t.tsv"),
    ]);
    assert_eq!(code, 0);
    assert_eq!(
        f.text("t.tsv"),
        "rna1\tchr1\t11\t41\t+\t2\t11-25,36-41\t11-25,36-41\t21\t21\tgene1\talpha\n"
    );
}

#[test]
fn a_command_that_cannot_be_run_says_so_and_stops() {
    let f = Folder::new("errors");
    let bad = |args: &[&str]| f.run(args);
    assert_eq!(bad(&[]), 1);
    assert_eq!(bad(&[&f.p("a.gff3"), "--nothing"]), 1);
    assert_eq!(
        bad(&[&f.p("a.gff3"), "-w", &f.p("w.fa")]),
        1,
        "-w needs the genome"
    );
    assert_eq!(bad(&[&f.p("a.gff3"), "-R"]), 1, "-R goes with -r");
    assert_eq!(bad(&[&f.p("a.gff3"), "-N"]), 1, "-N reads the sequence");
    assert_eq!(bad(&[&f.p("missing.gff3")]), 1);
    assert_eq!(bad(&[&f.p("a.gff3"), "-i", "many"]), 1);
    assert_eq!(f.run(&["--help"]), 0);
    // the sequence of the annotation is not in the genome
    std::fs::write(f.0.join("other.fa"), ">chrX\nACGT\n").expect("a file");
    assert_eq!(
        bad(&[&f.p("a.gff3"), "-g", &f.p("other.fa"), "-w", &f.p("w.fa")]),
        1
    );
}

/// A plasmid's annotation as a bacterial annotation program writes it: a gene and its CDS on each strand, a tRNA, and one feature that
/// is cut off by the end of the sequence.
const PLASMID_GFF: &str = "##gff-version 3
pGM\tProkka\tgene\t101\t400\t.\t+\t.\tID=g1;Name=repA
pGM\tProkka\tCDS\t101\t400\t.\t+\t0\tID=c1;Parent=g1;Name=repA;product=replication protein
pGM\tProkka\tCDS\t501\t700\t.\t-\t0\tID=c2;product=hypothetical protein
pGM\tProkka\ttRNA\t801\t875\t.\t+\t.\tID=t1;product=tRNA-Gly
pGM\tProkka\tCDS\t950\t1100\t.\t+\t0\tID=c3
";

#[test]
fn an_annotation_becomes_features_of_the_sequence_with_the_strand_the_gene_and_the_product() {
    let a = parse(PLASMID_GFF).expect("parsed");
    let done = features_for(&a, "pGM", 1000).expect("features");
    assert_eq!(done.seqname, "pGM");
    // the CDS that runs past 1000 is left out, and said
    assert_eq!(done.skipped, 1);
    let names: Vec<(&str, &str, bool)> = done
        .features
        .iter()
        .map(|f| (f.name.as_str(), f.kind.as_str(), f.reverse))
        .collect();
    assert_eq!(
        names,
        [
            ("repA", "CDS", false),
            ("c2", "CDS", true),
            ("t1", "tRNA", false)
        ]
    );
    // positions are counted from 0 and end where the feature does
    assert_eq!(done.features[0].spans, vec![(100, 400)]);
    let q = |f: &crate::seq::features::Feature, k: &str| {
        f.qualifiers
            .iter()
            .find(|(key, _)| key == k)
            .map(|(_, v)| v.clone())
    };
    assert_eq!(
        q(&done.features[0], "product").as_deref(),
        Some("replication protein")
    );
    assert_eq!(q(&done.features[0], "gene").as_deref(), Some("repA"));
    assert_eq!(
        q(&done.features[1], "product").as_deref(),
        Some("hypothetical protein")
    );
}

#[test]
fn a_cds_with_exons_beyond_it_gives_an_mrna_too_and_a_spliced_cds_keeps_its_segments() {
    let a = parse(GFF3_PLUS).expect("parsed");
    let done = features_for(&a, "chr1", 100).expect("features");
    let kinds: Vec<(&str, Vec<(usize, usize)>)> = done
        .features
        .iter()
        .map(|f| (f.kind.as_str(), f.spans.clone()))
        .collect();
    // the exons and the CDS are the same stretches here, so there is the CDS only
    assert_eq!(kinds, [("CDS", vec![(10, 25), (35, 41)])]);
    let with_utr = "##gff-version 3\nc\t.\tmRNA\t1\t60\t.\t+\t.\tID=m\nc\t.\texon\t1\t60\t.\t+\t.\tParent=m\nc\t.\tCDS\t11\t40\t.\t+\t0\tParent=m\n";
    let done = features_for(&parse(with_utr).expect("parsed"), "c", 100).expect("features");
    let kinds: Vec<&str> = done.features.iter().map(|f| f.kind.as_str()).collect();
    assert_eq!(kinds, ["mRNA", "CDS"]);
}

#[test]
fn the_sequence_of_an_annotation_is_found_by_name_by_the_first_word_of_the_name_or_as_the_only_one()
{
    let a = parse(PLASMID_GFF).expect("parsed");
    assert!(features_for(&a, "pGM some description", 2000).is_ok());
    // another name, but the annotation has one sequence only: that is the one
    assert_eq!(
        features_for(&a, "my plasmid", 2000)
            .expect("the only one")
            .seqname,
        "pGM"
    );
    // two sequences, and none is this one: they are named
    let two = format!("{PLASMID_GFF}other\t.\ttRNA\t1\t70\t.\t+\t.\tID=x\n");
    let said = features_for(&parse(&two).expect("parsed"), "my plasmid", 2000).expect_err("none");
    assert!(
        said.contains("pGM") && said.contains("other") && said.contains("my plasmid"),
        "{said}"
    );
}

#[test]
fn the_ending_of_the_name_makes_a_file_an_annotation() {
    use crate::ui::platform::route::{kind_of, Kind};
    for name in ["x.gff3", "x.GFF", "x.gtf", "y.gff2"] {
        assert_eq!(
            kind_of(std::path::Path::new(name)),
            Kind::Annotation,
            "{name}"
        );
    }
}

#[test]
fn on_the_minus_strand_the_phase_that_counts_is_that_of_the_segment_at_the_right() {
    // the CDS is TATGAAA + CCCTAA: the first base belongs to a codon before it (phase 1), and the second segment goes on in frame
    let plus = "GGGGGGGGGGTATGAAAGTAAGTAAGCCCTAAGGG".to_string();
    let length = plus.len();
    // the first segment is TATGAAA, the second CCCTAA, with introns of other bases between them and flanks around
    let (a, b) = (
        plus.find("TATGAAA").expect("found") + 1,
        plus.find("CCCTAA").expect("found") + 1,
    );
    let (a_end, b_end) = (a + 6, b + 5);
    let mirror = |s: usize, e: usize| (length + 1 - e, length + 1 - s);
    let (ma, mb) = (mirror(a, a_end), mirror(b, b_end));
    let chromosome = String::from_utf8(crate::seq::dna::revcom(plus.as_bytes())).expect("text");
    let gff = format!(
        "##gff-version 3\nc\t.\tmRNA\t{}\t{}\t.\t-\t.\tID=t\nc\t.\texon\t{}\t{}\t.\t-\t.\tParent=t\nc\t.\texon\t{}\t{}\t.\t-\t.\tParent=t\n\
         c\t.\tCDS\t{}\t{}\t.\t-\t0\tParent=t\nc\t.\tCDS\t{}\t{}\t.\t-\t1\tParent=t\n",
        mb.0, ma.1, mb.0, mb.1, ma.0, ma.1, mb.0, mb.1, ma.0, ma.1
    );
    let t = &parse(&gff).expect("parsed").transcripts[0];
    let genome = genome_of(&[("c", chromosome)]);
    // the segment at the right is the first of the transcript, and its phase 1 skips the T
    assert_eq!(cds_sequence(t, &genome), Ok(Some(b"ATGAAACCCTAA".to_vec())));
    assert_eq!(protein(t, &genome, b'.'), Ok(Some(b"MKP".to_vec())));
}

#[test]
fn the_kind_of_a_gene_that_has_cds_under_it_is_mrna_and_without_it_a_transcript() {
    let a = parse(PLASMID_GFF).expect("parsed");
    assert_eq!(a.transcripts[0].kind, "mRNA");
    assert_eq!(a.transcripts[0].gene_id.as_deref(), Some("g1"));
    let bare =
        "##gff-version 3\nc\t.\tgene\t1\t50\t.\t+\t.\tID=g\nc\t.\texon\t1\t50\t.\t+\t.\tParent=g\n";
    assert_eq!(
        parse(bare).expect("parsed").transcripts[0].kind,
        "transcript"
    );
}

#[test]
fn the_three_splice_pairs_are_allowed_and_no_other_and_a_tiny_intron_is_not_one() {
    let canonical = Filters {
        canonical_splice: true,
        ..Filters::default()
    };
    let a = parse(GFF3_PLUS).expect("parsed");
    let with = |intron: &str| {
        genome_of(&[(
            "chr1",
            format!("{FLANK_LEFT}{EXON_1}{intron}{EXON_2}{FLANK_RIGHT}"),
        )])
    };
    for ok in ["GTCCCCCCAG", "GCCCCCCCAG", "ATCCCCCCAC"] {
        assert_eq!(kept(&a, &canonical, Some(&with(ok))), ["rna1"], "{ok}");
    }
    for bad in ["ATCCCCCCAG", "GTCCCCCCAC", "CTCCCCCCAC"] {
        assert!(kept(&a, &canonical, Some(&with(bad))).is_empty(), "{bad}");
    }
    // an intron of one base: no pair of sites fits in it, and nothing is cut out of range
    let tiny = "##gff-version 3\nchr1\t.\tmRNA\t11\t41\t.\t+\t.\tID=r\nchr1\t.\texon\t11\t25\t.\t+\t.\tParent=r\nchr1\t.\texon\t27\t41\t.\t+\t.\tParent=r\n";
    let t = parse(tiny).expect("parsed");
    assert!(kept(&t, &canonical, Some(&with("GTCCCCCCAG"))).is_empty());
}

#[test]
fn a_complete_cds_is_whole_codons_and_one_that_ends_in_a_partial_codon_is_not() {
    let complete = Filters {
        complete_cds: true,
        ..Filters::default()
    };
    // ATG AAA TAG + CC: the codons look complete, but the CDS is not a whole number of them
    let genome = genome_of(&[("c", "ATGAAATAGCC".to_string())]);
    let whole =
        "##gff-version 3\nc\t.\tmRNA\t1\t9\t.\t+\t.\tID=t\nc\t.\tCDS\t1\t9\t.\t+\t0\tParent=t\n";
    assert_eq!(
        kept(&parse(whole).expect("parsed"), &complete, Some(&genome)),
        ["t"]
    );
    let partial = whole
        .replace("\t9\t.\t+\t.\tID", "\t11\t.\t+\t.\tID")
        .replace("CDS\t1\t9", "CDS\t1\t11");
    assert!(kept(&parse(&partial).expect("parsed"), &complete, Some(&genome)).is_empty());
}

#[test]
fn sequences_are_written_in_lines_of_seventy() {
    let f = Folder::new("wrap");
    let long = "ACGT".repeat(40);
    std::fs::write(f.0.join("genome.fa"), format!(">c\n{long}\n")).expect("a file");
    std::fs::write(
        f.0.join("a.gff3"),
        "##gff-version 3\nc\t.\tmRNA\t1\t160\t.\t+\t.\tID=t\nc\t.\texon\t1\t160\t.\t+\t.\tParent=t\n",
    )
    .expect("a file");
    assert_eq!(
        f.run(&[&f.p("a.gff3"), "-g", &f.p("genome.fa"), "-w", &f.p("w.fa")]),
        0
    );
    let lengths: Vec<usize> = f.text("w.fa").lines().skip(1).map(str::len).collect();
    assert_eq!(lengths, [70, 70, 20]);
}

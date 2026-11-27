use super::*;

fn exon() -> Vec<String> {
    vec!["exon".to_string()]
}

const GTF: &str = "# a comment\n\
chr1\tsrc\tgene\t1\t900\t.\t+\t.\tgene_id \"g1\";\n\
chr1\tsrc\texon\t101\t200\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n\
chr1\tsrc\texon\t301\t400\t.\t+\t.\tgene_id \"g1\"; transcript_id \"t1\";\n\
chr2\tsrc\texon\t5\t9\t.\t-\t.\tgene_id \"g2\"; transcript_id \"t2\";\n\
chr2\tsrc\tCDS\t5\t9\t.\t-\t0\tgene_id \"g2\"; transcript_id \"t2\";\n";

#[test]
fn the_exons_of_a_gtf_are_grouped_by_gene_in_the_order_each_gene_is_first_met() {
    let a = parse(GTF, &exon(), "gene_id").unwrap();
    assert_eq!(a.names, ["g1", "g2"]);
    assert_eq!(a.features.len(), 3);
    let f = &a.features[0];
    assert_eq!(
        (f.seqname.as_str(), f.start, f.end, f.meta),
        ("chr1", 100, 200, 0)
    );
    assert_eq!(f.strand, Strand::Plus);
    let g = &a.features[2];
    assert_eq!(
        (g.seqname.as_str(), g.start, g.end, g.meta, g.strand),
        ("chr2", 4, 9, 1, Strand::Minus)
    );
}

#[test]
fn the_types_and_the_attribute_are_those_that_are_asked_for() {
    let by_transcript = parse(GTF, &exon(), "transcript_id").unwrap();
    assert_eq!(by_transcript.names, ["t1", "t2"]);
    let both = parse(GTF, &["CDS".to_string(), "gene".to_string()], "gene_id").unwrap();
    assert_eq!(both.features.len(), 2);
    assert_eq!(both.names, ["g1", "g2"]);
}

#[test]
fn a_line_without_the_attribute_is_named_in_the_error_and_a_type_that_is_not_there_lists_what_is() {
    let e = parse(GTF, &exon(), "locus").unwrap_err();
    assert!(
        e.contains("exon at chr1:101-200") && e.contains("no attribute locus") && e.contains("-g"),
        "{e}"
    );
    let e = parse(GTF, &["mRNA".to_string()], "gene_id").unwrap_err();
    assert!(
        e.contains("mRNA") && e.contains("CDS, exon, gene") && e.contains("-t"),
        "{e}"
    );
}

#[test]
fn gff3_attributes_are_read_and_unescaped_and_a_fasta_section_is_left_alone() {
    let gff = "##gff-version 3\nc\t.\tCDS\t1\t30\t.\t+\t0\tID=cds1;locus_tag=b%3B1\nc\t.\tCDS\t40\t60\t.\t.\t0\tID=cds2;locus_tag=b2\n##FASTA\n>c\nACGT\n";
    let a = parse(gff, &["CDS".to_string()], "locus_tag").unwrap();
    assert_eq!(a.names, ["b;1", "b2"]);
    assert_eq!(a.features[1].strand, Strand::Unknown);
}

#[test]
fn a_saf_table_is_told_by_its_header() {
    let saf = "# made by hand\nGeneID\tChr\tStart\tEnd\tStrand\tExtra\ng1\tchr1\t101\t200\t+\tx\ng2\tchr1\t301\t400\t-\ng3\tchr2\t1\t1\t.\ng1\tchr1\t501\t600\t+\r\n";
    let a = parse(saf, &exon(), "gene_id").unwrap();
    assert_eq!(a.names, ["g1", "g2", "g3"]);
    assert_eq!(a.features.len(), 4);
    assert_eq!(
        (a.features[0].start, a.features[0].end, a.features[0].strand),
        (100, 200, Strand::Plus)
    );
    assert_eq!(
        (a.features[2].start, a.features[2].end, a.features[2].strand),
        (0, 1, Strand::Unknown)
    );
    assert_eq!(a.features[1].strand, Strand::Minus);
    assert_eq!(a.features[3].meta, 0);
}

#[test]
fn a_saf_table_with_a_fault_says_which_line() {
    let table = |rows: &str| {
        parse(
            &format!("GeneID\tChr\tStart\tEnd\tStrand\n{rows}"),
            &exon(),
            "gene_id",
        )
    };
    let e = table("g\tc\t1\t5\n").unwrap_err();
    assert!(e.contains("Line 2") && e.contains("4 columns"), "{e}");
    let e = table("g\tc\tx\t5\t+\n").unwrap_err();
    assert!(
        e.contains("Line 2") && e.contains("x for a position"),
        "{e}"
    );
    let e = table("g\tc\t0\t5\t+\n").unwrap_err();
    assert!(e.contains("0 for a position"), "{e}");
    let e = table("g\tc\t9\t5\t+\n").unwrap_err();
    assert!(e.contains("starts at 9, after its end 5"), "{e}");
    assert!(table("g\tc\t5\t5\t+\n").is_ok());
    assert!(table("").unwrap_err().contains("no rows"));
}

#[test]
fn a_gtf_whose_header_line_starts_like_a_saf_is_not_one() {
    // only the exact header, GeneID and a tab, is a SAF table
    assert!(!is_saf("GeneIDs\tChr\n"));
    assert!(is_saf("\n#x\nGeneID\tChr\n"));
    assert!(!is_saf(GTF));
}

//! The table of checks. Every line here is a claim to a person who has to decide whether to trust a result: write what a test
//! really compares, the tolerance it really uses, and what nobody has compared yet. A test that is left out of the usual run
//! (`#[ignore]`) is not listed under `tests`; say in `limits` that it exists and what it needs.

use super::{Check, Standard};

const CLONING: &str = "Sequences and cloning";
const SIGNALS: &str = "Sequencing signals and reads";
const ALIGNMENT: &str = "Alignment";
const TREES: &str = "Trees";
const SPECIES: &str = "Species trees";
const STRUCTURE: &str = "Structure";
const MOTIFS: &str = "Motifs and domains";
const PIPELINES: &str = "Pipelines and other programs";

pub const CHECKS: &[Check] = &[
    Check {
        group: CLONING,
        analysis: "Melting temperature (Tm)",
        standard: Standard::Reference,
        against: "The original program's own Tm procedure, on eleven oligos, one of them at a second setting of DNA and salt.",
        result: "Equal to the tenth of a degree that is shown.",
        limits: "Eleven oligos, not a random set. The number is as good as the original's method, which takes the DNA and salt concentrations and nothing else.",
        tests: &[("src/seq/dna.rs", "melting_temperature_matches_the_original")],
    },
    Check {
        group: CLONING,
        analysis: "Primer design (Primer3)",
        standard: Standard::Reference,
        against: "primer3 itself, built from its source: its melting temperature for ten oligos, and the pairs it designs around a target and anywhere in a template of 1,000 bases.",
        result: "Melting temperatures agree to 1e-8 °C; the pairs have the same places and product sizes, and penalties within 1e-5.",
        limits: "One template. Not ported from primer3: the thermodynamic dimer and hairpin checks (its older complementarity scores are used), template mispriming, internal oligos, excluded regions and several product-size ranges; with those options primer3 chooses differently.",
        tests: &[
            ("src/seq/primer3.rs", "melting_temperature_matches_primer3"),
            ("src/seq/primer3.rs", "pairs_around_a_target_match_primer3"),
            (
                "src/seq/primer3.rs",
                "the_best_pairs_anywhere_in_a_template_match_primer3",
            ),
        ],
    },
    Check {
        group: CLONING,
        analysis: "TaqMan probe design (between the primers)",
        standard: Standard::Independent,
        against: "A plain restatement of the rules, written apart from the program for the rules themselves (the melting temperature, the dimer scores, the enumeration and the weights are the program's own, so a misreading shared there would pass): every placement, strand and length between the primers of eight pairs, judged one by one, for four settings, and, for the best assays as a whole (the best pairs, each with its best probe, for upper and lower case), against the same search over forty pairs. The probes the program lists are exactly those that pass, with the penalties that the text of the module gives. The rules are Applied Biosystems' published ones, as the defaults are written (and a test pins the defaults): no G at the 5' end, no run of four G, G+C 30 to 80 %, at most 30 bases, a melting temperature about 10 °C above the primers, primers close to the probe without overlapping it. The guidance says only to avoid runs of other bases: the limit of four, the window of 6 to 14 °C around the 10, the smallest size of 18 and a gap of one base are this program's choices.",
        result: "The same set of probes and the same penalties (within 1e-9), for every pair and setting tried (and at least one probe passes in each setting); the G+C limits are tried at exactly 30 and 80 %.",
        limits: "Not compared with Primer Express, nor with primer3's internal oligo, on the same template. The melting temperature is primer3's (SantaLucia, 50 mM salt) and is placed against the primers' by a difference; Primer Express's own numbers differ, and the 10 °C is a guideline for its calculator. The weights of the penalty are this program's choice. Minor-groove-binder probes are not modelled; no check that the probe is unique in a genome or avoids a known variant. The preference for more C than G is widely quoted but was not in the Applied Biosystems documents read, so it is only a small penalty.",
        tests: &[
            (
                "src/seq/probe/tests.rs",
                "the_probes_for_a_pair_are_all_of_those_that_pass_and_the_best_comes_first",
            ),
            (
                "src/seq/probe/tests.rs",
                "every_probe_found_follows_the_rules_and_lies_between_its_primers",
            ),
            (
                "src/seq/probe/tests.rs",
                "each_rule_on_its_own_refuses_what_it_names_and_nothing_else",
            ),
            (
                "src/seq/probe/tests.rs",
                "the_assays_are_the_best_pairs_each_with_its_best_probe_for_any_case_of_the_template",
            ),
            (
                "src/seq/probe/tests.rs",
                "the_g_plus_c_limits_are_inclusive_at_thirty_and_eighty_percent",
            ),
            (
                "src/seq/probe/tests.rs",
                "the_defaults_are_the_published_guidelines_and_the_choices_the_module_names",
            ),
        ],
    },
    Check {
        group: CLONING,
        analysis: "Annotation files (GFF3 and GTF): transcripts, sequences, filters",
        standard: Standard::Known,
        against: "Cases worked out by hand: a gene of two exons and an intron on each strand (the same gene on the minus strand is the reverse complement with its coordinates mirrored), a CDS that begins in the middle of a codon (phase 1 and 2), a stop inside, the same gene as GFF3 and as GTF, UTR and codon lines that must be left out, exons made from CDS and UTR, a CDS with no parent and a gene with CDS under it (bacterial annotations), an exon with two parents, features in several lines with one ID, names with the characters GFF3 escapes. What is written (GFF3 and GTF) is read back to the same transcripts. Every filter, every error with its line, and the command with its outputs.",
        result: "Every sequence, protein, coordinate and refusal is what the test says by hand.",
        limits: "Not compared with gffread itself: its repository's examples (an annotation, a genome, the outputs its test script expects) would be the test, and fetching them is a download not yet asked for. Taken on trust from its documentation and not verified: that the stop codon at the end is left out of the protein, that the phase trims the first segment of the CDS, and the symbol for a stop inside. Only the standard genetic code. Merging transcripts into loci (-M), BED, TLF and the other options of gffread are not done. The features put on a sequence were not compared with another program's.",
        tests: &[
            (
                "src/seq/gff/tests.rs",
                "the_spliced_transcript_the_cds_and_the_protein_are_cut_from_the_genome_as_the_annotation_says",
            ),
            (
                "src/seq/gff/tests.rs",
                "a_gene_on_the_minus_strand_gives_the_same_transcript_and_protein",
            ),
            (
                "src/seq/gff/tests.rs",
                "a_cds_that_begins_in_the_middle_of_a_codon_loses_the_bases_its_phase_says",
            ),
            (
                "src/seq/gff/tests.rs",
                "gtf_gives_the_same_transcript_and_start_and_stop_codon_lines_are_ignored_when_there_is_cds",
            ),
            (
                "src/seq/gff/tests.rs",
                "what_is_written_is_read_back_the_same_as_gff3_and_as_gtf",
            ),
            ("src/seq/gff/tests.rs", "the_filters_keep_what_they_say"),
            (
                "src/seq/gff/tests.rs",
                "the_filters_that_read_the_genome_judge_splice_sites_and_stops",
            ),
            (
                "src/seq/gff/tests.rs",
                "an_annotation_becomes_features_of_the_sequence_with_the_strand_the_gene_and_the_product",
            ),
        ],
    },
    Check {
        group: ALIGNMENT,
        analysis: "K-mer counting and spectra (after KAT)",
        standard: Standard::Independent,
        against: "A plain string-based counter written apart from the program (every window of random sequences with letters that break windows, lower case and U, canonical by the smaller of the k-mer and its reverse complement as text), for k = 1, 2, 3, 5, 11, 21, 27, 31 and 32; cases worked out by hand (the spectrum, the G+C matrix, the comparison matrix, the section of a sequence); reads simulated from a random genome of 20,000 bases (30-fold, both strands, one base in 200 wrong) whose spectrum has its peak where the k-mer coverage is and whose size comes out within 10 % of the genome; synthetic spectra with one peak and two.",
        result: "Every count equals the plain counter's; the peak of the simulated reads is within 21 to 27 for a k-mer coverage of 24, and the size within 10 % of 20,000.",
        limits: "Not compared with KAT: its program and its example data were not run (that needs a download), so no claim of the same numbers. The genome size is the textbook estimate (k-mers from the valley on over the peak), not KAT's fit of distributions, whose method the pages read do not give; where the first rise after the errors is only a small bump on a long error tail the valley may be put after it. The default k of 27 is the manual's example, not KAT's documented default. Counting is in memory (not a human genome); no compressed input; filter, blob, cvg and the plots are not done.",
        tests: &[
            (
                "src/seq/kmer/tests.rs",
                "the_counts_are_those_of_the_plain_way_for_every_length_with_other_letters_and_lower_case",
            ),
            (
                "src/seq/kmer/tests.rs",
                "the_spectrum_the_gc_matrix_and_the_comparison_are_what_they_are_worked_out_to_be",
            ),
            (
                "src/seq/kmer/tests.rs",
                "reads_of_a_simulated_genome_give_a_peak_at_the_coverage_of_the_k_mers_and_the_size_of_the_genome",
            ),
            (
                "src/seq/kmer/tests.rs",
                "a_spectrum_of_errors_and_one_peak_gives_the_valley_the_peak_and_the_size",
            ),
        ],
    },
    Check {
        group: ALIGNMENT,
        analysis: "Aligning DNA by its translation (codons kept whole)",
        standard: Standard::Known,
        against: "Sequences written by hand from a protein with one codon for each amino acid (a second codon for alanine in two of them), with two residues missing in one, one changed and two extra in another; translations worked out by hand for three reading frames, for the standard and the vertebrate mitochondrial code (TGA is a stop in one and tryptophan in the other), and with an ambiguous base, lower case, U and gaps; 300 random cases of random sequences with random gaps put in their proteins by a stand-in aligner.",
        result: "The proteins are those worked out; the alignment of the DNA has no base lost or moved (without its gaps every row is the sequence given), every gap is three, the two missing residues are six gaps, heads line up and tails are padded; an aligner that changes, shortens or lengthens a sequence is refused with the sequence named; a protein, a sequence with no whole codon and a frame or code that does not exist are refused with the reason.",
        limits: "Not compared with pal2nal, TranslatorX or MACSE: none was run. One reading frame for all the sequences, no frame shifts or introns, the frame is not found, and a stop codon is aligned as X. The progressive method was run on the proteins of one small family; MAFFT and the workstation path were not run on proteins in the tests (the workstation sends the proteins and the program puts the codons back).",
        tests: &[
            (
                "src/seq/codonalign/tests.rs",
                "the_alignment_of_the_proteins_is_carried_back_as_whole_codons",
            ),
            (
                "src/seq/codonalign/tests.rs",
                "for_random_sequences_no_base_is_lost_or_moved_and_every_gap_is_a_whole_codon",
            ),
            (
                "src/seq/codonalign/tests.rs",
                "what_the_aligner_changes_or_leaves_out_is_refused_with_the_sequence_named",
            ),
            (
                "src/seq/codonalign/tests.rs",
                "the_reading_frame_leaves_a_head_and_a_tail_and_the_genetic_code_changes_the_protein",
            ),
            (
                "src/ui/platform/alignment/tests.rs",
                "aligned_by_translation_the_dna_keeps_its_codons_whole_and_the_record_says_how_they_were_read",
            ),
            (
                "src/ui/platform/alignment/tests.rs",
                "the_translation_is_kept_with_the_page_and_run_again_from_the_record",
            ),
        ],
    },
    Check {
        group: ALIGNMENT,
        analysis: "Primers and probes picked from an alignment",
        standard: Standard::Known,
        against: "Alignments made by hand from a template of 1,000 bases that primer3's own tests use: rows that agree everywhere (the pairs must be the ones primer3's rules give for the sequence, with the same penalties), one row with one changed base (no oligo may lie on that column, and the longest run of usable columns is counted), a lower share of agreeing rows (the oligo that lies on the variable column is counted as fitting 4 of 5 rows exactly and all 5 with a mismatch), a mismatch placed at each of the last five bases of a primer's 3' end, for the left and the right primer, and in a probe (not tolerated in a primer, tolerated in a probe), gaps counted as mismatches, a majority of gaps that leaves a column out of the consensus, a gap inside the amplicon (the product is 3 bases shorter in that row), positions counted in the first row's own bases, RNA and lower case, a region, and rows that share nothing.",
        result: "The pairs, penalties, counts, sizes, positions and messages are those worked out by hand; every oligo is the consensus or its reverse complement as it is ordered; the table has one line to an oligo.",
        limits: "Not compared with another program that designs primers from an alignment: none was run. The melting temperature and the penalties are primer3's for the consensus oligo, not for the oligo with the mismatches of a row; the number of rows that an oligo fits is a count of mismatches, not a prediction of whether it primes. No degenerate (IUPAC) oligos, and nothing is checked against a genome for other binding sites. The windows are looked at as text and shapes by the tests, not by a person.",
        tests: &[
            (
                "src/seq/alnprimers/tests.rs",
                "rows_that_agree_everywhere_give_the_pairs_primer3_gives_for_the_sequence",
            ),
            (
                "src/seq/alnprimers/tests.rs",
                "a_column_where_the_rows_differ_is_in_no_oligo_when_all_must_agree",
            ),
            (
                "src/seq/alnprimers/tests.rs",
                "a_mismatch_in_the_last_five_bases_of_a_primer_is_not_tolerated_but_one_in_a_probe_is",
            ),
            (
                "src/seq/alnprimers/tests.rs",
                "the_product_is_counted_in_the_bases_of_each_row",
            ),
            (
                "src/seq/alnprimers/tests.rs",
                "a_probe_lies_between_the_primers_and_on_columns_all_the_rows_share",
            ),
            (
                "src/ui/windows/aln_primers_tests.rs",
                "the_consensus_opens_with_the_primers_as_features_in_the_right_places",
            ),
        ],
    },
    Check {
        group: TREES,
        analysis: "Tree files with annotations (NEXUS, BEAST and FigTree notes)",
        standard: Standard::Known,
        against: "Files written by hand in the way BEAST's TreeAnnotator, MrBayes and FigTree write them: a translate table with quoted names, comments before statements, a semicolon inside a quoted name, case-insensitive block names, `[&posterior=0.98,height_95%_HPD={0.5,1.5},note=\"a, b\"]` after a node, after a tip and after the closing bracket, `[&R]` and NHX comments that are no annotation.",
        result: "Every name, branch length, annotation and range is as written; a tree written with its annotations reads back equal; a file with no trees block or no tree statement is refused with a reason.",
        limits: "Not compared with FigTree or TreeAnnotator on real files: neither was run. Annotations of a taxon block (colours of taxa) and FigTree's own `figtree` block are not read; NHX is read as a comment. Annotations are kept for the session, not in the session's saved results.",
        tests: &[
            (
                "src/seq/newick.rs",
                "annotations_in_square_brackets_belong_to_the_node_before_them",
            ),
            (
                "src/seq/newick.rs",
                "a_tree_written_with_its_annotations_reads_back_the_same",
            ),
            (
                "src/seq/nexus/tests.rs",
                "trees_are_read_with_their_names_translated_tip_names_and_annotations",
            ),
            (
                "src/seq/nexus/tests.rs",
                "case_comments_quotes_and_a_semicolon_in_a_name_do_not_confuse_the_reader",
            ),
        ],
    },
    Check {
        group: TREES,
        analysis: "Editing a tree (root here, rotate, order, collapse weak branches) and support from other trees",
        standard: Standard::Known,
        against: "Trees of six tips whose branches carry lengths and a posterior: after rooting on the branch above every kind of node (tip, inner node, child of a root that has two or three children) the set of branches seen from outside, their lengths added up where the old root joined two, and their values must equal the original's; replicates built by hand for which the share of trees with each branch and with its strongest rival is counted on paper; the value on the same branch of a tree that is rooted elsewhere.",
        result: "The branches, their total lengths and values are equal before and after an edit; the shares are 0.75, 0.5 and 0.75 with rivals at 0.25; the root's place does not change what is counted.",
        limits: "Not compared with TreeGraph 2 or FigTree on the same trees: neither was run, and their documentation does not state every rule (how a rooting inside a branch shares its length, how a missing branch is valued), so these are this program's choices, stated in the code. Branches are matched on the tips that every tree has (at least four). A conflict is searched among the 5,000 branches that most of the other trees have. Rooting at the midpoint keeps node labels but not other annotations.",
        tests: &[
            (
                "src/seq/treeedit/tests.rs",
                "rooting_on_a_branch_keeps_every_branch_with_its_length_and_value",
            ),
            (
                "src/seq/treeedit/tests.rs",
                "branches_below_a_value_are_collapsed_into_their_parent",
            ),
            (
                "src/seq/treeedit/tests.rs",
                "how_often_a_branch_is_among_the_replicates_and_how_often_a_rival_is",
            ),
            (
                "src/seq/treeedit/tests.rs",
                "the_value_on_the_same_branch_of_another_tree_is_found_wherever_it_is_rooted",
            ),
        ],
    },
    Check {
        group: TREES,
        analysis: "Drawing trees: radial layout, collapsed clades, annotation colours, axis",
        standard: Standard::Known,
        against: "Layouts worked out by hand on small trees (a branch is as long on the page as its length in the radial layout; four equal tips sit at 45, 135, 225 and 315 degrees; a clade of three tips takes three times the wedge of one), the rows and the numbers of the nodes after a collapsed clade, the gradient at its ends and middle, the tick numbers of an axis from the root, from the tips and in dates, and the page drawn without a window and read as its text.",
        result: "The numbers are those worked out; every layout draws every tip name; a collapsed clade is one triangle that says its tips; a picked clade still lists the names of all its tips; the legend and the dates are on the page.",
        limits: "Looked at by the tests as text and shapes, not by a person, on a screen. The radial layout is the equal-angle method as the literature describes it, not compared with FigTree's drawing. No fish-eye zoom, no polar root angle, no pie charts or text brackets of TreeGraph 2, and the node bars assume the range is a height before the tips.",
        tests: &[
            (
                "src/seq/treelayout/folded_tests.rs",
                "the_radial_layout_puts_the_root_in_the_middle_and_every_branch_its_own_length_long",
            ),
            (
                "src/seq/treelayout/folded_tests.rs",
                "a_collapsed_clade_is_one_tip_that_stands_for_all_of_its_tips",
            ),
            (
                "src/seq/treestyle/tests.rs",
                "a_number_gets_a_colour_from_the_gradient_by_its_place_between_the_least_and_the_greatest",
            ),
            (
                "src/seq/treestyle/tests.rs",
                "ticks_are_at_round_numbers_in_the_range_of_the_axis_and_written_with_the_decimals_they_need",
            ),
            (
                "src/ui/platform/trees/tests.rs",
                "every_layout_draws_every_tip_name",
            ),
            (
                "src/ui/platform/trees/tests.rs",
                "the_names_of_a_collapsed_clade_are_still_those_of_its_tips",
            ),
        ],
    },
    Check {
        group: CLONING,
        analysis: "PCR (primer sites and products)",
        standard: Standard::Reference,
        against: "The original's own PCR code, cut out of its source and run unchanged under Tcl, on 140 random templates and primers, and its copying of the template's features into a product on 240 random features.",
        result: "Every site, every product and every copied feature is the same.",
        limits: "Not ported: the primer database, the filter for primers that are unique in the template, the features of a product that runs through the origin of a circle, and leaving out a primer feature that the template already has under that name.",
        tests: &[
            (
                "src/seq/pcr/tests.rs",
                "every_scan_finds_the_sites_the_originals_own_code_finds",
            ),
            (
                "src/seq/pcr/tests.rs",
                "every_product_is_what_the_originals_own_code_makes",
            ),
            (
                "src/seq/pcr/tests.rs",
                "the_template_features_a_product_copies_are_those_the_originals_own_code_copies",
            ),
        ],
    },
    Check {
        group: CLONING,
        analysis: "Find, with mismatches",
        standard: Standard::Reference,
        against: "The original's own search on 300 random sequences and patterns: mismatches, either strand, wrapping round a circle, each action of the dialog.",
        result: "The same thing is found and selected in every case.",
        limits: "Short random sequences; a kind of pattern that the generator never made is untried.",
        tests: &[(
            "src/seq/findmm/tests.rs",
            "every_search_does_what_the_originals_own_code_does",
        )],
    },
    Check {
        group: CLONING,
        analysis: "Open reading frames (Find Next ORF)",
        standard: Standard::Reference,
        against: "The original's own ORF code on 240 random sequences: either strand, linear and circular, each start rule and minimum length.",
        result: "The same stretch is selected in every case.",
        limits: "Random sequences, not annotated genomes: it finds what the original finds, which is every frame that is open, not which of them are genes.",
        tests: &[(
            "src/seq/nextorf/tests.rs",
            "find_next_selects_what_the_originals_own_code_selects",
        )],
    },
    Check {
        group: CLONING,
        analysis: "Translate",
        standard: Standard::Reference,
        against: "The text the original's Translate window builds, from its own code run under Tcl, for 22 cases of letters, spacing, numbering, strand and width.",
        result: "The same text, character for character.",
        limits: "22 cases, all with the standard genetic code.",
        tests: &[(
            "src/seq/translation/tests.rs",
            "the_text_equals_what_the_originals_own_code_builds",
        )],
    },
    Check {
        group: CLONING,
        analysis: "Copy Special",
        standard: Standard::Reference,
        against: "What the original's own Copy Special code puts on the clipboard, for 46 cases of its six conversions.",
        result: "The same text in every case.",
        limits: "46 cases, all with the standard genetic code.",
        tests: &[(
            "src/seq/copyspecial/tests.rs",
            "every_conversion_equals_what_the_originals_own_code_copies",
        )],
    },
    Check {
        group: CLONING,
        analysis: "Restriction sites and List Enzymes",
        standard: Standard::Reference,
        against: "The original's own site finder and List Enzymes window on 120 random sequences, 48 of them dense with sites that overlap, and the enzyme file it reads.",
        result: "The same sites in the same order, the same text of the window, and the same 116 enzyme names and patterns.",
        limits: "The default set of 116 enzymes; another enzyme file is not compared. A wider comparison, on a larger file that the driver makes with other seeds, is left out of the usual run and is run by hand (HELIXLOOM_ENZLIST_GOLDEN).",
        tests: &[
            (
                "src/seq/enzlist/tests.rs",
                "lists_equal_the_originals_on_every_case",
            ),
            (
                "src/seq/enzlist/tests.rs",
                "the_sites_behind_the_lists_equal_the_originals_in_number_and_order",
            ),
            (
                "src/seq/enzlist/tests.rs",
                "lists_and_sites_equal_the_originals_on_the_dense_cases",
            ),
            (
                "src/seq/enzlist/tests.rs",
                "the_enzyme_names_and_patterns_equal_the_enzyme_file_the_original_reads",
            ),
        ],
    },
    Check {
        group: CLONING,
        analysis: "Silent sites and diagnostic sites",
        standard: Standard::Reference,
        against: "The original's own code for each, on 80 and 30 random sequences: the stretch selected, the enzymes chosen, and the text of its window.",
        result: "The same sites and the same report in every case.",
        limits: "Random sequences; a choice of enzymes that the generator never made is untried.",
        tests: &[
            (
                "src/seq/silent/tests.rs",
                "reports_equal_the_originals_on_every_case",
            ),
            (
                "src/seq/diagnostic/tests.rs",
                "reports_equal_the_originals_on_every_case",
            ),
        ],
    },
    Check {
        group: CLONING,
        analysis: "Text Map",
        standard: Standard::Reference,
        against: "The original's own Text Map code, run against real Tk text widgets, on 70 random sequences with features, enzymes and each setting of the dialog: the text of its window and the place and colour of every mark in it.",
        result: "The same page in every case.",
        limits: "70 cases; a selection through the origin of a circle stops with the original's error, as the original does. A wider comparison, on any file the driver makes, is left out of the usual run and is run by hand (HELIXLOOM_TEXTMAP_GOLDEN).",
        tests: &[
            (
                "src/seq/textmap/tests.rs",
                "the_page_is_the_one_the_originals_code_makes_in_every_case",
            ),
            (
                "src/seq/textmap/tests.rs",
                "a_selection_through_the_origin_is_the_error_the_original_stops_with",
            ),
        ],
    },
    Check {
        group: CLONING,
        analysis: "Molecular reaction calculator and the calculator",
        standard: Standard::Reference,
        against: "The original's own procedures, run under Tcl on about 3,800 lines of random inputs: molar ratios, master mixes, unit conversions, stock dilutions and lines of arithmetic.",
        result: "The same numbers, written as the original writes them.",
        limits: "Where the original is wrong (it divides moles per litre by the volume where it should multiply), HelixLoom gives the right number and so differs from the original on purpose. Some lines are not compared: master mixes whose numbers are all typed as whole microlitres (at most a third of the 150), and a few of the 3,000 lines of arithmetic (fewer than 20) whose numbers are too big or come back as they were typed.",
        tests: &[
            (
                "src/seq/molcalc/tests.rs",
                "molar_ratio_reactions_equal_what_the_original_works_out",
            ),
            (
                "src/seq/molcalc/tests.rs",
                "master_mixes_equal_what_the_original_works_out",
            ),
            (
                "src/seq/molcalc/tests.rs",
                "conversions_equal_what_the_original_works_out_except_where_it_is_wrong",
            ),
            (
                "src/seq/molcalc/tests.rs",
                "stock_dilutions_equal_what_the_original_works_out",
            ),
            (
                "src/seq/molcalc/tests.rs",
                "lines_of_arithmetic_equal_what_the_original_gives",
            ),
        ],
    },
    Check {
        group: CLONING,
        analysis: "Digest (fragment lengths)",
        standard: Standard::Known,
        against: "Small cases whose fragments are plain to count: a circle whose fragments must add up to the whole, a line with its two ends, one cut in a circle.",
        result: "The lengths expected.",
        limits: "Lengths only: where a band runs on the drawn gel is not compared with a real gel or with the original's picture.",
        tests: &[
            (
                "src/seq/enzymes.rs",
                "circular_digest_gives_fragments_that_sum_to_the_whole",
            ),
            ("src/seq/enzymes.rs", "linear_digest_adds_start_and_end"),
            (
                "src/seq/enzymes.rs",
                "a_single_cut_in_a_circle_gives_one_full_length_fragment",
            ),
        ],
    },
    Check {
        group: CLONING,
        analysis: "Contig assembly of Sanger reads",
        standard: Standard::Known,
        against: "Made-up reads cut from a known sequence, on both strands, with sequencing errors put in and with reads that lie inside others.",
        result: "The known sequence comes back, and the errors are voted away.",
        limits: "Never compared with an assembler (phrap, CAP3) on real reads; the original has no assembler to compare with.",
        tests: &[
            (
                "src/seq/assemble.rs",
                "assembles_mixed_strand_reads_into_one_contig",
            ),
            ("src/seq/assemble.rs", "votes_away_sequencing_errors"),
            (
                "src/seq/assemble.rs",
                "closes_overlapping_fragments_into_a_circle",
            ),
        ],
    },
    Check {
        group: SIGNALS,
        analysis: "Sanger traces (.ab1)",
        standard: Standard::Known,
        against: "Files made for the tests and laid out as the sequencers write them: the bases, their qualities, the peaks, the four channels, the dates and the signal strengths.",
        result: "Read back exactly; a damaged file is refused in plain words and never makes the reader fail.",
        limits: "The comparisons with Biopython's reader and with the original's own reader on three real files (from 310, 3100 and 3730 sequencers) exist but are left out of the usual run: the files are not in the repository, so they are run by hand with HELIXLOOM_AB1_DIR set.",
        tests: &[
            (
                "src/seq/abif/tests.rs",
                "a_trace_is_read_with_its_bases_qualities_peaks_and_channels",
            ),
            (
                "src/seq/abif/tests.rs",
                "run_dates_times_and_signal_strengths_are_read_as_the_original_reads_them",
            ),
            (
                "src/seq/abif/tests.rs",
                "damaged_pointers_do_not_make_the_reader_fail",
            ),
        ],
    },
    Check {
        group: SIGNALS,
        analysis: "Sanger traces (.scf), read and written",
        standard: Standard::Reference,
        against: "The original's own writer and reader, run under Tcl on a made-up trace of forty bases, and two files of version 2.",
        result: "The file written is the original's byte for byte, and each file reads as the original reads it.",
        limits: "One made-up trace and two small files, none from a sequencer.",
        tests: &[
            (
                "src/seq/scf/tests.rs",
                "the_export_is_byte_for_byte_what_the_originals_writer_makes",
            ),
            (
                "src/seq/scf/tests.rs",
                "the_originals_own_export_reads_as_its_own_reader_reads_it",
            ),
            (
                "src/seq/scf/tests.rs",
                "version_two_files_with_one_and_two_byte_samples_read_as_the_original_reads_them",
            ),
        ],
    },
    Check {
        group: SIGNALS,
        analysis: "Trimming a trace by quality",
        standard: Standard::Known,
        against: "Made-up qualities whose best stretch is plain to see.",
        result: "The stretch expected, and the share of good bases counted in percent.",
        limits: "Not compared with another program's trimming of the same trace.",
        tests: &[
            (
                "src/seq/abif/tests.rs",
                "the_good_stretch_is_the_run_of_bases_most_likely_right",
            ),
            (
                "src/seq/abif/tests.rs",
                "the_share_of_good_bases_is_counted_in_percent",
            ),
        ],
    },
    Check {
        group: SIGNALS,
        analysis: "Sequencing reads (FASTQ): statistics and cleaning",
        standard: Standard::Known,
        against: "Small reads made for the tests, cut and filtered by the rules read from fastp's source: its defaults, its adapter match and its sliding window.",
        result: "As the rules say.",
        limits: "Never run side by side with fastp or FastQC on the same file. Not ported: paired reads, finding the adapter (one is given), fastp's fallback for an adapter with one gap, poly-G and poly-X trimming, the low-complexity filter, removing duplicates.",
        tests: &[
            ("src/seq/fastq.rs", "the_default_filters_are_fastps"),
            ("src/seq/fastq.rs", "the_sliding_window_cuts_follow_fastp"),
            ("src/seq/fastq.rs", "statistics_count_bases_quality_and_gc"),
        ],
    },
    Check {
        group: ALIGNMENT,
        analysis: "Alignment (built in)",
        standard: Standard::Independent,
        against: "The plain aligner that the fast one replaced, kept as its oracle, on hundreds of pairs (DNA and protein letters, changes and gaps, repeats, free ends or not), and a whole alignment of many sequences made both ways.",
        result: "The same alignment and the same score to the last bit.",
        limits: "This shows that the fast code is the plain method, not how good the method is: a centre-star alignment was never scored on a benchmark set (BAliBASE) nor against MAFFT on distant sequences. One test compares with MAFFT where it is installed, on sequences with no gaps to place.",
        tests: &[
            (
                "src/seq/align.rs",
                "the_aligner_equals_the_reference_on_every_pair_whole",
            ),
            (
                "src/seq/align.rs",
                "eight_at_a_time_equals_the_reference_on_every_pair_whole",
            ),
            (
                "src/seq/msa.rs",
                "the_alignments_eight_at_a_time_give_the_same_msa_as_one_pair_at_a_time",
            ),
        ],
    },
    Check {
        group: ALIGNMENT,
        analysis: "Alignment (progressive, after ClustalW)",
        standard: Standard::Known,
        against: "Sequences evolved on a computer along a balanced tree of 8 leaves from one ancestor of 110 residues, with substitutions and insertions and deletions, so that the true alignment is known (proteins change in the way BLOSUM62 favours). Per branch 2 % (nucleotides) or 4 % (proteins) of the residues change, so two leaves are 2, 4 or 6 branches apart, about 5 on average. The share of the truly homologous pairs of residues that the alignment finds, and the share of the pairs it finds that are true, on four seeds. Also: every alignment, made with either guide tree, with end gaps free or costly, from lower case, RNA, ambiguity codes, stops and one-letter sequences, is its sequences with gaps in rows of one length and no empty column; refining never lowers the score and does raise it for some of 24 simulated sets; the score of an alignment worked out by hand; the distance of two proteins and of RNA, and that two sequences with one letter in common at their ends are far apart; BLOSUM62 against its published entries.",
        result: "The tests demand that the plain method finds 0.88 (nucleotides, 0.02 per branch) and 0.82 (proteins, 0.04) of the true pairs, and that refining loses no more than 0.02. A printout of the ignored comparison test, on 8 other seeds (the ones the gap costs were chosen on, so it is not a fresh test), gives for nucleotides at 0.02 per branch 0.93 plain, 0.95 refined (3 passes) and 0.95 for the centre-star method; at 0.04: 0.82, 0.84 and 0.81; at 0.07: 0.51, 0.49 and 0.53. For proteins at 0.04: 0.88, 0.91 and 0.92; at 0.06: 0.77, 0.81 and 0.82; at 0.1: 0.62, 0.64 and 0.61. So the method is about as accurate as the centre-star method on these: a little behind at low divergence, a little ahead at some higher ones, and refining raises the score and mostly the accuracy but at 0.07 it lowered it.",
        limits: "Not compared with ClustalW, ClustalX, MAFFT or MUSCLE themselves (none is run here for this), and not on a benchmark of real alignments (BAliBASE): the method follows ClustalW's published description and is not a port of its code, so its alignments are not ClustalW's. The gap costs were chosen on simulated sequences (the same 8 seeds as the printout), which have more insertions and deletions than most real data. Not done from ClustalW: several matrices by divergence (BLOSUM62 only), residue-specific gap penalties, the hydrophobic-stretch rule and delaying the divergent sequences. Not tested: groups of more than 60 sequences (they are not refined), 2,000 sequences, and the automatic change to word distances beyond its rule. Ambiguity codes score -1 against themselves. Two ignored tests print the numbers.",
        tests: &[
            (
                "src/seq/clustal/tests.rs",
                "the_alignment_is_close_to_the_true_one_for_sequences_that_evolved_with_substitutions_and_indels",
            ),
            (
                "src/seq/clustal/tests.rs",
                "every_alignment_is_its_sequences_with_gaps_in_rows_of_one_length_with_no_empty_column",
            ),
            (
                "src/seq/clustal/tests.rs",
                "lower_case_rna_ambiguity_codes_stops_and_one_letter_sequences_are_aligned_without_harm",
            ),
            (
                "src/seq/clustal/tests.rs",
                "refinement_never_lowers_the_score_and_never_changes_what_a_row_is",
            ),
            (
                "src/seq/clustal/tests.rs",
                "refinement_improves_the_score_of_some_alignments_and_not_only_leaves_it_alone",
            ),
            (
                "src/seq/clustal/tests.rs",
                "the_score_of_an_alignment_is_what_the_costs_say_worked_out_by_hand",
            ),
            (
                "src/seq/clustal/tests.rs",
                "unrelated_sequences_that_share_a_letter_at_their_ends_are_far_apart_not_identical",
            ),
            (
                "src/seq/clustal/tests.rs",
                "the_distance_of_proteins_is_kimuras_correction_of_the_share_that_differs",
            ),
            ("src/seq/clustal/tests.rs", "blosum62_is_the_published_matrix"),
            (
                "src/seq/clustal/tests.rs",
                "the_fast_distances_come_from_shared_words_and_still_give_an_alignment",
            ),
        ],
    },
    Check {
        group: ALIGNMENT,
        analysis: "Clustal format (.aln) and its conservation marks",
        standard: Standard::Known,
        against: "The groups of residues that ClustalW documents for its marks (nine strong and eleven weak, as the UCSF Chimera pages give them), written again in a test apart from the program's own lists, and columns made by hand; a write then a read of an alignment of three blocks (nucleotides, with gaps) and of one block (a protein); text laid out as ClustalW writes it (with its version, counts and marks), and the headers of MAFFT, Clustal Omega, MUSCLE, Kalign and ProbCons over a plain body; and refusals: text that is not a Clustal alignment, rows of different lengths, blocks out of order, and a count that shows that a block is missing.",
        result: "Every group gives its mark, every count at the end of a line is the residues so far, and every refusal names its line, as each test says by hand.",
        limits: "Not compared with a file that ClustalW or any of the other programs wrote, nor read by them: the layouts are written from their documentation and from memory, so a program that lays the format out in another way may be refused or misread. A name that holds a space gets an underscore on writing. Nucleotides get the star only. The count at the end of a line is checked against the residues or the columns so far, so a program that counts in another way would be refused.",
        tests: &[
            (
                "src/seq/clustal/aln/tests.rs",
                "proteins_get_clustalws_strong_and_weak_groups_and_a_gap_blanks_the_column",
            ),
            (
                "src/seq/clustal/aln/tests.rs",
                "every_group_that_clustalw_documents_gives_its_mark_and_a_letter_outside_it_takes_it_away",
            ),
            (
                "src/seq/clustal/aln/tests.rs",
                "nucleotides_get_stars_where_every_row_agrees_and_blanks_otherwise",
            ),
            (
                "src/seq/clustal/aln/tests.rs",
                "a_text_is_written_in_blocks_of_sixty_with_counts_and_marks_under_each",
            ),
            (
                "src/seq/clustal/aln/tests.rs",
                "what_is_written_is_read_back_whole",
            ),
            (
                "src/seq/clustal/aln/tests.rs",
                "a_file_of_clustalw_with_its_version_marks_counts_and_blank_lines_is_read",
            ),
            (
                "src/seq/clustal/aln/tests.rs",
                "files_of_other_programs_are_read_with_no_counts_with_no_marks_and_with_their_own_headers",
            ),
            (
                "src/seq/clustal/aln/tests.rs",
                "a_text_that_is_not_a_clustal_alignment_is_refused_in_words",
            ),
            (
                "src/seq/clustal/aln/tests.rs",
                "the_counts_at_the_ends_of_lines_catch_a_block_that_is_missing_and_a_count_of_columns_is_taken_too",
            ),
        ],
    },
    Check {
        group: ALIGNMENT,
        analysis: "Conservation scores",
        standard: Standard::Independent,
        against: "A separate, plain implementation of the same published formulas (Capra and Singh 2007, with the sequence weights of Henikoff and Henikoff 1994), written for the purpose, on eight proteins of 22 columns; and cases worked by hand.",
        result: "Every column within 0.0001.",
        limits: "The accuracy that the paper reports for finding functional residues was not reproduced: no benchmark data, and the authors' program was not run. Treat a score as a guide to where to look. A second implementation checks the arithmetic, not the data: it once passed a table with three background frequencies swapped, which a test of its own now catches.",
        tests: &[
            (
                "src/seq/conservation/tests.rs",
                "the_jensen_shannon_divergence_of_the_protein_example_matches_the_reference",
            ),
            (
                "src/seq/conservation/tests.rs",
                "the_entropy_of_the_protein_example_matches_the_reference",
            ),
            (
                "src/seq/conservation/tests.rs",
                "sequence_weights_follow_henikoff_and_henikoff_by_hand",
            ),
            (
                "src/seq/conservation/tests.rs",
                "a_column_of_one_residue_scores_what_its_background_frequency_says_whatever_order_the_table_is_kept_in",
            ),
        ],
    },
    Check {
        group: TREES,
        analysis: "Distance trees (neighbour-joining, UPGMA)",
        standard: Standard::Known,
        against: "The distance formulas worked by hand, a matrix that fits one tree exactly (neighbour-joining must give that tree back), and the plain neighbour-joining that the fast one replaced.",
        result: "The formulas and the tree come out as expected; the fast joins equal the plain ones to the last bit.",
        limits: "Not compared with MEGA, PHYLIP or another tree program on a real alignment.",
        tests: &[
            ("src/seq/tree.rs", "distance_models_follow_their_formulas"),
            (
                "src/seq/tree.rs",
                "neighbour_joining_reproduces_an_additive_matrix",
            ),
            (
                "src/seq/tree.rs",
                "upgma_joins_the_closest_pairs_and_ends_every_leaf_at_the_same_height",
            ),
            (
                "src/seq/tree/joins.rs",
                "the_joins_equal_the_plain_way_to_the_last_bit_on_every_kind_of_matrix",
            ),
        ],
    },
    Check {
        group: TREES,
        analysis: "Maximum-likelihood trees",
        standard: Standard::Independent,
        against: "The likelihood of a small tree worked out separately, in Python, by a sum over every state of its inner nodes (the number is kept in the test); and trees whose answer is known: the search must find the true split, and put right a tree that is one swap off.",
        result: "The likelihood agrees with the separate sum to 1e-8; the known trees are found.",
        limits: "Neither the likelihoods nor the trees were compared with IQ-TREE, RAxML or PhyML on the same alignment. Four DNA models (JC69, K80, F81, HKY85) with no rate variation between sites, and a search that takes the first swap that helps, so it can stop at a tree that is not the best.",
        tests: &[
            (
                "src/seq/likelihood.rs",
                "the_likelihood_matches_a_brute_force_sum",
            ),
            (
                "src/seq/likelihood.rs",
                "the_full_search_returns_the_right_tree_and_a_model_by_bic",
            ),
            (
                "src/seq/likelihood.rs",
                "a_swap_across_an_inner_edge_changes_the_topology_and_the_search_undoes_a_wrong_one",
            ),
            (
                "src/seq/likelihood.rs",
                "from_a_tree_one_swap_off_the_search_with_the_scan_puts_it_right_as_the_search_without_it_does",
            ),
            (
                "src/seq/likelihood.rs",
                "the_scan_of_one_branch_gives_the_likelihood_of_the_whole_tree_for_every_length",
            ),
        ],
    },
    Check {
        group: TREES,
        analysis: "Bayesian trees",
        standard: Standard::Known,
        against: "Two cases whose answer is known: with no data the sampler must give back its priors, and with clear data it must find the true split in every run.",
        result: "Both hold.",
        limits: "Not compared with MrBayes or BEAST on the same alignment; whether the chains have run long enough is judged only by the runs agreeing.",
        tests: &[
            (
                "src/seq/likelihood/bayes.rs",
                "with_no_data_the_sampler_returns_its_priors",
            ),
            (
                "src/seq/likelihood/bayes.rs",
                "the_sampler_finds_the_true_split_and_the_runs_agree",
            ),
        ],
    },
    Check {
        group: TREES,
        analysis: "Bootstrap support",
        standard: Standard::Known,
        against: "An alignment of four sequences with one clear split: the split must be supported in most of 40 resamples, and the same seed must give the same numbers.",
        result: "The split gets more than 60 %, and the numbers repeat.",
        limits: "A weak check of the numbers themselves: it asks for a support above 60 %, not for a value worked out another way. Not compared with another program's bootstrap of the same alignment.",
        tests: &[(
            "src/seq/tree.rs",
            "bootstrap_supports_a_clear_split_and_repeats_with_the_same_seed",
        )],
    },
    Check {
        group: TREES,
        analysis: "Trees on the GPUs (distances, likelihood)",
        standard: Standard::Independent,
        against: "The processor's own result for the same alignment and tree, on pools of one, two and three devices.",
        result: "The distances are the processor's to the last digit; the likelihood agrees to two parts in a million.",
        limits: "Run on the processor standing in for a card and on a Mac's GPU, never on a CUDA card. These tests are in the build with GPU support, which is the usual one. The same comparisons through PyTorch (the tests under src/accel/torch and src/seq/likelihood/torch) need a Python with PyTorch (HELIXLOOM_TORCH_PYTHON): without one they pass without checking, so they are not listed here.",
        tests: &[
            (
                "src/accel/distance/tests.rs",
                "the_gpu_distances_are_those_of_the_processor_to_the_last_digit",
            ),
            (
                "src/seq/likelihood/gpu/tests.rs",
                "a_tree_has_the_same_likelihood_on_the_cards_as_on_the_processor_for_every_model",
            ),
        ],
    },
    Check {
        group: TREES,
        analysis: "PhyloCNN (epidemic parameters from a tree)",
        standard: Standard::Reference,
        against: "The authors' own encoder, their Keras layers with their trained weights, and their pipeline from tree to estimate, on the Zurich HIV tree.",
        result: "The encoding is theirs; every network's output agrees to one part in ten thousand, and the estimates within 0.3 %.",
        limits: "One tree. The trained weights are a local add-on that is not in the repository: on a computer without them these tests pass without checking. How well the networks estimate is the authors' claim and was not measured here.",
        tests: &[
            (
                "src/seq/phylocnn/tests.rs",
                "the_encoding_matches_the_authors_encoder",
            ),
            (
                "src/seq/phylocnn/tests.rs",
                "every_network_matches_keras_on_the_same_input",
            ),
            (
                "src/seq/phylocnn/tests.rs",
                "the_estimate_agrees_with_the_authors_pipeline_within_three_parts_in_a_thousand",
            ),
        ],
    },
    Check {
        group: TREES,
        analysis: "Tree tools (rooting, ladderizing, distances, Robinson-Foulds)",
        standard: Standard::Reference,
        against: "ETE 3 (ladderizing, the Robinson-Foulds distance, distances between tips) and DendroPy (midpoint rooting), on five random trees of ten tips.",
        result: "The same trees and the same distances.",
        limits: "Five small trees.",
        tests: &[
            ("src/seq/treetools.rs", "ladderizing_matches_ete"),
            ("src/seq/treetools.rs", "midpoint_rooting_matches_dendropy"),
            ("src/seq/treetools.rs", "patristic_distances_match_ete"),
            ("src/seq/treetools.rs", "robinson_foulds_matches_ete"),
        ],
    },
    Check {
        group: SPECIES,
        analysis: "Species tree from gene trees",
        standard: Standard::Reference,
        against: "ASTRAL 5.7.8, run on twelve simulated data sets of 6 to 80 species with polytomies and missing species: its scores of given trees, what it writes on each branch, and the best scores its searches found.",
        result: "Quartet scores are equal; branch lengths and quartet supports agree to 1e-9 and posterior probabilities to 1e-6, which is as close as ASTRAL's own are to the formula; the search reaches ASTRAL's score.",
        limits: "Simulated data of up to 80 species: no real data set and none of hundreds of species was compared. ASTRAL's exact search was compared on the six data sets that are small enough for it.",
        tests: &[
            (
                "src/seq/astral/golden_tests.rs",
                "scores_and_annotations_equal_astrals_on_every_species_tree_it_scored",
            ),
            (
                "src/seq/astral/golden_tests.rs",
                "the_exact_search_finds_the_optimum_that_astrals_exact_search_found",
            ),
            (
                "src/seq/astral/golden_tests.rs",
                "the_search_reaches_astrals_scores_with_and_without_the_extra_clusters",
            ),
        ],
    },
    Check {
        group: SPECIES,
        analysis: "Weighted species trees (by support and branch length)",
        standard: Standard::Independent,
        against: "The definition: the weight of every quartet, found by walking small random trees one quartet at a time; and with every weight at one, the plain score and tree.",
        result: "The fast counts are the sums of the quartets' weights; weights of one give the plain result.",
        limits: "Not compared with wastral (ASTER), and its weighted branch support is not built.",
        tests: &[
            (
                "src/seq/astral/wtests.rs",
                "the_fast_weighted_counts_are_the_weights_of_the_quartets_one_by_one",
            ),
            (
                "src/seq/astral/wtests.rs",
                "weights_of_one_give_the_tree_and_the_score_of_astral",
            ),
        ],
    },
    Check {
        group: SPECIES,
        analysis: "Gene families (duplications and losses)",
        standard: Standard::Independent,
        against: "The definition: the quartet classes of small random family trees, counted one at a time; and without duplications, the plain score and tree.",
        result: "The score is the number of classes a species tree agrees with; families without duplications give the plain result.",
        limits: "Not compared with ASTRAL-Pro, and there are no branch supports with families.",
        tests: &[
            (
                "src/seq/astral/ftests.rs",
                "the_score_of_a_species_tree_is_the_number_of_classes_it_agrees_with",
            ),
            (
                "src/seq/astral/ftests.rs",
                "without_duplications_the_families_give_astrals_score_and_tree",
            ),
            (
                "src/seq/astral/ftests.rs",
                "the_species_tree_is_found_from_families_that_have_duplications_and_losses",
            ),
        ],
    },
    Check {
        group: SPECIES,
        analysis: "Bootstrap of species trees, with the replicate trees and their consensus",
        standard: Standard::Known,
        against: "Genes whose answer is known (all alike must give full support; split 40 against 35, a support between none and full; the same seed the same numbers), and a greedy consensus made the plain way on random trees.",
        result: "Full support where the genes agree, a support between 20 % and 95 % where they are split, and the same consensus as the plain way.",
        limits: "The support of split genes is checked to lie in a range, not against a value worked out another way. Not compared with ASTRAL's own bootstrap over genes and sites, nor with its consensus.",
        tests: &[
            (
                "src/seq/astral/btests.rs",
                "genes_that_all_say_the_same_give_every_branch_full_support",
            ),
            (
                "src/seq/astral/btests.rs",
                "genes_that_are_split_give_a_support_between_none_and_full_that_is_the_same_for_the_same_seed",
            ),
            (
                "src/seq/astral/ctests.rs",
                "the_choice_is_the_same_as_the_plain_one_and_follows_the_rules_of_a_greedy_consensus_on_random_trees",
            ),
        ],
    },
    Check {
        group: SPECIES,
        analysis: "Polytomy test",
        standard: Standard::Known,
        against: "The chi-squared statistic and its tail worked by hand for given quartet frequencies (Sayyari and Mirarab 2018).",
        result: "The statistic and the p-value expected.",
        limits: "Not compared with ASTRAL's own test (its option -t 10) on the same trees.",
        tests: &[
            (
                "src/seq/astral/polytomy.rs",
                "the_statistic_and_its_tail_are_the_ones_of_a_chi_squared_on_two_degrees",
            ),
            (
                "tests/astral.rs",
                "the_polytomy_test_labels_each_branch_with_the_p_value_the_chi_squared_of_its_frequencies_gives",
            ),
        ],
    },
    Check {
        group: STRUCTURE,
        analysis: "Reading a structure (PDB)",
        standard: Standard::Known,
        against: "Ubiquitin (PDB 1UBQ): its residues, its chain, its sequence and the helices and strands the file lists.",
        result: "Read as the file says.",
        limits: "The first model and the first conformer only are read, and water is left out.",
        tests: &[
            (
                "src/seq/pdb/tests.rs",
                "ubiquitin_is_read_with_its_residues_chain_sequence_and_annotated_elements",
            ),
            (
                "src/seq/pdb/tests.rs",
                "only_the_first_model_and_the_first_conformer_are_read_and_water_is_left_out",
            ),
        ],
    },
    Check {
        group: STRUCTURE,
        analysis: "Secondary structure",
        standard: Standard::Reference,
        against: "PyDSSP 0.9.1 on ubiquitin (PDB 1UBQ), and the helices and strands the file itself lists.",
        result: "At most 2 of the 76 residues differ from PyDSSP.",
        limits: "One protein, and three states only: helix, strand and coil.",
        tests: &[
            (
                "src/seq/secstruct/tests.rs",
                "ubiquitin_comes_out_as_pydssp_has_it",
            ),
            (
                "src/seq/secstruct/tests.rs",
                "ubiquitin_matches_the_helices_and_strands_the_file_lists",
            ),
        ],
    },
    Check {
        group: STRUCTURE,
        analysis: "Interface between chains (buried surface, contacts, hydrogen bonds)",
        standard: Standard::Independent,
        against: "A separate implementation (its own reader, the same points and radii) on ubiquitin whole and cut in two; the exact areas of one atom and of two that meet; and a real antibody with its antigen (PDB 1VFB).",
        result: "Surface areas agree within 0.05 % and the buried area within 0.2 %; in 1VFB the interface falls on the antibody's CDRs and on the two stretches of lysozyme that are known for this complex.",
        limits: "Not compared with NACCESS, FreeSASA or PISA. A hydrogen bond is a candidate by distance, with no angle; an ion is given the radius of carbon; the area that the literature gives for 1VFB was not looked up.",
        tests: &[
            (
                "src/seq/interface/tests.rs",
                "the_total_surface_of_ubiquitin_is_what_an_independent_program_gives",
            ),
            (
                "src/seq/interface/tests.rs",
                "ubiquitin_cut_in_two_gives_the_interface_an_independent_program_finds",
            ),
            (
                "src/seq/interface/tests.rs",
                "two_atoms_that_meet_bury_the_cap_that_geometry_says",
            ),
            (
                "src/seq/interface/tests.rs",
                "an_antibody_and_its_antigen_from_the_pdb",
            ),
        ],
    },
    Check {
        group: STRUCTURE,
        analysis: "Structure prediction (a folding program)",
        standard: Standard::Untested,
        against: "Nothing here: a predicted structure is the folding program's own.",
        result: "HelixLoom shows the confidence the program reports and adds no check of its own.",
        limits: "No comparison of predicted with known structures has been run here. Treat a prediction as a hypothesis, the more so where the confidence is low.",
        tests: &[],
    },
    Check {
        group: MOTIFS,
        analysis: "Motif Finder (scan, discovery, comparison)",
        standard: Standard::Known,
        against: "Motifs planted in random sequences (a word with one base in ten changed, on either strand), and the p-value of a score against scoring every word of a short motif one by one.",
        result: "The planted motif is found where it was put, and the p-value is the share of words that score as much.",
        limits: "Never compared with MEME on real promoter sets, so it has no claim to MEME's accuracy. The MEME text it writes was read back by HelixLoom, not by FIMO or Tomtom. A big set is searched less thoroughly, to bound the time.",
        tests: &[
            (
                "src/seq/motif/tests.rs",
                "a_motif_planted_in_random_sequences_on_either_strand_is_found_where_it_was_put",
            ),
            (
                "src/seq/motif/tests.rs",
                "the_p_value_of_a_score_is_the_share_of_all_words_that_score_so_much_and_the_threshold_for_a_p_value_is_the_lowest_that_fits",
            ),
            (
                "tests/motif.rs",
                "a_planted_motif_comes_out_as_meme_text_and_a_table_of_the_places",
            ),
        ],
    },
    Check {
        group: MOTIFS,
        analysis: "Protein domains (HMMER, Pfam)",
        standard: Standard::Known,
        against: "In every run, tables made to the format in the HMMER 3.4 User's Guide. By hand, with HMMER 3.4 and Pfam 39.0 (November 2026): the protein of HMMER's own tutorial against its two tutorial families, compared with the table that hmmscan writes when it is run directly; and the four proteins of the structures that come with the program (ubiquitin, the two variable domains of an antibody, lysozyme) against the whole of Pfam.",
        result: "The ten domains of the tutorial protein with the places and scores hmmscan gives, and the seven of them that pass Pfam's own thresholds. The four proteins are found in the families they are known for (PF00240, PF07686 twice, PF00062), each over more than 85 % of its length; the search took about two seconds on the Mac it was run on.",
        limits: "The two runs with a real HMMER are left out of the usual run, because HMMER and Pfam are not in the repository: `HELIXLOOM_HMMER_TUTORIAL=<the tutorial folder of HMMER> cargo test --bin helixloom -- real_hmmer real_pfam --ignored` on a computer that has them (src/remote/domains/real.rs). Five proteins in all, none of them large. Pfam's clan rule is not built, so every family of a clan that matches is listed: in that run ubiquitin showed three that overlap, which no test asks for. The tests that run a stand-in script are compiled on Unix only, so they are not listed here.",
        tests: &[
            (
                "src/seq/hmmer/tests.rs",
                "the_domains_of_a_table_are_read_with_their_places_scores_and_descriptions",
            ),
            (
                "src/pipeline/domains/tests.rs",
                "a_cutoff_is_checked_and_a_pressed_database_has_its_four_files",
            ),
        ],
    },
    Check {
        group: PIPELINES,
        analysis: "Pipeline files (run, cache, export)",
        standard: Standard::Known,
        against: "A pipeline run as a job, twice: the second run must find the answers the first kept. The Snakemake and Nextflow files it exports are compared with the text expected.",
        result: "As expected. The Nextflow file was also run once by hand, on three loci, and gave the species tree that HelixLoom's own run gives.",
        limits: "The Snakefile has never been run or read by Snakemake, which is not on the computer it was built on.",
        tests: &[
            (
                "tests/remote.rs",
                "a_pipeline_file_runs_in_a_job_and_a_second_run_finds_the_answers_the_first_kept",
            ),
            (
                "src/pipeline/flow/export/tests.rs",
                "a_pipeline_of_other_programs_is_a_snakefile_that_calls_them_as_the_engine_does",
            ),
            (
                "src/pipeline/flow/export/tests.rs",
                "the_same_pipeline_is_a_nextflow_script",
            ),
        ],
    },
    Check {
        group: PIPELINES,
        analysis: "Jobs on the workstation",
        standard: Standard::Known,
        against: "A real server started for the test: a job runs in a process of its own, its result waits for the person, and results survive a restart.",
        result: "The result comes back and is kept until it is fetched, across a restart too. The function that computes a job there is the one that computes it here.",
        limits: "Tried on a Mac with the server on the same computer: never on Linux or Windows, and never with CUDA cards.",
        tests: &[
            (
                "tests/remote.rs",
                "a_job_runs_in_its_own_process_and_its_result_waits_for_the_person",
            ),
            (
                "tests/remote.rs",
                "results_and_waiting_jobs_survive_a_restart_of_the_server",
            ),
        ],
    },
    Check {
        group: PIPELINES,
        analysis: "MAFFT, IQ-TREE and the other programs HelixLoom runs",
        standard: Standard::Untested,
        against: "Nothing here: an alignment by MAFFT or a tree by IQ-TREE is that program's own answer.",
        result: "HelixLoom hands over the input and reads the output; `helixloom tools` says which version it found.",
        limits: "Their accuracy is theirs to state: read their papers. HelixLoom does not check again what they return.",
        tests: &[],
    },
];

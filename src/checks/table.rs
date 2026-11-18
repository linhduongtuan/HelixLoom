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

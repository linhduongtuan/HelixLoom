use super::*;
use crate::pipeline::flow::Pipeline;

fn pipeline(text: &str) -> Pipeline {
    Pipeline::parse(text).expect("a pipeline file")
}

const FULL: &str = "helixloom-pipeline 1\nname: Align, gene trees, then a species tree\ninput: sequences\nstep align engine=mafft\nstep tree engine=iqtree bootstrap=1000\nstep species\n";

fn source() -> Source {
    Source {
        folder: "/data/loci".into(),
        ext: "fasta".into(),
    }
}

const SNAKE_FULL: &str = "\
# Written by HelixLoom from the pipeline \"Align, gene trees, then a species tree\". Run it with:  snakemake --cores 4
# It needs mafft, iqtree3 and helixloom on the PATH. Sequence names must be plain (letters, digits, - _ .) for IQ-TREE.

SAMPLES, = glob_wildcards(\"/data/loci/{sample}.fasta\")

wildcard_constraints:
    sample=\"[^/]+\"

rule all:
    input: \"results/species.tre\"

rule align:
    input: \"/data/loci/{sample}.fasta\"
    output: \"results/aligned/{sample}.fasta\"
    shell: \"mafft --auto --quiet {input} > {output}\"

rule tree:
    input: \"results/aligned/{sample}.fasta\"
    output: \"results/trees/{sample}.treefile\"
    shell: \"iqtree3 -s {input} -m MFP -B 1000 -seed 1 -redo -quiet --prefix results/trees/{wildcards.sample}\"

rule gather:
    input: expand(\"results/trees/{sample}.treefile\", sample=SAMPLES)
    output: \"results/gene_trees.tre\"
    shell: \"cat {input} > {output}\"

rule species:
    input: \"results/gene_trees.tre\"
    output: \"results/species.tre\"
    shell: \"helixloom astral -i {input} -m 4 -o {output}\"
";

#[test]
fn a_pipeline_of_other_programs_is_a_snakefile_that_calls_them_as_the_engine_does() {
    let text = snakefile(&pipeline(FULL), &source(), &Programs::default()).expect("a Snakefile");
    assert_eq!(text, SNAKE_FULL);
}

#[test]
fn the_same_pipeline_is_a_nextflow_script() {
    let text = nextflow(&pipeline(FULL), &source(), &Programs::default()).expect("a script");
    let expected = "\
// Written by HelixLoom from the pipeline \"Align, gene trees, then a species tree\". Run it with:  nextflow run main.nf
// It needs mafft, iqtree3 and helixloom on the PATH. Sequence names must be plain (letters, digits, - _ .) for IQ-TREE.
nextflow.enable.dsl = 2

params.input = \"/data/loci/*.fasta\"
params.outdir = \"${projectDir}/results\"

process ALIGN {
    publishDir \"${params.outdir}/aligned\", mode: 'copy'
    input:
    path locus
    output:
    path \"${locus.baseName}.aligned.fasta\"
    script:
    \"\"\"
    mafft --auto --quiet \"${locus}\" > \"${locus.baseName}.aligned.fasta\"
    \"\"\"
}

process TREE {
    publishDir \"${params.outdir}/trees\", mode: 'copy'
    input:
    path alignment
    output:
    path \"${alignment.baseName}.treefile\"
    script:
    \"\"\"
    iqtree3 -s \"${alignment}\" -m MFP -B 1000 -seed 1 -redo -quiet --prefix \"${alignment.baseName}\"
    \"\"\"
}

process SPECIES {
    publishDir params.outdir, mode: 'copy'
    input:
    path trees
    output:
    path \"species.tre\"
    script:
    \"\"\"
    cat ${trees} > gene_trees.tre
    helixloom astral -i gene_trees.tre -m 4 -o species.tre
    \"\"\"
}

workflow {
    loci = Channel.fromPath(params.input)
    aligned = ALIGN(loci)
    gene_trees = TREE(aligned)
    SPECIES(gene_trees.collect())
}
";
    assert_eq!(text, expected);
}

#[test]
fn a_pipeline_that_starts_from_alignments_has_no_align_rule_and_one_that_ends_in_trees_has_no_species(
) {
    let aligned = pipeline(
        "helixloom-pipeline 1\nname: Gene trees\ninput: alignment\nstep tree engine=iqtree\n",
    );
    let snake = snakefile(&aligned, &source(), &Programs::default()).expect("a Snakefile");
    assert!(!snake.contains("rule align"), "{snake}");
    assert!(
        snake.contains(
            "input: \"/data/loci/{sample}.fasta\"\n    output: \"results/trees/{sample}.treefile\""
        ),
        "{snake}"
    );
    assert!(
        snake.contains(
            "rule all:\n    input: expand(\"results/trees/{sample}.treefile\", sample=SAMPLES)"
        ),
        "{snake}"
    );
    assert!(!snake.contains("rule species") && !snake.contains("gather"));
    assert!(!snake.contains("-B"), "no bootstrap was asked for");
    let flow = nextflow(&aligned, &source(), &Programs::default()).expect("a script");
    assert!(
        !flow.contains("process ALIGN") && !flow.contains("process SPECIES"),
        "{flow}"
    );
    assert!(flow.contains("TREE(loci)"), "{flow}");
}

#[test]
fn a_folder_of_tree_files_gets_a_species_tree_from_each_and_no_gene_tree_is_left_out() {
    let trees = pipeline("helixloom-pipeline 1\nname: Species trees\ninput: trees\nstep species\n");
    let src = Source {
        folder: "/data/trees".into(),
        ext: "tre".into(),
    };
    let snake = snakefile(&trees, &src, &Programs::default()).expect("a Snakefile");
    assert!(snake.contains("rule species:\n    input: \"/data/trees/{sample}.tre\"\n    output: \"results/species/{sample}.tre\"\n    shell: \"helixloom astral -i {input} -o {output}\""), "{snake}");
    assert!(
        !snake.contains("-m 4"),
        "the leaf rule is for gathered gene trees only: {snake}"
    );
    assert!(snake.contains(
        "rule all:\n    input: expand(\"results/species/{sample}.tre\", sample=SAMPLES)"
    ));
    let flow = nextflow(&trees, &src, &Programs::default()).expect("a script");
    assert!(
        flow.contains("helixloom astral -i \"${trees}\" -o \"${trees.baseName}.species.tre\""),
        "{flow}"
    );
}

#[test]
fn the_names_of_the_programs_are_the_ones_given() {
    let programs = Programs {
        mafft: "/opt/bin/mafft".into(),
        iqtree: "iqtree2".into(),
        helixloom: "/usr/local/bin/helixloom".into(),
    };
    let snake = snakefile(&pipeline(FULL), &source(), &programs).expect("a Snakefile");
    assert!(snake.contains("shell: \"/opt/bin/mafft --auto"), "{snake}");
    assert!(snake.contains("shell: \"iqtree2 -s"), "{snake}");
    assert!(
        snake.contains("shell: \"/usr/local/bin/helixloom astral"),
        "{snake}"
    );
    assert!(
        snake.contains(
            "# It needs /opt/bin/mafft, iqtree2 and /usr/local/bin/helixloom on the PATH"
        ),
        "{snake}"
    );
    // a name that could add to a command is refused, by the workflows as well
    let bad = Programs {
        iqtree: "iq;tree".into(),
        ..Programs::default()
    };
    assert!(snakefile(&pipeline(FULL), &source(), &bad).is_err());
    assert!(nextflow(&pipeline(FULL), &source(), &bad).is_err());
    let odd_folder = Source {
        folder: "/data/my loci".into(),
        ext: "fasta".into(),
    };
    assert!(snakefile(&pipeline(FULL), &odd_folder, &Programs::default()).is_err());
    assert!(nextflow(&pipeline(FULL), &odd_folder, &Programs::default()).is_err());
    for bad in ["mafft; rm -rf", "a b", "iq\"tree", "", "$(x)", "a`b`"] {
        let said = Programs {
            mafft: bad.into(),
            ..Programs::default()
        }
        .check()
        .unwrap_err();
        assert!(said.contains("program name"), "{bad:?}: {said}");
    }
}

#[test]
fn steps_that_only_this_program_can_do_are_refused_with_what_to_change() {
    let builtin = pipeline("helixloom-pipeline 1\nname: Built in\ninput: sequences\nstep align\nstep tree engine=iqtree\nstep species\n");
    let said = snakefile(&builtin, &source(), &Programs::default()).unwrap_err();
    assert!(
        said.contains("Step 1 (align)")
            && said.contains("engine=mafft")
            && said.contains("pipeline run"),
        "{said}"
    );
    let tree = pipeline(
        "helixloom-pipeline 1\nname: Built in\ninput: alignment\nstep tree\nstep species\n",
    );
    let said = nextflow(&tree, &source(), &Programs::default()).unwrap_err();
    assert!(
        said.contains("Step 1 (tree)") && said.contains("engine=iqtree"),
        "{said}"
    );
    // a pipeline that does not hold together is refused the way a run refuses it, not exported
    let mut broken = pipeline(FULL);
    broken.steps.swap(0, 1);
    assert!(snakefile(&broken, &source(), &Programs::default()).is_err());
}

fn scratch(label: &str) -> std::path::PathBuf {
    crate::pipeline::tools::work_folder(label).expect("a folder")
}

#[test]
fn the_extension_of_the_files_is_found_and_a_folder_that_has_two_is_refused() {
    let dir = scratch("export-source");
    let wanted = |kind| source_of(&dir, kind);
    assert!(wanted(Kind::Alignment).unwrap_err().contains("no file"));
    std::fs::write(dir.join("a.fasta"), ">x\nAC\n").expect("a file");
    std::fs::write(dir.join("b.fasta"), ">x\nAC\n").expect("a file");
    std::fs::write(dir.join("notes.txt"), "not read").expect("a file");
    let found = wanted(Kind::Alignment).expect("a source");
    assert_eq!(found.ext, "fasta");
    assert_eq!(
        found.folder,
        std::fs::canonicalize(&dir)
            .expect("a path")
            .display()
            .to_string()
    );
    std::fs::write(dir.join("c.fa"), ">x\nAC\n").expect("a file");
    let said = wanted(Kind::Alignment).unwrap_err();
    assert!(
        said.contains("fa") && said.contains("fasta") && said.contains("one ending"),
        "{said}"
    );
    // trees are another set of endings
    assert!(wanted(Kind::Trees).unwrap_err().contains("no file"));
    let _ = std::fs::remove_dir_all(&dir);
}

#[test]
fn a_folder_whose_name_could_add_to_a_command_or_to_a_pattern_is_refused() {
    for bad in [
        "/data/a b",
        "/data/{x}",
        "/data/a\"b",
        "/data/a$b",
        "/data/a`b",
        "/data/a*b",
        "/data/a\\b",
        "/data/a[b]",
    ] {
        let said = check_folder(bad).unwrap_err();
        assert!(said.contains("plain"), "{bad}: {said}");
    }
    assert!(check_folder("/data/loci_v2/genes-1.5").is_ok());
}

#[test]
fn the_files_are_written_into_a_folder_and_say_what_they_are() {
    let dir = scratch("export-write");
    let out = dir.join("workflow");
    let files =
        write_folder(&out, &pipeline(FULL), &source(), &Programs::default()).expect("files");
    let names: Vec<String> = files
        .iter()
        .map(|f| {
            f.file_name()
                .expect("a name")
                .to_string_lossy()
                .into_owned()
        })
        .collect();
    assert_eq!(names, ["Snakefile", "main.nf"]);
    assert_eq!(
        std::fs::read_to_string(out.join("Snakefile")).expect("a file"),
        SNAKE_FULL
    );
    assert!(std::fs::read_to_string(out.join("main.nf"))
        .expect("a file")
        .contains("process SPECIES"));
    // nothing is written when the pipeline cannot be exported
    let other = dir.join("none");
    let built_in = Pipeline::unaligned_loci_to_species();
    assert!(write_folder(&other, &built_in, &source(), &Programs::default()).is_err());
    assert!(!other.exists());
    let _ = std::fs::remove_dir_all(&dir);
}

//! The same workflow as files other workflow engines run: a Snakefile for Snakemake and a main.nf for Nextflow,
//! each calling MAFFT, IQ-TREE and primer3 with the settings the window shows, with the input sequences and
//! primer3's input written beside them.

use super::steps::{fasta_text, primer3_input, Named};
use crate::seq::primer3::Settings;
use std::path::{Path, PathBuf};

/// What the exported workflow contains and how its programs are called.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Plan {
    pub align: bool,
    pub tree: bool,
    pub primers: bool,
    /// Ultrafast bootstrap replicates for the tree; 0 for none (IQ-TREE wants at least 1000).
    pub bootstrap: usize,
    pub mafft: String,
    pub iqtree: String,
    pub primer3: String,
}

impl Default for Plan {
    fn default() -> Self {
        Plan {
            align: true,
            tree: true,
            primers: true,
            bootstrap: 0,
            mafft: "mafft".into(),
            iqtree: "iqtree3".into(),
            primer3: "primer3_core".into(),
        }
    }
}

impl Plan {
    /// IQ-TREE's options (what `steps::tree` runs), without the alignment and the output name.
    fn iqtree_options(&self) -> String {
        let bootstrap = if self.bootstrap > 0 {
            format!(" -B {}", self.bootstrap.max(1000))
        } else {
            String::new()
        };
        format!("-m MFP{bootstrap} -seed 1 -redo -quiet")
    }

    /// The files each step makes, in the folder the workflow runs in.
    fn outputs(&self) -> Vec<&'static str> {
        let mut out = Vec::new();
        if self.align {
            out.push("results/aligned.fasta");
        }
        if self.tree {
            out.push("results/tree.treefile");
        }
        if self.primers {
            out.push("results/primers.txt");
        }
        out
    }
}

pub fn snakefile(p: &Plan) -> String {
    let mut t = String::from(
        "# Written by HelixLoom. Run it in this folder with:  snakemake --cores 1\n\
         # It needs the programs named in the shell lines on the PATH.\n\n",
    );
    let targets: Vec<String> = p
        .outputs()
        .iter()
        .map(|o| format!("        \"{o}\","))
        .collect();
    t.push_str(&format!("rule all:\n    input:\n{}\n", targets.join("\n")));
    if p.align || p.tree {
        t.push_str(&format!(
            "\nrule align:\n    input: \"sequences.fasta\"\n    output: \"results/aligned.fasta\"\n    shell: \"{} --auto --quiet {{input}} > {{output}}\"\n",
            p.mafft
        ));
    }
    if p.tree {
        t.push_str(&format!(
            "\nrule tree:\n    input: \"results/aligned.fasta\"\n    output: \"results/tree.treefile\"\n    shell: \"{} -s {{input}} {} --prefix results/tree\"\n",
            p.iqtree,
            p.iqtree_options()
        ));
    }
    if p.primers {
        t.push_str(&format!(
            "\nrule primers:\n    input: \"primer3_input.txt\"\n    output: \"results/primers.txt\"\n    shell: \"{} {{input}} > {{output}}\"\n",
            p.primer3
        ));
    }
    t
}

pub fn nextflow(p: &Plan) -> String {
    let mut t = String::from(
        "// Written by HelixLoom. Run it in this folder with:  nextflow run main.nf\n\
         // It needs the programs named in the script lines on the PATH.\n\
         nextflow.enable.dsl = 2\n\n\
         params.sequences = \"${projectDir}/sequences.fasta\"\n\
         params.primer_input = \"${projectDir}/primer3_input.txt\"\n\
         params.outdir = \"${projectDir}/results\"\n",
    );
    let process = |name: &str, input: &str, output: &str, command: String| {
        format!(
            "\nprocess {name} {{\n    publishDir params.outdir, mode: 'copy'\n    input:\n    path {input}\n    output:\n    path \"{output}\"\n    script:\n    \"\"\"\n    {command}\n    \"\"\"\n}}\n"
        )
    };
    if p.align || p.tree {
        t.push_str(&process(
            "ALIGN",
            "sequences",
            "aligned.fasta",
            format!("{} --auto --quiet ${{sequences}} > aligned.fasta", p.mafft),
        ));
    }
    if p.tree {
        t.push_str(&process(
            "TREE",
            "alignment",
            "tree.treefile",
            format!(
                "{} -s ${{alignment}} {} --prefix tree",
                p.iqtree,
                p.iqtree_options()
            ),
        ));
    }
    if p.primers {
        t.push_str(&process(
            "PRIMERS",
            "spec",
            "primers.txt",
            format!("{} ${{spec}} > primers.txt", p.primer3),
        ));
    }
    t.push_str("\nworkflow {\n");
    if p.align || p.tree {
        t.push_str("    aligned = ALIGN(Channel.fromPath(params.sequences))\n");
    }
    if p.tree {
        t.push_str("    TREE(aligned)\n");
    }
    if p.primers {
        t.push_str("    PRIMERS(Channel.fromPath(params.primer_input))\n");
    }
    t.push_str("}\n");
    t
}

/// What primer3 is asked: the sequence, the stretch (first base from 0, length) a product must contain, and
/// the settings.
pub struct PrimerSpec<'a> {
    pub sequence: &'a Named,
    pub target: Option<(usize, usize)>,
    pub settings: &'a Settings,
}

/// Write the workflow folder: the sequences, primer3's input, the Snakefile and main.nf. Returns the files made.
pub fn write_folder(
    dir: &Path,
    seqs: &[Named],
    plan: &Plan,
    primer: Option<PrimerSpec>,
) -> Result<Vec<PathBuf>, String> {
    let mut files = vec![
        ("sequences.fasta", fasta_text(&clean_names(seqs))),
        ("Snakefile", snakefile(plan)),
        ("main.nf", nextflow(plan)),
    ];
    if plan.primers {
        let spec = primer.ok_or("The primer step needs a sequence.")?;
        files.push((
            "primer3_input.txt",
            primer3_input(
                &clean_names(std::slice::from_ref(spec.sequence))[0].0,
                &spec.sequence.1,
                spec.target,
                spec.settings,
            ),
        ));
    }
    std::fs::create_dir_all(dir).map_err(|e| format!("Could not make {}: {e}", dir.display()))?;
    files
        .into_iter()
        .map(|(name, text)| {
            let path = dir.join(name);
            std::fs::write(&path, text)
                .map(|()| path.clone())
                .map_err(|e| format!("Could not write {}: {e}", path.display()))
        })
        .collect()
}

/// Names without spaces or characters that FASTA tools and Newick trees treat specially.
fn clean_names(seqs: &[Named]) -> Vec<Named> {
    seqs.iter()
        .map(|(name, bases)| {
            let safe: String = name
                .chars()
                .map(|c| {
                    if c.is_ascii_alphanumeric() || "-_.".contains(c) {
                        c
                    } else {
                        '_'
                    }
                })
                .collect();
            (safe, bases.clone())
        })
        .collect()
}

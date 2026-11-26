//! A pipeline file as the workflow files other engines run: a `Snakefile` for Snakemake and a `main.nf` for Nextflow, each calling
//! MAFFT, IQ-TREE or a HelixLoom command with the settings the engine uses. Export requires a command for each step.

use super::{accepted_files, extensions, Kind, Pipeline, Step};
use std::path::{Path, PathBuf};

/// How the programs are named in the commands.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Programs {
    pub mafft: String,
    pub iqtree: String,
    pub helixloom: String,
}

impl Default for Programs {
    fn default() -> Self {
        Programs {
            mafft: "mafft".into(),
            iqtree: "iqtree3".into(),
            helixloom: "helixloom".into(),
        }
    }
}

/// Whether `text` is made of characters that mean nothing to a shell, a Python string or a glob.
fn plain(text: &str) -> bool {
    !text.is_empty()
        && text
            .chars()
            .all(|c| c.is_ascii_alphanumeric() || "/._+-~@,:=".contains(c))
}

impl Programs {
    /// A name goes into a command line, so it has to be plain.
    pub fn check(&self) -> Result<(), String> {
        for name in [&self.mafft, &self.iqtree, &self.helixloom] {
            if !plain(name) {
                return Err(format!(
                    "The program name {name:?} is not plain: use letters, digits and / . _ + - only."
                ));
            }
        }
        Ok(())
    }
}

/// Where the workflow reads its files from: the folder, as written into the workflow, and the one ending the files have.
#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Source {
    pub folder: String,
    pub ext: String,
}

/// A folder name goes into a pattern and a string of the workflow, so it has to be plain.
pub fn check_folder(path: &str) -> Result<(), String> {
    if plain(path) {
        Ok(())
    } else {
        Err(format!(
            "The folder {path:?} is not plain enough to write into a workflow: it may hold letters, digits and / . _ + - ~ @ , : = only (no space, quote, brace, bracket, star, dollar or backslash). Move or rename it."
        ))
    }
}

/// The folder (as an absolute path) and the ending of the files in it that a pipeline taking `kind` reads. A workflow finds its files
/// by one ending, so a folder that has two is refused.
pub fn source_of(dir: &Path, kind: Kind) -> Result<Source, String> {
    let files =
        accepted_files(dir, kind).map_err(|e| format!("Could not read {}: {e}", dir.display()))?;
    let mut ends: Vec<String> = files
        .iter()
        .filter_map(|p| p.extension().and_then(|x| x.to_str()).map(str::to_string))
        .collect();
    ends.sort();
    ends.dedup();
    let ext = match ends.as_slice() {
        [] => {
            return Err(format!(
                "There is no file in {} with one of the endings {}.",
                dir.display(),
                extensions(kind).join(", ")
            ))
        }
        [one] => one.clone(),
        many => {
            return Err(format!(
                "The files of {} have more than one ending ({}): a workflow finds its files by one ending, so give them one ending.",
                dir.display(),
                many.join(", ")
            ))
        }
    };
    let folder = std::fs::canonicalize(dir)
        .map_err(|e| format!("Could not find {}: {e}", dir.display()))?
        .display()
        .to_string();
    check_folder(&folder)?;
    Ok(Source { folder, ext })
}

/// What is done for each file, in order, and what is done with the trees after.
struct Plan {
    each: Vec<Each>,
    species: Option<Species>,
}

#[derive(Clone, Copy)]
enum Each {
    Align,
    Tree { bootstrap: usize },
    GfaInspect,
}

#[derive(Clone, Copy, PartialEq)]
enum Species {
    /// The gene trees of all the files are put together and the species tree is made once.
    Gathered,
    /// Every file holds trees already, and gets a species tree of its own.
    PerFile,
}

impl Plan {
    fn of(p: &Pipeline, programs: &Programs) -> Result<Plan, String> {
        programs.check()?;
        p.check()?;
        for (k, step) in p.steps.iter().enumerate() {
            let (name, fix) = match step {
                Step::Align { mafft: false } => ("align", "engine=mafft"),
                Step::Tree { iqtree: false, .. } => ("tree", "engine=iqtree"),
                Step::Asm(a) => {
                    return Err(format!(
                        "Step {} ({}) is a step of the assembly routes, which `helixloom pipeline run` does and no Snakefile or main.nf is written for.",
                        k + 1,
                        a.name()
                    ))
                }
                _ => continue,
            };
            return Err(format!(
                "Step {} ({name}) is done by this program's own code, which no other workflow engine can run: use {fix} to export it, or run the pipeline itself with `helixloom pipeline run`.",
                k + 1
            ));
        }
        let each = p
            .steps
            .iter()
            .filter_map(|s| match s {
                Step::Align { .. } => Some(Each::Align),
                Step::Tree { bootstrap, .. } => Some(Each::Tree {
                    bootstrap: *bootstrap,
                }),
                Step::Species => None,
                Step::GfaInspect => Some(Each::GfaInspect),
                // refused above
                Step::Asm(_) => None,
            })
            .collect();
        let species = p.steps.iter().any(|s| matches!(s, Step::Species)).then(|| {
            if p.input == Kind::Trees {
                Species::PerFile
            } else {
                Species::Gathered
            }
        });
        Ok(Plan { each, species })
    }

    fn has_tree(&self) -> bool {
        self.each.iter().any(|e| matches!(e, Each::Tree { .. }))
    }

    /// The programs the workflow calls, in the order they are used.
    fn needs(&self, programs: &Programs) -> String {
        let mut names = Vec::new();
        if self.each.iter().any(|e| matches!(e, Each::Align)) {
            names.push(programs.mafft.as_str());
        }
        if self.has_tree() {
            names.push(programs.iqtree.as_str());
        }
        if self.species.is_some() || self.each.iter().any(|e| matches!(e, Each::GfaInspect)) {
            names.push(programs.helixloom.as_str());
        }
        match names.as_slice() {
            [] => String::new(),
            [one] => (*one).to_string(),
            [init @ .., last] => format!("{} and {last}", init.join(", ")),
        }
    }

    fn header(&self, comment: &str, name: &str, run: &str, programs: &Programs) -> String {
        let mut text = format!(
            "{comment} Written by HelixLoom from the pipeline \"{name}\". Run it with:  {run}\n{comment} It needs {} on the PATH.",
            self.needs(programs)
        );
        if self.has_tree() {
            text.push_str(" Sequence names must be plain (letters, digits, - _ .) for IQ-TREE.");
        }
        text.push('\n');
        text
    }
}

/// IQ-TREE's options after the alignment, as the engine runs it (`steps::tree`).
fn iqtree_options(bootstrap: usize) -> String {
    let bootstrap = if bootstrap > 0 {
        // (a pipeline that asks for 1 to 999 replicates is refused before this, as IQ-TREE makes at least 1000)
        format!(" -B {bootstrap}")
    } else {
        String::new()
    };
    format!("-m MFP{bootstrap} -seed 1 -redo -quiet")
}

pub fn snakefile(p: &Pipeline, src: &Source, programs: &Programs) -> Result<String, String> {
    let plan = Plan::of(p, programs)?;
    check_folder(&src.folder)?;
    let mut t = plan.header("#", &p.name, "snakemake --cores 4", programs);
    t.push_str(&format!(
        "\nSAMPLES, = glob_wildcards(\"{}/{{sample}}.{}\")\n\nwildcard_constraints:\n    sample=\"[^/]+\"\n",
        src.folder, src.ext
    ));
    // the file each step reads and the one it writes, by rule
    let mut rules = String::new();
    let mut now = format!("{}/{{sample}}.{}", src.folder, src.ext);
    for e in &plan.each {
        let (rule, out, shell) = match e {
            Each::Align => (
                "align",
                "results/aligned/{sample}.fasta".to_string(),
                format!("{} --auto --quiet {{input}} > {{output}}", programs.mafft),
            ),
            Each::GfaInspect => (
                "gfa_inspect",
                "results/summaries/{sample}.tsv".to_string(),
                format!(
                    "{} gfa inspect {{input:q}} > {{output:q}}",
                    programs.helixloom
                ),
            ),
            Each::Tree { bootstrap } => (
                "tree",
                "results/trees/{sample}.treefile".to_string(),
                format!(
                    "{} -s {{input}} {} --prefix results/trees/{{wildcards.sample}}",
                    programs.iqtree,
                    iqtree_options(*bootstrap)
                ),
            ),
        };
        rules.push_str(&format!(
            "\nrule {rule}:\n    input: \"{now}\"\n    output: \"{out}\"\n    shell: \"{shell}\"\n"
        ));
        now = out;
    }
    let target = match plan.species {
        Some(Species::Gathered) => {
            rules.push_str(&format!(
                "\nrule gather:\n    input: expand(\"{now}\", sample=SAMPLES)\n    output: \"results/gene_trees.tre\"\n    shell: \"cat {{input}} > {{output}}\"\n"
            ));
            rules.push_str(&format!(
                "\nrule species:\n    input: \"results/gene_trees.tre\"\n    output: \"results/species.tre\"\n    shell: \"{} astral -i {{input}} -m 4 -o {{output}}\"\n",
                programs.helixloom
            ));
            "\"results/species.tre\"".to_string()
        }
        Some(Species::PerFile) => {
            rules.push_str(&format!(
                "\nrule species:\n    input: \"{now}\"\n    output: \"results/species/{{sample}}.tre\"\n    shell: \"{} astral -i {{input}} -o {{output}}\"\n",
                programs.helixloom
            ));
            "expand(\"results/species/{sample}.tre\", sample=SAMPLES)".to_string()
        }
        None => format!("expand(\"{now}\", sample=SAMPLES)"),
    };
    t.push_str(&format!("\nrule all:\n    input: {target}\n"));
    t.push_str(&rules);
    Ok(t)
}

pub fn nextflow(p: &Pipeline, src: &Source, programs: &Programs) -> Result<String, String> {
    let plan = Plan::of(p, programs)?;
    check_folder(&src.folder)?;
    let mut t = plan.header("//", &p.name, "nextflow run main.nf", programs);
    t.push_str(&format!(
        "nextflow.enable.dsl = 2\n\nparams.input = \"{}/*.{}\"\nparams.outdir = \"${{projectDir}}/results\"\n",
        src.folder, src.ext
    ));
    let process = |name: &str, publish: &str, input: &str, output: &str, command: &str| {
        format!(
            "\nprocess {name} {{\n    publishDir {publish}, mode: 'copy'\n    input:\n    path {input}\n    output:\n    path \"{output}\"\n    script:\n    \"\"\"\n{command}\n    \"\"\"\n}}\n"
        )
    };
    let mut channel = if p.input == Kind::Gfa {
        "graphs"
    } else {
        "loci"
    }
    .to_string();
    let mut calls = format!("    {channel} = Channel.fromPath(params.input)\n");
    for e in &plan.each {
        match e {
            Each::GfaInspect => {
                // Keep sample names out of shell commands. Nextflow stages a fixed input name
                // and assigns the per-sample name when publishing the completed summary.
                t.push_str(&format!(concat!(
                    "\nprocess GFA_INSPECT {{\n",
                    "    publishDir \"${{params.outdir}}/summaries\", mode: 'copy', saveAs: {{ \"${{sample}}.tsv\" }}\n",
                    "    input:\n    tuple val(sample), path('graph.gfa')\n",
                    "    output:\n    path 'summary.tsv'\n",
                    "    script:\n    \"\"\"\n    {} gfa inspect graph.gfa > summary.tsv\n    \"\"\"\n}}\n"
                ), programs.helixloom));
                calls.push_str(&format!(
                    "    GFA_INSPECT({channel}.map {{ graph -> tuple(graph.baseName, graph) }})\n"
                ));
            }
            Each::Align => {
                t.push_str(&process(
                    "ALIGN",
                    "\"${params.outdir}/aligned\"",
                    "locus",
                    "${locus.baseName}.aligned.fasta",
                    &format!(
                        "    {} --auto --quiet \"${{locus}}\" > \"${{locus.baseName}}.aligned.fasta\"",
                        programs.mafft
                    ),
                ));
                calls.push_str(&format!("    aligned = ALIGN({channel})\n"));
                channel = "aligned".into();
            }
            Each::Tree { bootstrap } => {
                t.push_str(&process(
                    "TREE",
                    "\"${params.outdir}/trees\"",
                    "alignment",
                    "${alignment.baseName}.treefile",
                    &format!(
                        "    {} -s \"${{alignment}}\" {} --prefix \"${{alignment.baseName}}\"",
                        programs.iqtree,
                        iqtree_options(*bootstrap)
                    ),
                ));
                calls.push_str(&format!("    gene_trees = TREE({channel})\n"));
                channel = "gene_trees".into();
            }
        }
    }
    match plan.species {
        Some(Species::Gathered) => {
            t.push_str(&process(
                "SPECIES",
                "params.outdir",
                "trees",
                "species.tre",
                &format!(
                    "    cat ${{trees}} > gene_trees.tre\n    {} astral -i gene_trees.tre -m 4 -o species.tre",
                    programs.helixloom
                ),
            ));
            calls.push_str(&format!("    SPECIES({channel}.collect())\n"));
        }
        Some(Species::PerFile) => {
            t.push_str(&process(
                "SPECIES",
                "\"${params.outdir}/species\"",
                "trees",
                "${trees.baseName}.species.tre",
                &format!(
                    "    {} astral -i \"${{trees}}\" -o \"${{trees.baseName}}.species.tre\"",
                    programs.helixloom
                ),
            ));
            calls.push_str(&format!("    SPECIES({channel})\n"));
        }
        None => {}
    }
    t.push_str(&format!("\nworkflow {{\n{calls}}}\n"));
    Ok(t)
}

/// Write the `Snakefile` and `main.nf` of `p` into `dir` (made if it is not there); the paths written. Nothing is written when the
/// pipeline cannot be exported.
pub fn write_folder(
    dir: &Path,
    p: &Pipeline,
    src: &Source,
    programs: &Programs,
) -> Result<Vec<PathBuf>, String> {
    let files = [
        ("Snakefile", snakefile(p, src, programs)?),
        ("main.nf", nextflow(p, src, programs)?),
    ];
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

#[cfg(test)]
mod tests;

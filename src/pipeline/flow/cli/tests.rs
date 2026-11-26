use super::*;
use crate::pipeline::flow::tests::locus;
use std::sync::atomic::{AtomicU32, Ordering};

static NEXT: AtomicU32 = AtomicU32::new(0);

struct Folder(PathBuf);

impl Folder {
    fn new() -> Folder {
        let path = std::env::temp_dir().join(format!(
            "helixloom-flowcli-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&path).expect("a folder");
        Folder(path)
    }

    fn path(&self, name: &str) -> String {
        self.0.join(name).to_string_lossy().into_owned()
    }

    fn write(&self, name: &str, text: &str) {
        std::fs::write(self.0.join(name), text).expect("a file");
    }

    /// A folder of loci and a pipeline file for them.
    fn loci(&self, n: u64) {
        std::fs::create_dir_all(self.0.join("loci")).unwrap();
        for k in 1..=n {
            let text = Value::Alignment(locus(k * 17, 150)).to_text();
            self.write(&format!("loci/gene{k}.fasta"), &text);
        }
        self.write("p.txt", &Pipeline::loci_to_species().to_text());
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

fn words(s: &str) -> Vec<String> {
    s.split_whitespace().map(str::to_string).collect()
}

fn say(args: &[String]) -> (String, bool) {
    work(args).expect("the command worked")
}

#[test]
fn the_examples_are_files_that_read_back() {
    let (text, ok) = say(&words("example"));
    assert!(ok);
    assert_eq!(Pipeline::parse(&text), Ok(Pipeline::loci_to_species()));
    let (text, _) = say(&words("example unaligned"));
    assert_eq!(
        Pipeline::parse(&text),
        Ok(Pipeline::unaligned_loci_to_species())
    );
    let (text, _) = say(&words("example gfa"));
    assert_eq!(Pipeline::parse(&text), Ok(Pipeline::inspect_graphs()));
    assert!(work(&words("example nonsense")).is_err());
}

#[test]
fn check_says_what_the_pipeline_does() {
    let f = Folder::new();
    f.write("p.txt", &Pipeline::unaligned_loci_to_species().to_text());
    let (text, ok) = say(&["check".to_string(), f.path("p.txt")]);
    assert!(ok);
    assert!(
        text.contains("1. align") && text.contains("3. species"),
        "{text}"
    );
    assert!(
        text.contains("Steps 1 to 2 are done for each file"),
        "{text}"
    );
    assert!(text.contains("ends with tree"), "{text}");
    f.write(
        "bad.txt",
        "helixloom-pipeline 1\ninput: alignment\nstep align\n",
    );
    let why = work(&["check".to_string(), f.path("bad.txt")]).unwrap_err();
    assert!(
        why.contains("bad.txt") && why.contains("Step 1 (align) takes sequences"),
        "{why}"
    );
}

#[test]
fn a_folder_is_run_written_out_and_the_second_run_finds_it_all_in_the_cache() {
    let f = Folder::new();
    f.loci(3);
    let args = |extra: &str| {
        let mut a = vec!["run".to_string(), f.path("p.txt"), f.path("loci")];
        a.extend(words(extra).into_iter().map(|w| match w.strip_prefix('@') {
            Some(name) => f.path(name),
            None => w,
        }));
        a
    };
    let flags = "--cache @cache --out @species.tre --out-dir @trees";
    let (first, ok) = say(&args(flags));
    assert!(ok, "{first}");
    assert!(
        first.contains("4 steps run, 0 found in the cache"),
        "{first}"
    );
    assert!(first.contains("gene1: tree: ran in"), "{first}");
    let species = std::fs::read_to_string(f.path("species.tre")).unwrap();
    assert!(
        species.contains("taxon1") && species.trim_end().ends_with(';'),
        "{species}"
    );
    for k in 1..=3 {
        assert!(
            std::fs::read_to_string(f.path(&format!("trees/gene{k}.tre")))
                .unwrap()
                .contains("taxon")
        );
    }
    let (second, ok) = say(&args(flags));
    assert!(ok);
    assert!(
        second.contains("0 steps run, 4 found in the cache"),
        "{second}"
    );
    assert_eq!(
        species,
        std::fs::read_to_string(f.path("species.tre")).unwrap()
    );
}

#[test]
fn a_run_with_a_failed_step_says_so_writes_what_it_can_and_is_not_ok() {
    let f = Folder::new();
    f.loci(3);
    f.write("loci/broken.fasta", ">a\nACGTACGT\n>b\nACGTACGA\n");
    let (text, ok) = say(&["run".to_string(), f.path("p.txt"), f.path("loci")]);
    assert!(!ok);
    assert!(text.contains("broken: tree: FAILED"), "{text}");
    assert!(text.contains("1 failed; left out: broken"), "{text}");
}

#[test]
fn one_file_may_be_given_instead_of_a_folder() {
    let f = Folder::new();
    f.loci(1);
    f.write(
        "one.txt",
        "helixloom-pipeline 1\ninput: alignment\nstep tree\n",
    );
    let (text, ok) = say(&[
        "run".to_string(),
        f.path("one.txt"),
        f.path("loci/gene1.fasta"),
        "--out".to_string(),
        f.path("t.tre"),
    ]);
    assert!(ok, "{text}");
    assert!(text.contains("gene1: tree: ran"), "{text}");
    assert!(std::fs::read_to_string(f.path("t.tre"))
        .unwrap()
        .contains("taxon"));
}

#[test]
fn clearing_the_cache_forgets_what_it_holds() {
    let f = Folder::new();
    f.loci(2);
    say(&[
        "run".to_string(),
        f.path("p.txt"),
        f.path("loci"),
        "--cache".to_string(),
        f.path("c"),
    ]);
    let (text, _) = say(&["clear-cache".to_string(), f.path("c")]);
    assert!(text.contains("Forgot 3 answers"), "{text}");
}

#[test]
fn a_command_that_cannot_be_understood_says_why_and_no_command_gives_the_usage() {
    assert!(work(&words("run onlyone"))
        .unwrap_err()
        .contains("pipeline file and an input"));
    assert!(work(&words("run a b --colour red"))
        .unwrap_err()
        .contains("no option --colour"));
    assert!(work(&words("run a b --cache"))
        .unwrap_err()
        .contains("needs a value"));
    assert!(work(&words("fly"))
        .unwrap_err()
        .contains("no pipeline command called fly"));
    assert!(work(&words("check")).is_err());
    let (text, ok) = say(&[]);
    assert!(ok && text.contains("helixloom pipeline run PIPELINE INPUT"));
    assert!(work(&words("run nowhere.txt x"))
        .unwrap_err()
        .contains("Could not read"));
}

const EXTERNAL: &str = "helixloom-pipeline 1\nname: Gene trees, then a species tree\ninput: alignment\nstep tree engine=iqtree\nstep species\n";

#[test]
fn a_pipeline_is_exported_as_two_workflows_that_read_the_folder() {
    let folder = Folder::new();
    folder.loci(3);
    folder.write("ext.txt", EXTERNAL);
    let (text, ok) = say(&words(&format!(
        "export {} {} --out-dir {} --iqtree iqtree2",
        folder.path("ext.txt"),
        folder.path("loci"),
        folder.path("wf")
    )));
    assert!(ok);
    assert!(
        text.contains("Snakefile") && text.contains("main.nf") && text.contains("*.fasta"),
        "{text}"
    );
    let snake = std::fs::read_to_string(folder.0.join("wf/Snakefile")).expect("a Snakefile");
    let canonical = std::fs::canonicalize(folder.0.join("loci")).expect("a path");
    assert!(
        snake.contains(&format!(
            "glob_wildcards(\"{}/{{sample}}.fasta\")",
            canonical.display()
        )),
        "{snake}"
    );
    assert!(snake.contains("shell: \"iqtree2 -s {input}"), "{snake}");
    assert!(folder.0.join("wf/main.nf").is_file());
    // without --out-dir the folder is called workflow (here, beside where the test runs: not written, so use the error paths)
    let (usage, _) = say(&[]);
    assert!(
        usage.contains("helixloom pipeline export PIPELINE INPUT"),
        "{usage}"
    );
}

#[test]
fn an_export_that_cannot_be_made_says_why_and_writes_nothing() {
    let folder = Folder::new();
    folder.loci(2);
    // the example pipeline is done by this program's own code
    let built_in = format!(
        "export {} {} --out-dir {}",
        folder.path("p.txt"),
        folder.path("loci"),
        folder.path("wf")
    );
    let said = work(&words(&built_in)).unwrap_err();
    assert!(
        said.contains("engine=iqtree") && said.contains("pipeline run"),
        "{said}"
    );
    assert!(!folder.0.join("wf").exists());
    folder.write("ext.txt", EXTERNAL);
    // a file is not a folder; two endings are refused; the options are checked
    let file = work(&words(&format!(
        "export {} {}",
        folder.path("ext.txt"),
        folder.path("p.txt")
    )))
    .unwrap_err();
    assert!(file.contains("not a folder"), "{file}");
    folder.write("loci/extra.fa", ">a\nAC\n");
    let mixed = work(&words(&format!(
        "export {} {}",
        folder.path("ext.txt"),
        folder.path("loci")
    )))
    .unwrap_err();
    assert!(mixed.contains("more than one ending"), "{mixed}");
    assert!(work(&words("export a"))
        .unwrap_err()
        .contains("needs a pipeline file and the folder"));
    assert!(work(&words("export a b --colour red"))
        .unwrap_err()
        .contains("no option --colour"));
    assert!(work(&words("export a b --mafft"))
        .unwrap_err()
        .contains("needs a value"));
    std::fs::remove_file(folder.0.join("loci/extra.fa")).expect("removed");
    let bad_name = work(&words(&format!(
        "export {} {} --iqtree a;b --out-dir {}",
        folder.path("ext.txt"),
        folder.path("loci"),
        folder.path("wf")
    )))
    .unwrap_err();
    assert!(bad_name.contains("program name"), "{bad_name}");
    assert!(!folder.0.join("wf").exists());
}

#[test]
fn the_routes_are_examples_that_read_back_and_an_unknown_name_lists_them() {
    for (name, _) in crate::pipeline::asm::routes::NAMES {
        let (text, ok) = say(&words(&format!("example {name}")));
        assert!(ok, "{name}");
        let p = Pipeline::parse(&text).unwrap_or_else(|e| panic!("{name}: {e}"));
        assert_eq!(Some(p), crate::pipeline::asm::routes::named(name), "{name}");
    }
    let e = work(&words("example nothing")).unwrap_err();
    assert!(
        e.contains("plasmid") && e.contains("nanopore-plasmid") && e.contains("minion-reads"),
        "{e}"
    );
    assert!(work(&words("example plasmid more"))
        .unwrap_err()
        .contains("takes one name"));
}

#[test]
fn check_says_what_a_route_does_step_by_step() {
    let f = Folder::new();
    f.write("p.txt", &crate::pipeline::asm::routes::plasmid().to_text());
    let (text, ok) = say(&["check".to_string(), f.path("p.txt")]);
    assert!(ok);
    assert!(text.contains("input: reads"), "{text}");
    assert!(
        text.contains(
            "1. trim adapter=AGATCGGAAGAGC quality=20 min-length=20 cut=tail  (reads to reads)"
        ),
        "{text}"
    );
    assert!(
        text.contains("2. assemble ks=auto cov-cutoff=auto min-qual=10  (reads to gfa)"),
        "{text}"
    );
    assert!(text.contains("3. hybrid long=long-reads"), "{text}");
    assert!(text.ends_with("ends with sequences.\n"), "{text}");
}

#[test]
fn a_folder_of_long_reads_is_selected_and_written_as_fastq_named_after_each_file() {
    let f = Folder::new();
    std::fs::create_dir_all(f.0.join("reads")).unwrap();
    let read = |len: usize, q: char| {
        let seq = "ACGT".repeat(len / 4);
        format!("@r\n{seq}\n+\n{}\n", q.to_string().repeat(seq.len()))
    };
    f.write(
        "reads/barcode01.fastq",
        &format!("{}{}", read(3000, 'I'), read(500, 'I')),
    );
    f.write("reads/barcode02.fastq", &read(2400, '5'));
    f.write(
        "p.txt",
        &crate::pipeline::asm::routes::minion_reads().to_text(),
    );
    let (text, ok) = say(&[
        "run".to_string(),
        f.path("p.txt"),
        f.path("reads"),
        "--out-dir".to_string(),
        f.path("out"),
    ]);
    assert!(ok, "{text}");
    let one = std::fs::read_to_string(f.0.join("out/barcode01.fastq")).unwrap();
    // the short read is gone, the long one is there
    assert_eq!(one.lines().count(), 4);
    assert!(std::fs::read_to_string(f.0.join("out/barcode02.fastq"))
        .unwrap()
        .starts_with("@r\n"));
}

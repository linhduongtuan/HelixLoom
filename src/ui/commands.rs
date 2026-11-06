//! Every command of the program in one list, for the command palette and the sheet of shortcuts: its group, its name, the
//! key that runs it, other words to find it by, and the `Action` that does it. The menus keep their own lists, as a menu is
//! drawn with its own gaps and sub-menus; the tests below check the keys written here against the keys the program reads
//! (`ui/app/tests.rs`) and that nothing is listed twice.

use super::app::Action;
use super::platform::Workspace;
use super::state::{ThemeMode, TransMode};

#[derive(Clone, Debug)]
pub struct Command {
    pub group: &'static str,
    pub title: &'static str,
    /// The shortcut as the menus write it (`Cmd+Shift+S`); empty for none.
    pub keys: &'static str,
    /// Other words that should find it.
    pub words: &'static str,
    pub action: Action,
}

const fn c(
    group: &'static str,
    title: &'static str,
    keys: &'static str,
    words: &'static str,
    action: Action,
) -> Command {
    Command {
        group,
        title,
        keys,
        words,
        action,
    }
}

/// The commands, in the order the palette lists them when nothing has been typed: the workspaces first, with the keys of the
/// rail (Cmd+Alt+1 and so on).
const COMMANDS: &[Command] = &[
    c(
        "Go to",
        "Home",
        "Cmd+Alt+1",
        "start dashboard overview",
        Action::Go(Workspace::Home),
    ),
    c(
        "Go to",
        "Sequences",
        "Cmd+Alt+2",
        "plasmid editor dna map enzymes primers",
        Action::Go(Workspace::Sequences),
    ),
    c(
        "Go to",
        "Signals",
        "Cmd+Alt+3",
        "sanger chromatogram trace ab1 reads quality",
        Action::Go(Workspace::Signals),
    ),
    c(
        "Go to",
        "Alignment",
        "Cmd+Alt+4",
        "msa mafft align sequences conservation",
        Action::Go(Workspace::Alignment),
    ),
    c(
        "Go to",
        "Phylogeny",
        "Cmd+Alt+5",
        "tree estimate bootstrap likelihood bayesian distance",
        Action::Go(Workspace::Phylogeny),
    ),
    c(
        "Go to",
        "Trees",
        "Cmd+Alt+6",
        "newick draw phylogram cladogram circle",
        Action::Go(Workspace::Trees),
    ),
    c(
        "Go to",
        "Structure",
        "Cmd+Alt+7",
        "protein 3d pdb esmfold fold ribbon",
        Action::Go(Workspace::Structure),
    ),
    c(
        "Go to",
        "Pipelines",
        "Cmd+Alt+8",
        "workflow snakemake nextflow steps run",
        Action::Go(Workspace::Pipelines),
    ),
    c(
        "Go to",
        "Workstation",
        "Cmd+Alt+9",
        "gpu server remote compute jobs",
        Action::Go(Workspace::Workstation),
    ),
    c(
        "Examples",
        "Try the alignment example",
        "",
        "demo sample",
        Action::Example(Workspace::Alignment),
    ),
    c(
        "Examples",
        "Try the phylogeny example",
        "",
        "demo sample",
        Action::Example(Workspace::Phylogeny),
    ),
    c(
        "Examples",
        "Try the tree example",
        "",
        "demo sample",
        Action::Example(Workspace::Trees),
    ),
    c(
        "Examples",
        "Try the Sanger trace example",
        "",
        "demo sample",
        Action::Example(Workspace::Signals),
    ),
    c(
        "Examples",
        "Try the protein structure example",
        "",
        "demo sample",
        Action::Example(Workspace::Structure),
    ),
    c(
        "Examples",
        "Try the pipeline example",
        "",
        "demo sample",
        Action::Example(Workspace::Pipelines),
    ),
    c("File", "New sequence", "Cmd+N", "create", Action::New),
    c(
        "File",
        "Open a file…",
        "Cmd+O",
        "import load fasta genbank",
        Action::Open,
    ),
    c(
        "File",
        "Duplicate selection",
        "",
        "copy new sequence part",
        Action::DuplicateSelection,
    ),
    c(
        "File",
        "Duplicate full sequence",
        "",
        "copy new sequence whole clone",
        Action::DuplicateFull,
    ),
    c("File", "Save", "Cmd+S", "write", Action::Save),
    c(
        "File",
        "Save as…",
        "Cmd+Shift+S",
        "write export",
        Action::SaveAs,
    ),
    c("File", "Export as FASTA…", "", "write", Action::ExportFasta),
    c("File", "Close", "Cmd+W", "", Action::CloseCurrent),
    c("Edit", "Undo", "Cmd+Z", "", Action::Undo),
    c("Edit", "Redo", "Cmd+Shift+Z", "", Action::Redo),
    c("Edit", "Cut", "Cmd+X", "", Action::Cut),
    c("Edit", "Copy", "Cmd+C", "clipboard", Action::Copy),
    c("Edit", "Paste", "Cmd+V", "clipboard", Action::Paste),
    c(
        "Edit",
        "Cut reverse-complement",
        "",
        "revcom clipboard",
        Action::CutRevCom,
    ),
    c(
        "Edit",
        "Copy reverse-complement",
        "",
        "revcom clipboard",
        Action::CopyRevCom,
    ),
    c(
        "Edit",
        "Paste reverse-complement",
        "",
        "revcom clipboard",
        Action::PasteRevCom,
    ),
    c("Edit", "Select all", "Cmd+A", "", Action::SelectAll),
    c(
        "Edit",
        "Select from-to…",
        "Cmd+Shift+A",
        "range",
        Action::SelectRange,
    ),
    c("Sequence", "Find…", "Cmd+F", "search motif", Action::Find),
    c("Sequence", "Find again", "Cmd+G", "next", Action::FindAgain),
    c(
        "Sequence",
        "Clear find highlighting",
        "Cmd+Shift+F",
        "search",
        Action::ClearFind,
    ),
    c(
        "Sequence",
        "Jump to a position…",
        "Cmd+J",
        "go to base",
        Action::Jump,
    ),
    c(
        "Sequence",
        "Reverse-complement",
        "Cmd+/",
        "revcom",
        Action::RevCom,
    ),
    c(
        "Sequence",
        "Convert to UPPERCASE",
        "",
        "case",
        Action::Upper,
    ),
    c(
        "Sequence",
        "Convert to lowercase",
        "",
        "case",
        Action::Lower,
    ),
    c(
        "Sequence",
        "Swap upper and lower case",
        "",
        "case",
        Action::SwapCase,
    ),
    c(
        "Sequence",
        "Switch circular and linear",
        "",
        "topology plasmid",
        Action::ToggleCircular,
    ),
    c(
        "Sequence",
        "Set origin",
        "Cmd+I",
        "rotate plasmid",
        Action::SetOrigin,
    ),
    c(
        "Sequence",
        "Linearize at the caret",
        "",
        "plasmid cut open",
        Action::Linearize,
    ),
    c(
        "Analyse",
        "Restriction enzymes…",
        "Cmd+E",
        "enzyme selector sites cut",
        Action::EnzymeSelector,
    ),
    c(
        "Analyse",
        "Digest…",
        "Cmd+R",
        "fragments gel restriction",
        Action::Digest,
    ),
    c(
        "Analyse",
        "Highlight restriction sites",
        "",
        "enzyme sites",
        Action::ToggleHighlight,
    ),
    c(
        "Analyse",
        "Silent restriction sites…",
        "",
        "silent mutation codon protein make site enzyme",
        Action::SilentSites,
    ),
    c(
        "Analyse",
        "Diagnostic restriction sites…",
        "",
        "mismatch primer mutation screen site enzyme",
        Action::DiagnosticSites,
    ),
    c(
        "Analyse",
        "ORF map and list…",
        "",
        "open reading frame genes",
        Action::OrfList,
    ),
    c(
        "Analyse",
        "Next ORF",
        "Cmd+Shift+.",
        "open reading frame",
        Action::OrfNext,
    ),
    c(
        "Analyse",
        "Previous ORF",
        "Cmd+Shift+,",
        "open reading frame",
        Action::OrfPrev,
    ),
    c(
        "Analyse",
        "Translate…",
        "Cmd+T",
        "protein amino acids",
        Action::Translate,
    ),
    c(
        "Analyse",
        "Align sequences…",
        "Cmd+L",
        "pairwise",
        Action::Align,
    ),
    c(
        "Analyse",
        "PCR products and primers…",
        "",
        "amplify",
        Action::Pcr,
    ),
    c(
        "Analyse",
        "Primer3 pair design…",
        "",
        "primers design",
        Action::Primer3,
    ),
    c(
        "Analyse",
        "Read quality and trimming…",
        "",
        "fastq reads clean",
        Action::Reads,
    ),
    c(
        "Features",
        "New feature…",
        "Cmd+.",
        "annotate",
        Action::NewFeature,
    ),
    c(
        "Features",
        "Edit the chosen feature…",
        "",
        "annotate",
        Action::EditFeature,
    ),
    c(
        "Features",
        "Delete the chosen feature",
        "",
        "annotate remove",
        Action::DeleteFeature,
    ),
    c(
        "Features",
        "Clear all features",
        "",
        "annotate remove",
        Action::ClearFeatures,
    ),
    c(
        "Features",
        "Features from the found matches",
        "",
        "find annotate",
        Action::FeaturesFromFound,
    ),
    c(
        "Features",
        "Annotate using the library",
        "Cmd+Shift+K",
        "features automatic",
        Action::Annotate,
    ),
    c(
        "Features",
        "Open a feature library and annotate…",
        "",
        "features automatic",
        Action::OpenLibrary,
    ),
    c(
        "Assembly",
        "Contig assembler…",
        "",
        "reads assemble consensus",
        Action::OpenAssembler,
    ),
    c(
        "Assembly",
        "Gibson assembly…",
        "",
        "cloning overlap",
        Action::GibsonAssembler,
    ),
    c(
        "Assembly",
        "Add the open sequences to the assembly",
        "",
        "reads",
        Action::AddOpenToAssembly,
    ),
    c(
        "Assembly",
        "Assemble",
        "",
        "contig run",
        Action::AssembleNow,
    ),
    c(
        "Assembly",
        "Load demo reads",
        "",
        "example",
        Action::DemoReads,
    ),
    c(
        "Tools",
        "Calculator…",
        "Cmd+1",
        "arithmetic maths",
        Action::Calculator,
    ),
    c(
        "Tools",
        "Molecular reaction calculator…",
        "Cmd+2",
        "molar ratio master mix molarity dilution",
        Action::MolecularCalculator,
    ),
    c(
        "Tools",
        "Tm and GC calculator…",
        "",
        "melting temperature primer",
        Action::TmCalculator,
    ),
    c(
        "Tools",
        "Phylogenetic tree…",
        "",
        "neighbour joining",
        Action::Tree,
    ),
    c(
        "Tools",
        "PhyloCNN epidemic estimates…",
        "",
        "outbreak r0 birth death",
        Action::PhyloCnn,
    ),
    c(
        "Tools",
        "Tree tools (Newick)…",
        "",
        "root midpoint distance",
        Action::TreeTools,
    ),
    c(
        "View",
        "Show or hide the map",
        "",
        "plasmid",
        Action::ToggleMap,
    ),
    c(
        "View",
        "Show or hide the feature table",
        "",
        "",
        Action::ToggleFeatures,
    ),
    c(
        "View",
        "Show or hide the complement strand",
        "",
        "",
        Action::ToggleSecondStrand,
    ),
    c(
        "View",
        "Translation in 3 frames",
        "",
        "protein",
        Action::Translation(TransMode::Forward),
    ),
    c(
        "View",
        "Translation in 6 frames",
        "",
        "protein",
        Action::Translation(TransMode::Six),
    ),
    c(
        "View",
        "Theme: same as the system",
        "",
        "appearance dark light",
        Action::Theme(ThemeMode::System),
    ),
    c(
        "View",
        "Theme: light",
        "",
        "appearance",
        Action::Theme(ThemeMode::Light),
    ),
    c(
        "View",
        "Theme: dark",
        "",
        "appearance night",
        Action::Theme(ThemeMode::Dark),
    ),
    c(
        "Help",
        "Command palette",
        "Cmd+K",
        "search commands",
        Action::Palette,
    ),
    c(
        "Help",
        "Keyboard shortcuts",
        "",
        "keys hotkeys cheat sheet",
        Action::Shortcuts,
    ),
    c(
        "Help",
        "About HelixLoom",
        "",
        "version credits",
        Action::About,
    ),
];

/// All the commands, in the order the palette lists them when nothing has been typed.
pub fn all() -> &'static [Command] {
    COMMANDS
}

/// Whether `action` reads or changes the sequence that is open, so that the person has to be looking at it: running it from
/// another workspace (the palette and the menus are on every page) takes them to the Sequences workspace first.
pub fn acts_on_the_sequence(action: &Action) -> bool {
    !matches!(
        action,
        Action::Go(_)
            | Action::Example(_)
            | Action::Open
            | Action::CloseTab(_)
            | Action::Palette
            | Action::Shortcuts
            | Action::Theme(_)
            | Action::About
            | Action::Calculator
            | Action::TmCalculator
            | Action::MolecularCalculator
            | Action::Tree
            | Action::PhyloCnn
            | Action::TreeTools
            | Action::Workflow
            | Action::Reads
    )
}

/// How well `command` answers `query` (every word of it must be found, as a start of a word of the title scoring best):
/// none when a word is missing.
fn score(command: &Command, query: &str) -> Option<i32> {
    let title = command.title.to_lowercase();
    let other = format!("{} {}", command.group, command.words).to_lowercase();
    let mut total = 0;
    for word in query.split_whitespace() {
        let word = word.to_lowercase();
        total += if title.starts_with(&word) {
            100
        } else if title
            .split(|ch: char| !ch.is_alphanumeric())
            .any(|w| w.starts_with(&word))
        {
            70
        } else if title.contains(&word) {
            45
        } else if other.split_whitespace().any(|w| w.starts_with(&word)) {
            25
        } else if other.contains(&word) {
            10
        } else {
            return None;
        };
    }
    // a shorter name is nearer to what was typed
    Some(total * 4 - (command.title.len() / 8) as i32)
}

/// The places in `commands` of those that answer `query`, best first, at most `limit`. With nothing typed it is the first
/// `limit` of the list.
pub fn search(commands: &[Command], query: &str, limit: usize) -> Vec<usize> {
    if query.trim().is_empty() {
        return (0..commands.len().min(limit)).collect();
    }
    let mut found: Vec<(i32, usize)> = commands
        .iter()
        .enumerate()
        .filter_map(|(i, command)| score(command, query).map(|s| (s, i)))
        .collect();
    // the better score first, and the earlier in the list when they tie
    found.sort_by(|a, b| b.0.cmp(&a.0).then(a.1.cmp(&b.1)));
    found.into_iter().take(limit).map(|(_, i)| i).collect()
}

#[cfg(test)]
mod tests {
    use super::*;

    fn titles(commands: &[Command], query: &str) -> Vec<&'static str> {
        search(commands, query, 8)
            .into_iter()
            .map(|i| commands[i].title)
            .collect()
    }

    #[test]
    fn a_name_is_found_by_the_start_of_a_word_and_the_best_match_comes_first() {
        let commands = all();
        assert_eq!(titles(commands, "align")[0], "Alignment");
        let tree = titles(commands, "tree");
        assert_eq!(tree[0], "Trees");
        assert!(tree.contains(&"Phylogenetic tree…") && tree.contains(&"Tree tools (Newick)…"));
        assert_eq!(titles(commands, "mol")[0], "Molecular reaction calculator…");
    }

    #[test]
    fn the_commands_of_the_edit_menu_are_found_by_what_people_type() {
        let commands = all();
        for (typed, title) in [
            ("copy", "Copy"),
            ("paste", "Paste"),
            ("cut", "Cut"),
            ("revcom", "Reverse-complement"),
            ("linear", "Linearize at the caret"),
            ("previous orf", "Previous ORF"),
            ("clear features", "Clear all features"),
            ("assemble", "Assemble"),
        ] {
            assert!(titles(commands, typed).contains(&title), "{typed}");
        }
    }

    #[test]
    fn other_words_find_a_command_too_and_every_word_must_match() {
        let commands = all();
        assert!(titles(commands, "sanger").contains(&"Signals"));
        assert!(titles(commands, "melting").contains(&"Tm and GC calculator…"));
        assert!(titles(commands, "gpu server").contains(&"Workstation"));
        assert!(titles(commands, "gpu zebra").is_empty());
        assert!(titles(commands, "zzzz").is_empty());
    }

    #[test]
    fn with_nothing_typed_the_first_commands_are_listed_in_order() {
        let commands = all();
        let first = search(commands, "  ", 5);
        assert_eq!(first, vec![0, 1, 2, 3, 4]);
        assert_eq!(commands[0].title, "Home");
    }

    #[test]
    fn the_keys_are_written_the_way_the_menus_write_them_and_no_two_commands_share_a_key_or_a_name()
    {
        let commands = all();
        let (mut keys_seen, mut names_seen) = (
            std::collections::HashSet::new(),
            std::collections::HashSet::new(),
        );
        for command in commands {
            assert!(
                names_seen.insert(command.title),
                "{} is listed twice",
                command.title
            );
            if command.keys.is_empty() {
                continue;
            }
            assert!(command.keys.starts_with("Cmd+"), "{}", command.keys);
            assert!(
                keys_seen.insert(command.keys),
                "{} is used twice",
                command.keys
            );
        }
        // the workspaces keep the keys of the rail
        let go: Vec<&str> = commands
            .iter()
            .filter(|c| c.group == "Go to")
            .map(|c| c.keys)
            .collect();
        assert_eq!(go.len(), Workspace::ALL.len());
        assert!(go.iter().all(|k| k.starts_with("Cmd+Alt+")));
    }

    #[test]
    fn what_works_on_the_sequence_takes_the_person_to_it_and_what_does_not_leaves_them_where_they_are(
    ) {
        for action in [
            Action::Undo,
            Action::Find,
            Action::RevCom,
            Action::Upper,
            Action::Save,
            Action::New,
            Action::OrfNext,
        ] {
            assert!(acts_on_the_sequence(&action), "{action:?}");
        }
        for action in [
            Action::Go(Workspace::Trees),
            Action::Example(Workspace::Alignment),
            Action::Open,
            Action::Palette,
            Action::Theme(ThemeMode::Dark),
            Action::Calculator,
            Action::PhyloCnn,
        ] {
            assert!(!acts_on_the_sequence(&action), "{action:?}");
        }
    }
}

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
fn the_keys_are_written_the_way_the_menus_write_them_and_no_two_commands_share_a_key_or_a_name() {
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
fn what_works_on_the_sequence_takes_the_person_to_it_and_what_does_not_leaves_them_where_they_are()
{
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
        Action::SpeciesTree,
    ] {
        assert!(!acts_on_the_sequence(&action), "{action:?}");
    }
}

use super::*;

#[test]
fn every_name_in_the_list_is_a_route_and_each_has_a_word_about_it() {
    for (name, about) in NAMES {
        assert!(named(name).is_some(), "{name}");
        assert!(about.len() > 20, "{name}");
    }
    assert_eq!(NAMES.len(), 5);
}

#[test]
fn the_long_read_route_names_the_files_the_documentation_names() {
    let text = nanopore_plasmid().to_text();
    assert!(text.contains("rotate genes=start-genes.fasta"), "{text}");
    assert!(text.contains("medaka reads=reads"), "{text}");
    assert!(plasmid().to_text().contains("hybrid long=long-reads"));
}

#[test]
fn the_routes_say_in_their_names_what_they_follow() {
    assert!(plasmid()
        .name
        .contains("Plasmids from short and long reads"));
    assert!(nanopore_plasmid().name.contains("Flye") && nanopore_plasmid().name.contains("Medaka"));
    assert!(minion_reads().name.contains("2 kb, 500 Mbp"));
}

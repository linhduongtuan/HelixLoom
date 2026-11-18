use super::*;
use crate::seq::doc::Doc;
use crate::seq::features::Feature;
use crate::seq::inspect::Subject;
use crate::ui::app::Startup;
use crate::ui::platform::Workspace;
use crate::ui::testing::Screen;

/// An app with one plasmid-ish sequence: a gene, and a promoter that overlaps it.
fn app_with_gene() -> App {
    let mut app = App::new(Startup::default());
    let mut doc = Doc::new("p", b"CCATGAAATAAGGCCCCCCCCCC".to_vec());
    doc.features.push(Feature::new("gene", "CDS", 2, 11, false));
    doc.features
        .push(Feature::new("prom", "promoter", 8, 16, false));
    app.docs.clear();
    app.add_doc(doc);
    app.workspace = Workspace::Sequences;
    app
}

#[test]
fn it_says_nothing_until_something_is_selected_and_then_follows_the_selection() {
    let mut app = app_with_gene();
    assert!(app.inspection(0).is_none());

    app.docs[0].select_feature(0);
    let i = app.inspection(0).unwrap();
    assert_eq!(i.subject, Subject::Feature(0));
    assert_eq!(
        i.overlaps,
        vec![1],
        "the promoter shares bases with the gene"
    );

    // dragging off the feature leaves it, and the panel reads the bases that are selected instead
    app.docs[0].select(0, 5);
    assert_eq!(app.inspection(0).unwrap().subject, Subject::Range(0, 5));

    app.docs[0].place_caret(3);
    assert!(app.inspection(0).is_none());
}

#[test]
fn it_reads_with_the_genetic_code_chosen_in_the_preferences() {
    let mut app = app_with_gene();
    app.docs[0].doc.seq[2..11].copy_from_slice(b"ATGTGATAA");
    app.docs[0].select_feature(0);
    assert_eq!(app.inspection(0).unwrap().protein.unwrap().shown, "M**");
    // vertebrate mitochondrial: TGA is Trp
    app.prefs.code = 1;
    assert_eq!(app.inspection(0).unwrap().protein.unwrap().shown, "MW*");
}

#[test]
fn the_panel_says_what_is_selected_in_every_state() {
    let mut app = app_with_gene();
    let mut screen = Screen::new();
    let mut text = |app: &mut App| {
        screen.settle(app);
        screen.text()
    };
    let t = text(&mut app);
    assert!(
        t.contains("INSPECTOR") && t.contains("Nothing selected."),
        "{t}"
    );

    app.docs[0].select_feature(0); // a feature, with its protein and what overlaps it
    let t = text(&mut app);
    assert!(t.contains("9 bp · 3 codons"), "{t}");
    assert!(t.contains("MK*"), "{t}");
    assert!(t.contains("prom (promoter)"), "{t}");

    app.docs[0].select(0, 5); // bases without a feature
    let t = text(&mut app);
    assert!(t.contains("Selection") && t.contains("5 bp"), "{t}");

    app.map_linear = true; // the linear map beside it
    assert!(text(&mut app).contains("Selection"));
    app.prefs.show_map = false; // the inspector alone
    assert!(text(&mut app).contains("Selection"));
    app.prefs.show_inspector = false; // neither
    let t = text(&mut app);
    assert!(!t.contains("INSPECTOR") && !t.contains("Selection"), "{t}");
}

#[test]
fn the_toggle_action_flips_the_setting_back_and_forth() {
    let ctx = egui::Context::default();
    let mut app = app_with_gene();
    assert!(app.prefs.show_inspector);
    app.run(crate::ui::app::Action::ToggleInspector, &ctx);
    assert!(!app.prefs.show_inspector);
    app.run(crate::ui::app::Action::ToggleInspector, &ctx);
    assert!(app.prefs.show_inspector);
}

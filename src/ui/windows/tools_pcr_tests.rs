use super::tools::PcrResult;
use crate::seq::io::to_genbank;
use crate::ui::app::{App, Startup};

/// The demo plasmid with the PCR window's example primers: the forward one starts at base 151 and
/// the reverse one ends at base 506.
fn app_with_a_product() -> App {
    App::new(Startup {
        view: "pcr".to_string(),
        ..Startup::default()
    })
}

#[test]
fn the_example_primers_make_one_product_between_them() {
    let app = app_with_a_product();
    assert_eq!(app.pcr_win.products.len(), 1);
    let p = &app.pcr_win.products[0].product;
    assert_eq!(p.length(), 506 - 150);
    assert_eq!((p.fwd, p.rev), ((150, 180), (480, 506)));
    assert!(p.warnings.is_empty());
}

#[test]
fn opening_a_product_gives_the_document_the_original_makes() {
    let mut app = app_with_a_product();
    let found = &app.pcr_win.products[0];
    let result = PcrResult {
        product: found.product.clone(),
        carried: found.carried.clone(),
    };
    app.open_product(&result);
    let doc = &app.docs.last().expect("the product's tab").doc;
    assert_eq!(doc.name, "PCR pDemo Forward Reverse");
    assert!(!doc.circular && !doc.methylated);
    assert_eq!(doc.seq, result.product.seq);
    assert!(doc.comment.starts_with("PCR:\ntemplate: pDemo\nForward "));
    assert!(doc.comment.ends_with("Product length: 356"));
    // the template's features inside the product (the Reporter ORF is cut off at its end) and the two primers
    let names: Vec<&str> = doc.features.iter().map(|f| f.name.as_str()).collect();
    for expected in [
        "T7 promoter",
        "lac operator",
        "MCS",
        "FLAG",
        "Reporter ORF (demo)",
    ] {
        assert!(
            names.contains(&expected),
            "{expected} is missing from {names:?}"
        );
    }
    let primers: Vec<_> = doc
        .features
        .iter()
        .filter(|f| f.kind == "primer_bind")
        .collect();
    assert_eq!(primers.len(), 2);
    assert_eq!(primers[0].spans, vec![(0, 30)]);
    assert_eq!(primers[1].spans, vec![(330, 356)]);
    let promoter = doc
        .features
        .iter()
        .find(|f| f.name == "T7 promoter")
        .expect("the T7 promoter");
    // bases 151 to 170 of the template, and the product starts at base 151
    assert_eq!(promoter.spans, vec![(0, 20)]);
    // saved as GenBank, the comment travels with the file
    let text = to_genbank(doc);
    assert!(text.contains("COMMENT     PCR:\nCOMMENT     template: pDemo\n"));
}

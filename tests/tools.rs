//! `helixloom tools` as a person runs it: the real program, started with an argument or none.

use std::process::{Command, Output};

fn stderr(out: &Output) -> String {
    String::from_utf8_lossy(&out.stderr).into_owned()
}

fn stdout(out: &Output) -> String {
    String::from_utf8_lossy(&out.stdout).into_owned()
}

#[test]
fn the_tools_command_lists_every_program_with_its_licence_and_how_to_install_it() {
    let out = Command::new(env!("CARGO_BIN_EXE_helixloom"))
        .arg("tools")
        .output()
        .expect("the program starts");
    assert!(out.status.success(), "{}", stderr(&out));
    let text = stdout(&out);
    for label in [
        "MAFFT",
        "IQ-TREE",
        "HMMER",
        "primer3",
        "fastp",
        "curl",
        "Snakemake",
        "Nextflow",
    ] {
        assert!(text.contains(&format!("{label} — ")), "{label} in {text}");
    }
    assert!(
        text.starts_with("HelixLoom looks for these programs"),
        "{text}"
    );
    assert!(
        text.contains("licence: BSD-3-Clause") && text.contains("install: brew install"),
        "{text}"
    );
    assert!(text.contains(" of 8 found."), "{text}");
    // an argument it does not take is said, and nothing is listed
    let bad = Command::new(env!("CARGO_BIN_EXE_helixloom"))
        .args(["tools", "--all"])
        .output()
        .expect("the program starts");
    assert!(!bad.status.success());
    assert!(stderr(&bad).contains("takes no argument") && stdout(&bad).is_empty());
}

#[test]
fn tools_install_refuses_what_is_no_program_of_the_registry_and_says_what_is_installed_already() {
    let run = |args: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_helixloom"))
            .arg("tools")
            .args(args)
            .output()
            .expect("the program starts")
    };
    let out = run(&["install", "blast"]);
    assert_eq!(out.status.code(), Some(1));
    assert!(
        stderr(&out).starts_with(
            "blast is not a program HelixLoom calls; it calls: MAFFT, IQ-TREE, HMMER,"
        ),
        "{}",
        stderr(&out)
    );
    assert_eq!(stdout(&out), "");

    let out = run(&["install"]);
    assert_eq!(out.status.code(), Some(1));
    assert!(
        stderr(&out).starts_with("Which program?"),
        "{}",
        stderr(&out)
    );

    let out = run(&["install", "mafft", "--with", "apt"]);
    assert_eq!(out.status.code(), Some(1));
    assert!(
        stderr(&out).contains("apt is no package manager HelixLoom knows"),
        "{}",
        stderr(&out)
    );

    // an argument that is neither `install` nor nothing is still refused
    let out = run(&["everything"]);
    assert_eq!(out.status.code(), Some(1));
    assert!(
        stderr(&out).contains("`tools` takes no argument but `install NAME`"),
        "{}",
        stderr(&out)
    );

    // the report says how one is installed
    let out = run(&[]);
    assert!(
        stdout(&out)
            .contains("it installs one only when asked to (`helixloom tools install NAME`)"),
        "{}",
        stdout(&out)
    );
}

#[test]
fn the_icon_command_writes_the_mark_as_a_png_of_the_size_asked_and_refuses_a_size_no_icon_has() {
    let dir = std::env::temp_dir().join(format!("helixloom-icon-test-{}", std::process::id()));
    std::fs::create_dir_all(&dir).expect("a folder");
    let file = dir.join("icon.png");
    let run = |args: &[&str]| {
        Command::new(env!("CARGO_BIN_EXE_helixloom"))
            .arg("icon")
            .args(args)
            .output()
            .expect("the program starts")
    };
    let out = run(&[file.to_str().expect("a path"), "64"]);
    assert!(out.status.success(), "{}", stderr(&out));
    let bytes = std::fs::read(&file).expect("the icon");
    assert_eq!(&bytes[..8], b"\x89PNG\r\n\x1a\n");
    // 64 by 64 pixels, with a colour and how see-through it is for each
    assert_eq!(&bytes[16..24], [0, 0, 0, 64, 0, 0, 0, 64]);
    assert_eq!(bytes[25], 6, "RGBA");

    for wrong in ["5", "4000", "big"] {
        let out = run(&[file.to_str().expect("a path"), wrong]);
        assert_eq!(out.status.code(), Some(1), "{wrong}");
        assert!(stderr(&out).contains("from 16 to 2048"), "{}", stderr(&out));
    }
    let out = run(&[]);
    assert_eq!(out.status.code(), Some(1));
    assert!(
        stderr(&out).starts_with("Usage: helixloom icon FILE.png [SIZE]"),
        "{}",
        stderr(&out)
    );
    let _ = std::fs::remove_dir_all(dir);
}

use super::json::{self, Json};
use super::*;
use std::collections::HashMap;

/// A table of what addresses answer, in place of curl.
struct Table(HashMap<String, Reply>);

impl Table {
    fn with(entries: &[(&str, u16, &str)]) -> Table {
        Table(
            entries
                .iter()
                .map(|(u, s, b)| {
                    (
                        u.to_string(),
                        Reply {
                            status: *s,
                            body: b.to_string(),
                        },
                    )
                })
                .collect(),
        )
    }
}

impl Fetcher for Table {
    fn get(&self, url: &str) -> Result<Reply, String> {
        self.0
            .get(url)
            .cloned()
            .ok_or_else(|| format!("no answer for {url}"))
    }
}

const BASE: &str = "https://github.com/linhduongtuan/HelixLoom/releases/download/v0.3.0/";

fn release_json(tag: &str) -> String {
    format!(
        r#"{{"tag_name":"{tag}","name":"HelixLoom {tag}","draft":false,"prerelease":false,
        "html_url":"https://github.com/linhduongtuan/HelixLoom/releases/tag/{tag}",
        "body":"What is new\r\n- k-mers \u00e9 \ud83d\ude00",
        "assets":[
          {{"name":"HelixLoom-{tag}-macos-arm64.dmg","size":9683082,"browser_download_url":"{BASE}HelixLoom-{tag}-macos-arm64.dmg"}},
          {{"name":"HelixLoom-{tag}-windows-x64-NOT-TRIED.zip","size":9689590,"browser_download_url":"{BASE}HelixLoom-{tag}-windows-x64-NOT-TRIED.zip"}},
          {{"name":"SHA256SUMS","size":208,"browser_download_url":"{BASE}SHA256SUMS"}},
          {{"name":"no-url"}}
        ] }}"#
    )
}

#[test]
fn sha256_gives_the_published_test_values() {
    // FIPS 180-4 / NIST examples
    assert_eq!(
        sha256::hex(b""),
        "e3b0c44298fc1c149afbf4c8996fb92427ae41e4649b934ca495991b7852b855"
    );
    assert_eq!(
        sha256::hex(b"abc"),
        "ba7816bf8f01cfea414140de5dae2223b00361a396177a9cb410ff61f20015ad"
    );
    assert_eq!(
        sha256::hex(b"abcdbcdecdefdefgefghfghighijhijkijkljklmklmnlmnomnopnopq"),
        "248d6a61d20638b8e5c026930c3e6039a33ce45964ff2167f6ecedd419db06c1"
    );
    // a message of one million letters a, and lengths either side of the block edges
    assert_eq!(
        sha256::hex(&vec![b'a'; 1_000_000]),
        "cdc76e5c9914fb9281a1c7e284d73e67f1809a48a497200e046d39ccc7112cd0"
    );
    let distinct: std::collections::HashSet<String> =
        (54..70).map(|n| sha256::hex(&vec![b'x'; n])).collect();
    assert_eq!(distinct.len(), 16);
}

#[test]
fn json_is_read_with_its_escapes_numbers_and_nesting_and_refused_when_it_is_not_json() {
    let j = json::parse(
        r#" {"a":[1,-2.5e2,true,null,"x\n\"\u00e9\ud83d\ude00"],"b":{"c":{}},"d":[]} "#,
    )
    .expect("json");
    let a = j.get("a").and_then(Json::as_array).expect("an array");
    assert_eq!(a[0].as_f64(), Some(1.0));
    assert_eq!(a[1].as_f64(), Some(-250.0));
    assert_eq!(a[2].as_bool(), Some(true));
    assert_eq!(a[3], Json::Null);
    assert_eq!(a[4].as_str(), Some("x\n\"é😀"));
    assert_eq!(
        j.get("b").and_then(|b| b.get("c")),
        Some(&Json::Object(vec![]))
    );
    assert_eq!(j.get("missing"), None);
    for bad in [
        "",
        "{",
        r#"{"a":}"#,
        r#"{"a" 1}"#,
        "[1,]",
        "[1 2]",
        r#""abc"#,
        r#""\x""#,
        r#""\ud83d""#,
        "tru",
        "1 2",
        "{} x",
        "\"\u{1}\"",
        "-",
        r#"{1:2}"#,
    ] {
        assert!(json::parse(bad).is_err(), "{bad:?}");
    }
    // nested without end is refused, not followed to the end of the stack
    assert!(json::parse(&"[".repeat(500)).is_err());
}

#[test]
fn a_release_is_read_from_the_answer_of_github() {
    let r = parse_release(&release_json("v0.3.0")).expect("a release");
    assert_eq!(r.tag, "v0.3.0");
    assert_eq!(
        r.page,
        "https://github.com/linhduongtuan/HelixLoom/releases/tag/v0.3.0"
    );
    assert!(
        r.notes.starts_with("What is new\r\n- k-mers é 😀"),
        "{:?}",
        r.notes
    );
    // the asset with no address is passed over
    let names: Vec<&str> = r.assets.iter().map(|a| a.name.as_str()).collect();
    assert_eq!(
        names,
        [
            "HelixLoom-v0.3.0-macos-arm64.dmg",
            "HelixLoom-v0.3.0-windows-x64-NOT-TRIED.zip",
            "SHA256SUMS"
        ]
    );
    assert_eq!(r.assets[0].size, 9_683_082);
    // a draft or a prerelease is no release to update to, and a text with no tag is no release
    for refused in [
        r#"{"tag_name":"v1","draft":true}"#,
        r#"{"tag_name":"v1","prerelease":true}"#,
        "{}",
        "[]",
        "not json",
    ] {
        assert!(parse_release(refused).is_err(), "{refused}");
    }
    // very long notes are cut
    let long = format!(r#"{{"tag_name":"v1.0.0","body":"{}"}}"#, "x".repeat(10_000));
    assert_eq!(parse_release(&long).expect("a release").notes.len(), 4000);
}

#[test]
fn a_version_is_read_from_a_tag_and_compared_by_its_numbers() {
    assert_eq!(version_of("v0.2.1"), Ok((0, 2, 1)));
    assert_eq!(version_of("0.10.0-rc1"), Ok((0, 10, 0)));
    assert_eq!(version_of("V1.4"), Ok((1, 4, 0)));
    for bad in ["", "latest", "v1", "v1.2.3.4", "v.1.2", "vx.y"] {
        assert!(version_of(bad).is_err(), "{bad}");
    }
    assert_eq!(is_newer("v0.2.1", "0.2.0"), Ok(true));
    // by number and not by letters: 10 is after 9
    assert_eq!(is_newer("v0.10.0", "0.9.9"), Ok(true));
    assert_eq!(is_newer("v0.2.0", "0.2.0"), Ok(false));
    assert_eq!(is_newer("v0.1.9", "0.2.0"), Ok(false));
    assert_eq!(is_newer("v1.0.0", "0.99.99"), Ok(true));
    assert!(is_newer("nightly", "0.2.0").is_err());
}

#[test]
fn the_file_for_a_computer_is_chosen_by_its_system_and_processor() {
    let r = parse_release(&release_json("v0.3.0")).expect("a release");
    let name = |os, arch| asset_for(&r, os, arch).map(|a| a.name.as_str());
    assert_eq!(
        name("macos", "aarch64"),
        Some("HelixLoom-v0.3.0-macos-arm64.dmg")
    );
    assert_eq!(
        name("windows", "x86_64"),
        Some("HelixLoom-v0.3.0-windows-x64-NOT-TRIED.zip")
    );
    // no file for these: the person is told to build from the source
    assert_eq!(name("macos", "x86_64"), None);
    assert_eq!(name("linux", "x86_64"), None);
    assert_eq!(name("windows", "aarch64"), None);
}

#[test]
fn the_checksum_of_a_file_is_found_in_the_list_in_either_way_of_writing_it() {
    let a = "a".repeat(64);
    let b = "B".repeat(64);
    let sums = format!(
        "{a}  one.dmg\n{b} *two.zip\nnot a line\n{} three\n",
        "c".repeat(10)
    );
    assert_eq!(checksum_for(&sums, "one.dmg"), Some(a));
    assert_eq!(checksum_for(&sums, "two.zip"), Some("b".repeat(64)));
    // a hash that is not 64 digits, a name that is not there, a name that only ends like one
    assert_eq!(checksum_for(&sums, "three"), None);
    assert_eq!(checksum_for(&sums, "four"), None);
    assert_eq!(checksum_for(&sums, "ne.dmg"), None);
}

#[test]
fn the_check_says_whether_a_later_release_is_out_and_why_it_could_not_tell() {
    let url = latest_url();
    let reply = |status, body: &str| check(&Table::with(&[(&url, status, body)]), "0.2.0");
    match reply(200, &release_json("v0.3.0")).expect("a release") {
        Outcome::Available(r) => assert_eq!(r.tag, "v0.3.0"),
        other => panic!("{other:?}"),
    }
    assert_eq!(
        reply(200, &release_json("v0.2.0")),
        Ok(Outcome::UpToDate {
            latest: "v0.2.0".into()
        })
    );
    assert_eq!(
        reply(200, &release_json("v0.1.0")),
        Ok(Outcome::UpToDate {
            latest: "v0.1.0".into()
        })
    );
    assert!(reply(404, "{}").expect_err("none").contains("No release"));
    assert!(reply(403, "{}")
        .expect_err("limited")
        .contains("try again later"));
    assert!(reply(429, "{}")
        .expect_err("limited")
        .contains("try again later"));
    assert!(reply(500, "{}").expect_err("a failure").contains("500"));
    assert!(reply(200, "<html>").is_err());
    // curl that could not reach it says why
    assert!(check(&Table::with(&[]), "0.2.0")
        .expect_err("no answer")
        .contains("no answer"));
}

/// A folder of its own for a test, removed when the test ends.
struct Folder(PathBuf);

impl Folder {
    fn new(name: &str) -> Folder {
        let path =
            std::env::temp_dir().join(format!("helixloom-update-{}-{name}", std::process::id()));
        std::fs::remove_dir_all(&path).ok();
        Folder(path)
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        std::fs::remove_dir_all(&self.0).ok();
    }
}

const FILE: &[u8] = b"pretend this is a disk image";

#[test]
fn a_download_is_kept_when_its_checksum_matches_and_removed_when_it_does_not() {
    let release = parse_release(&release_json("v0.3.0")).expect("a release");
    let asset = asset_for(&release, "macos", "aarch64")
        .expect("an asset")
        .clone();
    let save = |_: &str, to: &Path| std::fs::write(to, FILE).map_err(|e| e.to_string());
    let sums = |hash: &str| {
        Table::with(&[(
            &format!("{BASE}SHA256SUMS"),
            200,
            &format!("{hash}  {}\n", asset.name),
        )])
    };
    // a match
    let f = Folder::new("match");
    let good =
        download(&release, &asset, &f.0, &save, &sums(&sha256::hex(FILE))).expect("downloaded");
    assert_eq!(
        good,
        Downloaded {
            path: f.0.join(&asset.name),
            verified: true
        }
    );
    assert_eq!(std::fs::read(&good.path).expect("the file"), FILE);
    assert!(!f.0.join(format!("{}.part", asset.name)).exists());
    // no match: nothing is left
    let g = Folder::new("mismatch");
    let said =
        download(&release, &asset, &g.0, &save, &sums(&"0".repeat(64))).expect_err("a mismatch");
    assert!(said.contains("does not match"), "{said}");
    assert!(
        std::fs::read_dir(&g.0)
            .expect("the folder")
            .next()
            .is_none(),
        "the file and its part are gone"
    );
    // no list of checksums, or none for this file: kept, and not said to be checked
    let h = Folder::new("nosums");
    let kept = download(&release, &asset, &h.0, &save, &Table::with(&[])).expect("downloaded");
    assert!(!kept.verified && kept.path.exists());
    let other = Table::with(&[(
        &format!("{BASE}SHA256SUMS"),
        200,
        &format!("{}  other.dmg\n", "a".repeat(64)),
    )]);
    let i = Folder::new("othersums");
    assert!(
        !download(&release, &asset, &i.0, &save, &other)
            .expect("downloaded")
            .verified
    );
    // a download that fails leaves nothing and says why
    let j = Folder::new("failed");
    let fail = |_: &str, _: &Path| Err("the network went away".to_string());
    assert!(download(&release, &asset, &j.0, &fail, &Table::with(&[]))
        .expect_err("failed")
        .contains("network"));
}

#[test]
fn only_a_download_of_a_release_of_this_repository_under_a_plain_name_is_taken() {
    let release = parse_release(&release_json("v0.3.0")).expect("a release");
    let save = |_: &str, to: &Path| std::fs::write(to, FILE).map_err(|e| e.to_string());
    let f = Folder::new("refused");
    let with = |name: &str, url: &str| Asset {
        name: name.into(),
        url: url.into(),
        size: 1,
    };
    let good_url = format!("{BASE}x.dmg");
    for asset in [
        with("../x.dmg", &good_url),
        with("a/b.dmg", &good_url),
        with("a\\b.dmg", &good_url),
        with(".hidden", &good_url),
        with("", &good_url),
        with("x.dmg", "https://example.org/x.dmg"),
        with(
            "x.dmg",
            "http://github.com/linhduongtuan/HelixLoom/releases/download/v1/x.dmg",
        ),
        with(
            "x.dmg",
            "https://github.com/someone-else/HelixLoom/releases/download/v1/x.dmg",
        ),
        with("x.dmg", "file:///etc/passwd"),
    ] {
        assert!(
            download(&release, &asset, &f.0, &save, &Table::with(&[])).is_err(),
            "{asset:?}"
        );
    }
    assert!(
        !f.0.exists()
            || std::fs::read_dir(&f.0)
                .expect("the folder")
                .next()
                .is_none()
    );
}

#[test]
fn what_to_do_with_a_file_depends_on_what_it_is() {
    assert!(next_step(Path::new("/x/HelixLoom-v1-macos-arm64.dmg"))
        .contains("drag HelixLoom to Applications"));
    assert!(next_step(Path::new("C:/x/HelixLoom-v1-windows-x64.zip")).contains("unzip"));
    assert!(next_step(Path::new("/x/other")).contains("install"));
}

#[test]
fn this_programs_own_version_is_a_version() {
    assert!(version_of(env!("CARGO_PKG_VERSION")).is_ok());
}

#[test]
fn the_command_says_what_it_is_for_and_refuses_what_it_does_not_know() {
    assert_eq!(cli::run(&["--help".to_string()]), 0);
    assert_eq!(cli::run(&["--unknown".to_string()]), 1);
    assert_eq!(cli::run(&["--dir".to_string()]), 1, "--dir needs a folder");
}

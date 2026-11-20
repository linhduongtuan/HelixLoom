//! Updating the program: asking GitHub for the latest release of HelixLoom, comparing it with this version, and, when the person says
//! so, downloading the file for this computer and checking it against the release's `SHA256SUMS`.
//!
//! **Nothing is asked of any server unless the person presses the button (or runs `helixloom update`)**: no check at start, none in the
//! background. The one request is a plain GET of `https://api.github.com/repos/linhduongtuan/HelixLoom/releases/latest` made by `curl`,
//! which sends the program's name and nothing of the person's work. The program does not replace itself: a disk image or a zip is
//! saved in the Downloads folder, checked, and shown to the person, who opens it (a program that rewrites its own files cannot be tried on
//! every system, and an ad hoc signed Mac app cannot be replaced behind Gatekeeper's back). The web version updates when its page is
//! loaded again.

use crate::connect::{Fetcher, Reply};
use std::path::{Path, PathBuf};

pub mod cli;
mod json;
mod sha256;
#[cfg(test)]
mod tests;

/// The repository whose releases are looked at.
pub const REPO: &str = "linhduongtuan/HelixLoom";
/// The name of the file of checksums in a release.
pub const SUMS: &str = "SHA256SUMS";
/// The longest part of the release notes that is kept.
const MOST_NOTES: usize = 4000;

pub fn latest_url() -> String {
    format!("https://api.github.com/repos/{REPO}/releases/latest")
}

#[derive(Clone, Debug, PartialEq)]
pub struct Asset {
    pub name: String,
    pub url: String,
    pub size: u64,
}

#[derive(Clone, Debug, PartialEq)]
pub struct Release {
    pub tag: String,
    pub notes: String,
    /// The release's page on GitHub.
    pub page: String,
    pub assets: Vec<Asset>,
}

#[derive(Clone, Debug, PartialEq)]
pub enum Outcome {
    UpToDate { latest: String },
    Available(Release),
}

/// A release from the answer of GitHub (`tag_name`, `body`, `html_url`, `assets`); a draft or a prerelease is not one to update to.
pub fn parse_release(text: &str) -> Result<Release, String> {
    let json = json::parse(text)?;
    let tag = json
        .get("tag_name")
        .and_then(json::Json::as_str)
        .ok_or("The answer of GitHub names no tag.")?
        .to_string();
    if json.get("draft").and_then(json::Json::as_bool) == Some(true)
        || json.get("prerelease").and_then(json::Json::as_bool) == Some(true)
    {
        return Err(format!("The release {tag} is a draft or a prerelease."));
    }
    let notes: String = json
        .get("body")
        .and_then(json::Json::as_str)
        .unwrap_or("")
        .chars()
        .take(MOST_NOTES)
        .collect();
    let page = json
        .get("html_url")
        .and_then(json::Json::as_str)
        .map_or_else(
            || format!("https://github.com/{REPO}/releases/tag/{tag}"),
            str::to_string,
        );
    let assets = json
        .get("assets")
        .and_then(json::Json::as_array)
        .unwrap_or(&[])
        .iter()
        .filter_map(|a| {
            Some(Asset {
                name: a.get("name")?.as_str()?.to_string(),
                url: a.get("browser_download_url")?.as_str()?.to_string(),
                size: a
                    .get("size")
                    .and_then(json::Json::as_f64)
                    .map_or(0, |s| s as u64),
            })
        })
        .collect();
    Ok(Release {
        tag,
        notes,
        page,
        assets,
    })
}

/// The numbers of a version such as `v0.2.1` or `0.10.0-rc1` (what follows the numbers is left off).
pub fn version_of(tag: &str) -> Result<(u64, u64, u64), String> {
    let bare = tag.trim().trim_start_matches(['v', 'V']);
    let numbers: &str = bare
        .split(|c: char| !(c.is_ascii_digit() || c == '.'))
        .next()
        .unwrap_or("");
    let parts: Vec<u64> = numbers
        .split('.')
        .map(str::parse)
        .collect::<Result<_, _>>()
        .map_err(|_| format!("{tag:?} is not a version number."))?;
    match parts.as_slice() {
        [a, b] => Ok((*a, *b, 0)),
        [a, b, c] => Ok((*a, *b, *c)),
        _ => Err(format!("{tag:?} is not a version number (such as v0.2.1).")),
    }
}

/// Whether `tag` is a later version than `current`.
pub fn is_newer(tag: &str, current: &str) -> Result<bool, String> {
    Ok(version_of(tag)? > version_of(current)?)
}

/// The file of a release that is for a computer of this system and processor (`std::env::consts::OS` and `ARCH`), if the release has one:
/// the disk image of a Mac with Apple silicon, the zip of Windows.
pub fn asset_for<'a>(release: &'a Release, os: &str, arch: &str) -> Option<&'a Asset> {
    let wants = |parts: &[&str], ending: &str| {
        release
            .assets
            .iter()
            .find(|a| parts.iter().any(|p| a.name.contains(p)) && a.name.ends_with(ending))
    };
    match (os, arch) {
        ("macos", "aarch64") => wants(&["macos-arm64"], ".dmg"),
        ("macos", "x86_64") => wants(&["macos-x64", "macos-x86_64"], ".dmg"),
        ("windows", "x86_64") => wants(&["windows-x64"], ".zip"),
        _ => None,
    }
}

/// The hash that a `SHA256SUMS` text gives for the file called `name` (`hash  name` or `hash *name` on each line).
pub fn checksum_for(sums: &str, name: &str) -> Option<String> {
    sums.lines().find_map(|line| {
        let (hash, file) = line.trim().split_once(char::is_whitespace)?;
        let file = file.trim_start().trim_start_matches('*');
        (file == name && hash.len() == 64 && hash.bytes().all(|b| b.is_ascii_hexdigit()))
            .then(|| hash.to_ascii_lowercase())
    })
}

/// Ask GitHub for the latest release and say whether it is later than `current`.
pub fn check(fetcher: &dyn Fetcher, current: &str) -> Result<Outcome, String> {
    let Reply { status, body } = fetcher.get(&latest_url())?;
    match status {
        200 => {}
        404 => return Err(format!("No release has been published for {REPO} yet.")),
        403 | 429 => {
            return Err("GitHub refused the request (a limit on the requests from this address): try again later.".into())
        }
        other => return Err(format!("GitHub answered {other}.")),
    }
    let release = parse_release(&body)?;
    if is_newer(&release.tag, current)? {
        Ok(Outcome::Available(release))
    } else {
        Ok(Outcome::UpToDate {
            latest: release.tag,
        })
    }
}

/// A file that was downloaded, and whether its checksum was found and matched.
#[derive(Clone, Debug, PartialEq)]
pub struct Downloaded {
    pub path: PathBuf,
    pub verified: bool,
}

/// The only addresses a file is taken from: the downloads of this repository's releases.
fn is_release_download(url: &str) -> bool {
    url.starts_with(&format!("https://github.com/{REPO}/releases/download/"))
}

/// Save `asset` in `dir` and check it against the release's `SHA256SUMS`. `get_file` saves an address in a file (curl, in the program);
/// `fetcher` gets the text of the checksums. A file that does not match is removed and refused; one with no checksum to be matched is
/// kept and `verified` says it is false. Only a download of this repository's releases is taken, under a name that has no folder in it.
pub fn download(
    release: &Release,
    asset: &Asset,
    dir: &Path,
    get_file: &dyn Fn(&str, &Path) -> Result<(), String>,
    fetcher: &dyn Fetcher,
) -> Result<Downloaded, String> {
    if asset.name.is_empty() || asset.name.contains(['/', '\\']) || asset.name.starts_with('.') {
        return Err(format!(
            "{:?} is not a file name to save under.",
            asset.name
        ));
    }
    if !is_release_download(&asset.url) {
        return Err(format!(
            "{} is not a download of a release of {REPO}.",
            asset.url
        ));
    }
    std::fs::create_dir_all(dir)
        .map_err(|e| format!("Could not make the folder {}: {e}", dir.display()))?;
    let path = dir.join(&asset.name);
    let part = dir.join(format!("{}.part", asset.name));
    get_file(&asset.url, &part)?;
    let bytes =
        std::fs::read(&part).map_err(|e| format!("Could not read what was downloaded: {e}"))?;
    let want = match release.assets.iter().find(|a| a.name == SUMS) {
        Some(sums) if is_release_download(&sums.url) => fetcher
            .get(&sums.url)
            .ok()
            .filter(|r| r.status == 200)
            .and_then(|r| checksum_for(&r.body, &asset.name)),
        _ => None,
    };
    if let Some(want) = &want {
        if *want != sha256::hex(&bytes) {
            let _ = std::fs::remove_file(&part);
            return Err(format!(
                "The file that was downloaded does not match its checksum ({SUMS} of the release): it was removed."
            ));
        }
    }
    std::fs::rename(&part, &path).map_err(|e| format!("Could not save {}: {e}", path.display()))?;
    Ok(Downloaded {
        path,
        verified: want.is_some(),
    })
}

/// The folder that downloads go to: the person's Downloads folder if there is one, else the temporary folder.
pub fn downloads_dir() -> PathBuf {
    let home =
        std::env::var_os(if cfg!(windows) { "USERPROFILE" } else { "HOME" }).map(PathBuf::from);
    home.map(|h| h.join("Downloads"))
        .filter(|d| d.is_dir())
        .unwrap_or_else(std::env::temp_dir)
}

/// Save an address in a file with `curl` (https only, at most five redirects, a few minutes, at most 300 MB).
pub fn curl_file(curl: &Path, url: &str, to: &Path) -> Result<(), String> {
    let output = std::process::Command::new(curl)
        .arg("-q")
        .args([
            "--silent",
            "--show-error",
            "--globoff",
            "--location",
            "--fail",
            "--max-redirs",
            "5",
            "--max-time",
            "600",
        ])
        .args([
            "--max-filesize",
            "300000000",
            "--proto",
            "=https",
            "--proto-redir",
            "=https",
            "--user-agent",
            "HelixLoom",
        ])
        .arg("--output")
        .arg(to)
        .arg("--")
        .arg(url)
        .output()
        .map_err(|e| format!("could not start curl: {e}"))?;
    if output.status.success() {
        Ok(())
    } else {
        let _ = std::fs::remove_file(to);
        let why = String::from_utf8_lossy(&output.stderr);
        Err(if why.trim().is_empty() {
            "curl gave no reason".to_string()
        } else {
            why.trim().to_string()
        })
    }
}

/// Show a file where it is (the Finder, the Explorer, or the folder in the file manager).
pub fn reveal(path: &Path) {
    let folder = path.parent().unwrap_or(path);
    let _ = if cfg!(target_os = "macos") {
        std::process::Command::new("open")
            .arg("-R")
            .arg(path)
            .spawn()
    } else if cfg!(windows) {
        std::process::Command::new("explorer")
            .arg(format!("/select,{}", path.display()))
            .spawn()
    } else {
        std::process::Command::new("xdg-open").arg(folder).spawn()
    };
}

/// What to do with a downloaded file, in a sentence, for this system.
pub fn next_step(file: &Path) -> String {
    let name = file
        .file_name()
        .map_or(String::new(), |n| n.to_string_lossy().into_owned());
    if name.ends_with(".dmg") {
        format!("Quit HelixLoom, open {name}, and drag HelixLoom to Applications (replace the old one). The first start needs a right-click and Open.")
    } else if name.ends_with(".zip") {
        format!("Quit HelixLoom, unzip {name}, and use the new program in place of the old one.")
    } else {
        format!("Quit HelixLoom and install {name}.")
    }
}

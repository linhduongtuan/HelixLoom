//! `helixloom update`: look for a later release, and with `--download` save the file for this computer, checked.

use super::{asset_for, check, curl_file, download, downloads_dir, next_step, Outcome};
use crate::connect::curl::Curl;
use crate::pipeline::tools::locate;
use std::path::PathBuf;

pub const USAGE: &str = "\
Usage: helixloom update [--download] [--dir FOLDER]

  Asks GitHub (one request, made by curl) whether a later release of HelixLoom than this version exists. Nothing is asked unless
  you run this. It says which version is the latest and where its page is.

  --download     also save the file for this computer (the disk image of a Mac with Apple silicon, the zip of Windows), check it against
                 the release's SHA256SUMS, and say what to do with it; the program does not replace itself
  --dir FOLDER   where to save it (the Downloads folder when left out)
  -h, --help     this text
";

pub fn run(args: &[String]) -> i32 {
    if args.iter().any(|a| a == "-h" || a == "--help") {
        print!("{USAGE}");
        return 0;
    }
    match run_inner(args) {
        Ok(()) => 0,
        Err(why) => {
            eprintln!("{why}");
            1
        }
    }
}

fn run_inner(args: &[String]) -> Result<(), String> {
    let (mut want_file, mut dir) = (false, None::<PathBuf>);
    let mut it = args.iter();
    while let Some(a) = it.next() {
        match a.as_str() {
            "--download" => want_file = true,
            "--dir" => dir = Some(PathBuf::from(it.next().ok_or("--dir needs a folder.")?)),
            other => {
                return Err(format!(
                    "{other} is not an option of `helixloom update`; --help lists them."
                ))
            }
        }
    }
    let curl_at =
        locate(&["curl"]).ok_or("curl was not found on this computer: it is what asks GitHub.")?;
    let current = env!("CARGO_PKG_VERSION");
    match check(&Curl::new(curl_at.clone()), current)? {
        Outcome::UpToDate { latest } => {
            println!("HelixLoom {current} is up to date (the latest release is {latest}).");
        }
        Outcome::Available(release) => {
            println!(
                "A later version is out: {} (this is {current}). {}",
                release.tag, release.page
            );
            if want_file {
                let Some(asset) = asset_for(&release, std::env::consts::OS, std::env::consts::ARCH)
                else {
                    return Err("This release has no file for this computer's system: build it from the source (docs/UBUNTU.md).".into());
                };
                let to = dir.unwrap_or_else(downloads_dir);
                let got = download(
                    &release,
                    asset,
                    &to,
                    &|url, file| curl_file(&curl_at, url, file),
                    &Curl::new(curl_at.clone()),
                )?;
                println!(
                    "Saved {}{}",
                    got.path.display(),
                    if got.verified {
                        " (checksum matches)."
                    } else {
                        " (no checksum to compare it with)."
                    }
                );
                println!("{}", next_step(&got.path));
            }
        }
    }
    Ok(())
}

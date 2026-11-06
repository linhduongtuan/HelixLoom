//! What the app remembers about the workstation: its address, the token and whether heavy jobs go there. It is kept in
//! one small file in the person's own configuration folder, which only the person can read (the token is a password).

use std::path::{Path, PathBuf};

const FILE: &str = "workstation.txt";

/// Where the heavy jobs are computed: the person's choice, kept with the connection.
#[derive(Clone, Copy, Debug, Default, PartialEq, Eq)]
pub enum Compute {
    /// On this computer, whatever the job.
    #[default]
    Local,
    /// On the workstation, whatever the job.
    Workstation,
    /// On the workstation when the job is big (`tasks::Weight`) or can use its GPU, on this computer otherwise.
    Automatic,
}

impl Compute {
    pub const ALL: [Compute; 3] = [Compute::Local, Compute::Workstation, Compute::Automatic];

    /// The word the settings file keeps.
    pub fn name(self) -> &'static str {
        match self {
            Compute::Local => "local",
            Compute::Workstation => "workstation",
            Compute::Automatic => "automatic",
        }
    }

    pub fn parse(name: &str) -> Option<Compute> {
        Compute::ALL.into_iter().find(|c| c.name() == name.trim())
    }

    /// What the choice is called on the pages.
    pub fn label(self) -> &'static str {
        match self {
            Compute::Local => "This computer",
            Compute::Workstation => "Workstation",
            Compute::Automatic => "Automatic",
        }
    }
}

#[derive(Clone, Default, PartialEq, Eq)]
pub struct Saved {
    pub address: String,
    pub token: String,
    /// Where the heavy jobs are computed.
    pub compute: Compute,
}

/// Everything but the token, which a log or a failed test must never print.
impl std::fmt::Debug for Saved {
    fn fmt(&self, f: &mut std::fmt::Formatter<'_>) -> std::fmt::Result {
        f.debug_struct("Saved")
            .field("address", &self.address)
            .field("token", &"<hidden>")
            .field("compute", &self.compute)
            .finish()
    }
}

impl Saved {
    pub fn parse(text: &str) -> Saved {
        let m = super::wire::meta(text);
        // a file from before the choice had three ways says `use=1` for "the workstation"
        let compute = m.get("compute").and_then(|c| Compute::parse(c)).unwrap_or(
            if m.get("use") == Some(&"1") {
                Compute::Workstation
            } else {
                Compute::Local
            },
        );
        Saved {
            address: m.get("address").copied().unwrap_or("").trim().to_string(),
            token: m.get("token").copied().unwrap_or("").trim().to_string(),
            compute,
        }
    }

    pub fn to_text(&self) -> String {
        format!(
            "address={}\ntoken={}\ncompute={}\n",
            super::clean_line(&self.address),
            super::clean_line(&self.token),
            self.compute.name()
        )
    }
}

/// The name the program had before it was called HelixLoom: its settings folder is still used while it is the only one.
const FORMER_NAME: &str = "ApECraft";

/// The folder named `name` for a program's settings on `os`, from the environment `var` gives.
fn settings_folder(os: &str, var: &impl Fn(&str) -> Option<String>, name: &str) -> Option<PathBuf> {
    let home = var("HOME").filter(|h| !h.is_empty()).map(PathBuf::from);
    match os {
        "windows" => var("APPDATA").map(|a| PathBuf::from(a).join(name)),
        "macos" => home.map(|h| h.join("Library/Application Support").join(name)),
        _ => var("XDG_CONFIG_HOME")
            .filter(|d| !d.is_empty())
            .map(|d| PathBuf::from(d).join(name.to_lowercase()))
            .or_else(|| home.map(|h| h.join(".config").join(name.to_lowercase()))),
    }
}

/// The folder for a program's settings on `os`, from the environment `var` gives.
pub fn default_dir(os: &str, var: impl Fn(&str) -> Option<String>) -> Option<PathBuf> {
    if let Some(dir) = var("HELIXLOOM_CONFIG_DIR").filter(|d| !d.is_empty()) {
        return Some(PathBuf::from(dir));
    }
    settings_folder(os, &var, "HelixLoom")
}

/// The folder this computer's program keeps its settings in: its own, or the one it had under its former name when that
/// is the only one there is (so that a saved workstation address and token are not lost by the new name).
/// Which of the two folders is in use: the current one, unless it was asked for by name or has no settings file while the
/// former one has. It is the settings file that counts, not the folder: one made only to hold model files hides nothing.
fn choose_dir(current: PathBuf, former: Option<PathBuf>, asked_for: bool) -> PathBuf {
    if asked_for || current.join(FILE).exists() {
        return current;
    }
    former
        .filter(|former| former.join(FILE).exists())
        .unwrap_or(current)
}

pub fn dir() -> Option<PathBuf> {
    // an empty setting is no setting
    let var = |k: &str| std::env::var(k).ok().filter(|v| !v.is_empty());
    let os = std::env::consts::OS;
    let current = default_dir(os, var)?;
    let former = settings_folder(os, &var, FORMER_NAME);
    Some(choose_dir(
        current,
        former,
        var("HELIXLOOM_CONFIG_DIR").is_some(),
    ))
}

pub fn load_from(dir: &Path) -> Saved {
    std::fs::read_to_string(dir.join(FILE))
        .map(|t| Saved::parse(&t))
        .unwrap_or_default()
}

// ponytail: the token is kept as plain text in a file only its owner can read (as an ssh key is); the system's own key store
// would take a crate for each operating system
pub fn save_to(dir: &Path, saved: &Saved) -> Result<(), String> {
    std::fs::create_dir_all(dir).map_err(|e| format!("Could not make {}: {e}", dir.display()))?;
    let path = dir.join(FILE);
    let text = saved.to_text();
    #[cfg(unix)]
    {
        use std::io::Write;
        use std::os::unix::fs::{OpenOptionsExt, PermissionsExt};
        let mut file = std::fs::OpenOptions::new()
            .write(true)
            .create(true)
            .truncate(true)
            .mode(0o600)
            .open(&path)
            .map_err(|e| format!("Could not save {}: {e}", path.display()))?;
        // a file that was there before keeps its old permissions unless they are set again
        let _ = std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o600));
        file.write_all(text.as_bytes())
            .map_err(|e| format!("Could not save {}: {e}", path.display()))
    }
    #[cfg(not(unix))]
    std::fs::write(&path, text).map_err(|e| format!("Could not save {}: {e}", path.display()))
}

#[cfg(test)]
mod tests {
    use super::*;
    use crate::pipeline::tools::work_folder;

    #[test]
    fn settings_read_back_as_written_and_missing_ones_are_empty() {
        for compute in Compute::ALL {
            let s = Saved {
                address: "lab-gpu:7878".into(),
                token: "abc123".into(),
                compute,
            };
            assert_eq!(Saved::parse(&s.to_text()), s, "{compute:?}");
        }
        assert_eq!(Saved::parse(""), Saved::default());
        assert_eq!(Saved::default().compute, Compute::Local);
        assert_eq!(Saved::parse("address=x\nuse=0\n").compute, Compute::Local);
        // a token or address with a line break in it cannot add another setting
        for tricky in ["a\nuse=1", "a\ncompute=workstation"] {
            let s = Saved {
                address: tricky.into(),
                ..Saved::default()
            };
            assert_eq!(
                Saved::parse(&s.to_text()).compute,
                Compute::Local,
                "{tricky}"
            );
        }
    }

    #[test]
    fn a_file_from_before_the_choice_had_three_ways_still_means_the_workstation() {
        // `use=1` was "send the heavy jobs to the workstation"
        let old = Saved::parse("address=lab-gpu\ntoken=abc\nuse=1\n");
        assert_eq!(old.compute, Compute::Workstation);
        // the new word wins over the old one, and a word nobody knows is not taken for a choice
        assert_eq!(
            Saved::parse("compute=automatic\nuse=1\n").compute,
            Compute::Automatic
        );
        assert_eq!(
            Saved::parse("compute=local\nuse=1\n").compute,
            Compute::Local
        );
        assert_eq!(Saved::parse("compute=cloud\n").compute, Compute::Local);
        // the words are the ones the file keeps, whatever the page calls the choices
        for compute in Compute::ALL {
            assert_eq!(Compute::parse(compute.name()), Some(compute));
            assert!(!compute.label().is_empty());
        }
    }

    #[test]
    fn the_settings_folder_follows_each_systems_habit() {
        let env = |pairs: &'static [(&'static str, &'static str)]| {
            move |k: &str| {
                pairs
                    .iter()
                    .find(|(n, _)| *n == k)
                    .map(|(_, v)| v.to_string())
            }
        };
        assert_eq!(
            default_dir("macos", env(&[("HOME", "/Users/a")])),
            Some(PathBuf::from(
                "/Users/a/Library/Application Support/HelixLoom"
            ))
        );
        assert_eq!(
            default_dir("linux", env(&[("HOME", "/home/a")])),
            Some(PathBuf::from("/home/a/.config/helixloom"))
        );
        assert_eq!(
            default_dir(
                "linux",
                env(&[("HOME", "/home/a"), ("XDG_CONFIG_HOME", "/cfg")])
            ),
            Some(PathBuf::from("/cfg/helixloom"))
        );
        assert_eq!(
            default_dir(
                "windows",
                env(&[("APPDATA", "C:\\Users\\a\\AppData\\Roaming")])
            ),
            Some(PathBuf::from("C:\\Users\\a\\AppData\\Roaming").join("HelixLoom"))
        );
        assert_eq!(
            default_dir("linux", env(&[("HELIXLOOM_CONFIG_DIR", "/x")])),
            Some(PathBuf::from("/x"))
        );
        assert_eq!(default_dir("linux", env(&[])), None);
    }

    #[test]
    fn the_former_folder_is_used_only_while_it_alone_has_the_saved_settings() {
        let root = work_folder("settings-fallback").expect("a folder");
        let (current, former) = (root.join("new"), root.join("old"));
        std::fs::create_dir_all(&current).expect("a folder");
        std::fs::create_dir_all(&former).expect("a folder");
        let pick = |asked_for| choose_dir(current.clone(), Some(former.clone()), asked_for);
        assert_eq!(pick(false), current, "neither has settings");
        std::fs::write(former.join(FILE), "address=h:1\n").expect("a file");
        assert_eq!(
            pick(false),
            former,
            "an empty new folder (made for model files) hides nothing"
        );
        assert_eq!(pick(true), current, "a folder asked for by name is the one");
        std::fs::write(current.join(FILE), "address=h:2\n").expect("a file");
        assert_eq!(
            pick(false),
            current,
            "the new folder wins once it has settings"
        );
        assert_eq!(choose_dir(current.clone(), None, false), current);
    }

    #[test]
    fn the_folder_of_the_former_name_is_found_on_each_system() {
        let var = |k: &str| (k == "HOME").then(|| "/Users/a".to_string());
        assert_eq!(
            settings_folder("macos", &var, FORMER_NAME),
            Some(PathBuf::from(
                "/Users/a/Library/Application Support/ApECraft"
            ))
        );
        assert_eq!(
            settings_folder("linux", &var, FORMER_NAME),
            Some(PathBuf::from("/Users/a/.config/apecraft"))
        );
    }

    #[test]
    fn saved_settings_load_again_and_only_the_owner_can_read_the_file() {
        let dir = work_folder("settings").expect("a folder").join("nested");
        assert_eq!(load_from(&dir), Saved::default());
        let s = Saved {
            address: "h:1".into(),
            token: "secret-token-1234".into(),
            compute: Compute::Automatic,
        };
        save_to(&dir, &s).expect("saved");
        assert_eq!(load_from(&dir), s);
        // saved again over a file that anyone could read
        #[cfg(unix)]
        {
            use std::os::unix::fs::PermissionsExt;
            let path = dir.join(FILE);
            std::fs::set_permissions(&path, std::fs::Permissions::from_mode(0o644))
                .expect("loosened");
            save_to(&dir, &s).expect("saved");
            let mode = std::fs::metadata(&path)
                .expect("a file")
                .permissions()
                .mode()
                & 0o777;
            assert_eq!(mode, 0o600);
        }
    }
}

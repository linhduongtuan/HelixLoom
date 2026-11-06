//! The GPUs of this machine, as the server hands them out: found with `nvidia-smi` (NVIDIA, on Linux or Windows) or as
//! the Apple GPU (Metal), or chosen by the person with `--gpus`. A job that needs a GPU gets one of its own by number,
//! and the program it runs sees only that one (`CUDA_VISIBLE_DEVICES`).

use std::process::Command;
use std::sync::mpsc;
use std::time::Duration;

#[derive(Clone, Debug, PartialEq, Eq)]
pub struct Gpu {
    pub id: usize,
    pub name: String,
}

/// How long `nvidia-smi` may take to list the GPUs.
const SMI_SECONDS: u64 = 10;

/// The GPU numbers a person wrote: `0,1,3`, a range `0-3`, or both (`0-1,3`); `none` (or nothing) for no GPU.
pub fn parse_ids(text: &str) -> Result<Vec<usize>, String> {
    let text = text.trim();
    if text.is_empty() || text.eq_ignore_ascii_case("none") {
        return Ok(Vec::new());
    }
    let bad = |part: &str| format!("\"{part}\" is not a GPU number or a range like 0-3.");
    let mut ids = Vec::new();
    for part in text.split(',').map(str::trim) {
        let (from, to) = match part.split_once('-') {
            Some((a, b)) => (
                a.trim().parse().map_err(|_| bad(part))?,
                b.trim().parse().map_err(|_| bad(part))?,
            ),
            None => {
                let n: usize = part.parse().map_err(|_| bad(part))?;
                (n, n)
            }
        };
        if from > to || to > 63 {
            return Err(bad(part));
        }
        ids.extend(from..=to);
    }
    ids.sort_unstable();
    ids.dedup();
    Ok(ids)
}

/// The GPUs in the output of `nvidia-smi -L`: lines like `GPU 0: NVIDIA RTX A6000 (UUID: GPU-...)`.
pub fn parse_nvidia_smi(output: &str) -> Vec<Gpu> {
    output
        .lines()
        .filter_map(|l| {
            let rest = l.strip_prefix("GPU ")?;
            let (id, name) = rest.split_once(':')?;
            let name = name.split(" (UUID").next().unwrap_or(name).trim();
            Some(Gpu {
                id: id.trim().parse().ok()?,
                name: name.to_string(),
            })
        })
        .collect()
}

/// Ask `nvidia-smi` for the GPUs, giving up after a few seconds if it does not answer.
fn nvidia() -> Vec<Gpu> {
    let (send, receive) = mpsc::channel();
    std::thread::spawn(move || {
        let found = Command::new("nvidia-smi")
            .arg("-L")
            .output()
            .ok()
            .filter(|o| o.status.success())
            .map(|o| parse_nvidia_smi(&String::from_utf8_lossy(&o.stdout)))
            .unwrap_or_default();
        let _ = send.send(found);
    });
    receive
        .recv_timeout(Duration::from_secs(SMI_SECONDS))
        .unwrap_or_default()
}

/// Every GPU this machine has that the program can use.
pub fn detect() -> Vec<Gpu> {
    let found = nvidia();
    if !found.is_empty() {
        return found;
    }
    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    if let Some((_, name)) = crate::seq::phylocnn::gpu::device() {
        return vec![Gpu { id: 0, name }];
    }
    Vec::new()
}

/// The GPUs a `--gpus` setting names: `auto` finds them, anything else is a list of numbers.
pub fn select(setting: &str) -> Result<Vec<Gpu>, String> {
    if setting.trim().eq_ignore_ascii_case("auto") {
        return Ok(detect());
    }
    let ids = parse_ids(setting)?;
    if ids.is_empty() {
        return Ok(Vec::new());
    }
    let known = nvidia();
    Ok(ids
        .into_iter()
        .map(|id| Gpu {
            id,
            name: known
                .iter()
                .find(|g| g.id == id)
                .map_or_else(|| format!("GPU {id}"), |g| g.name.clone()),
        })
        .collect())
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn gpu_numbers_are_read_as_lists_and_ranges() {
        assert_eq!(parse_ids("0,1,3"), Ok(vec![0, 1, 3]));
        assert_eq!(parse_ids("0-3"), Ok(vec![0, 1, 2, 3]));
        assert_eq!(parse_ids(" 2 , 0-1 ,2"), Ok(vec![0, 1, 2]));
        assert_eq!(parse_ids("none"), Ok(vec![]));
        assert_eq!(parse_ids(""), Ok(vec![]));
        for bad in ["a", "1,", "3-1", "-1", "0-99", "1.5", "0--2"] {
            assert!(parse_ids(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn the_listing_of_nvidia_smi_gives_numbers_and_names() {
        let out = "GPU 0: NVIDIA RTX A6000 (UUID: GPU-6c1f)\nGPU 1: NVIDIA RTX A6000 (UUID: GPU-77aa)\nnot a gpu line\nGPU x: broken\n";
        assert_eq!(
            parse_nvidia_smi(out),
            vec![
                Gpu {
                    id: 0,
                    name: "NVIDIA RTX A6000".into()
                },
                Gpu {
                    id: 1,
                    name: "NVIDIA RTX A6000".into()
                },
            ]
        );
        assert!(parse_nvidia_smi("").is_empty());
        assert!(parse_nvidia_smi(
            "NVIDIA-SMI has failed because it couldn't communicate with the driver"
        )
        .is_empty());
    }

    #[test]
    fn a_chosen_list_of_gpus_is_used_as_given() {
        let chosen = select("1,3").expect("a list");
        assert_eq!(chosen.iter().map(|g| g.id).collect::<Vec<_>>(), vec![1, 3]);
        assert!(select("none").expect("none").is_empty());
        assert!(select("x").is_err());
    }
}

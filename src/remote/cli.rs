//! The command line of the workstation side: `helixloom serve` (the server), `helixloom token` (a new token),
//! `helixloom status` (ask a server how it is) and `helixloom run-job` (the process each job runs in, started by the
//! server and not by hand).

use super::client::Remote;
use super::gpus;
use super::jobs::Config;
use super::run::{execute, Env};
use super::server;
use super::site::Site;
use super::tasks::{Kind, Request};
use super::wire::{meta, Container};
use super::DEFAULT_PORT;
use crate::pipeline::install;
use crate::pipeline::tools::{self, short_version};
use std::collections::hash_map::RandomState;
use std::hash::{BuildHasher, Hasher};
use std::net::{IpAddr, Ipv4Addr, SocketAddr, ToSocketAddrs};
use std::path::{Path, PathBuf};
use std::sync::Arc;

pub const USAGE: &str = "\
HelixLoom workstation server

  helixloom serve [options]        run the server (jobs sent by the app are computed here)
  helixloom token                  print a new random token
  helixloom status --server HOST:PORT --token-file FILE
                                  ask a server how it is
  helixloom astral -i GENE_TREES   the species tree of a set of gene trees, with no server (`astral --help`)
  helixloom motif -i SEQUENCES     motifs that DNA sequences share, in MEME's text format, with no server (`motif --help`)
  helixloom gff ANNOTATION         read, filter and write GFF3/GTF, and cut transcripts, CDS and proteins out of a genome (`gff --help`)
  helixloom kat COMMAND INPUT...   k-mer analysis in the manner of KAT: hist, gcp, comp, sect (`kat --help`)
  helixloom update [--download]    look for a later release of HelixLoom on GitHub and, if asked, save its file, checked (`update --help`)
  helixloom pipeline run FILE DIR  run a pipeline file on a folder, keeping what was done so a stopped run resumes (`pipeline --help`)
  helixloom checks [WORDS…]        how each analysis was checked: what it was compared with, how close it came, what was not checked
  helixloom tools                  the programs HelixLoom can call (MAFFT, IQ-TREE, HMMER, …): which are found, their licences, how to install them
  helixloom tools install NAME [--with brew|micromamba|mamba|conda] [--yes]
                                  install one of them with the package manager that is on this computer: the command is shown
                                  first and nothing runs without a yes (--yes answers it)
  helixloom gpu-check [--bench]    try the GPUs against the processor, and say from what size they pay (`gpu-check --help`)

Options of `serve`:
  --listen ADDR         where to listen (default 127.0.0.1:7878; use 0.0.0.0:7878 for other computers)
  --token-file FILE     file holding the token that clients must send (or set HELIXLOOM_TOKEN)
  --data-dir DIR        where jobs and results are kept (default ./helixloom-data)
  --gpus auto|none|LIST the GPUs to use, as numbers: 0,1,3 or 0-3 (default auto)
  --cpu-jobs N          jobs that use no GPU at once (default 2)
  --max-upload-mb N     largest request accepted (default 512)
  --keep-days N         days a finished job's result is kept (default 7)
  --job-hours N         hours one job may run (default 24)
  --fold-command CMD    program that folds a protein on a GPU, run as CMD input.fasta output.pdb
  --fold-gpus N         GPUs that program uses for one structure, given to it at once (default 1)
  --name NAME           the name shown to clients (default: this computer's)
  --web-dir DIR         also serve the web version found in DIR (what `trunk build --release` makes, docs/WEB.md): open
                        http://HOST:PORT/ in a browser and use the Workstation page. Read once, at the start
  --allow-origin ORIGIN let the web version hosted at ORIGIN (https://host or http://host:port) call this server from
                        another address; may be given more than once. Not needed for the page --web-dir serves

The server does not encrypt anything: reach it over an SSH tunnel or a private network (see DEPLOY.md).
";

/// The most a person may ask for with the options of `serve`.
const MOST_CPU_JOBS: u64 = 1024;
const MOST_UPLOAD_MB: u64 = 64 * 1024;
const MOST_KEEP_DAYS: u64 = 3650;
const MOST_JOB_HOURS: u64 = 24 * 365;
const MOST_FOLD_GPUS: u64 = 64;

/// A word after a flag, or why there is none.
fn value<'a>(args: &'a [String], flag: &str, i: usize) -> Result<&'a str, String> {
    args.get(i + 1)
        .map(String::as_str)
        .ok_or_else(|| format!("{flag} needs a value."))
}

fn number<T: std::str::FromStr>(text: &str, flag: &str) -> Result<T, String> {
    text.parse()
        .map_err(|_| format!("{flag} needs a whole number, not \"{text}\"."))
}

/// The token in a file: its first line, without spaces round it.
pub fn read_token_file(path: &Path) -> Result<String, String> {
    #[cfg(unix)]
    {
        use std::os::unix::fs::PermissionsExt;
        if let Ok(meta) = std::fs::metadata(path) {
            if meta.permissions().mode() & 0o077 != 0 {
                eprintln!(
                    "Warning: {} can be read by other users; run `chmod 600 {}`.",
                    path.display(),
                    path.display()
                );
            }
        }
    }
    let text = std::fs::read_to_string(path)
        .map_err(|e| format!("Could not read the token file {}: {e}", path.display()))?;
    let token = text.lines().next().unwrap_or("").trim().to_string();
    if token.is_empty() {
        return Err(format!("The token file {} is empty.", path.display()));
    }
    Ok(token)
}

/// A new token: 256 bits read from the system's source of random numbers where there is one (`/dev/urandom`).
// ponytail: elsewhere (Windows) the random keys `std` seeds its hash maps with are used, which are 128 bits at the most
// (the four hashes below share one key); a token from the operating system's generator would need a crate.
pub fn random_token() -> String {
    #[cfg(unix)]
    if let Some(bytes) = urandom(32) {
        return bytes.iter().map(|b| format!("{b:02x}")).collect();
    }
    let nanos = std::time::SystemTime::now()
        .duration_since(std::time::UNIX_EPOCH)
        .map_or(0, |d| d.as_nanos());
    (0..4)
        .map(|i| {
            let mut h = RandomState::new().build_hasher();
            h.write_u128(nanos);
            h.write_u32(std::process::id());
            h.write_usize(i);
            format!("{:016x}", h.finish())
        })
        .collect()
}

#[cfg(unix)]
fn urandom(count: usize) -> Option<Vec<u8>> {
    use std::io::Read;
    let mut bytes = vec![0u8; count];
    std::fs::File::open("/dev/urandom")
        .ok()?
        .read_exact(&mut bytes)
        .ok()?;
    Some(bytes)
}

/// What `helixloom serve` was asked to do.
pub struct ServeOptions {
    pub config: Config,
    pub listen: SocketAddr,
}

/// The options of `serve` as they were typed, before anything is checked.
struct Typed {
    listen: SocketAddr,
    token_file: Option<PathBuf>,
    data_dir: PathBuf,
    gpus: String,
    cpu_jobs: usize,
    upload_mb: usize,
    keep_days: u64,
    job_hours: u64,
    fold: Option<Vec<String>>,
    fold_gpus: usize,
    name: Option<String>,
    web_dir: Option<PathBuf>,
    origins: Vec<String>,
}

impl Default for Typed {
    fn default() -> Typed {
        Typed {
            listen: SocketAddr::new(IpAddr::V4(Ipv4Addr::LOCALHOST), DEFAULT_PORT),
            token_file: None,
            data_dir: PathBuf::from("helixloom-data"),
            gpus: "auto".to_string(),
            cpu_jobs: 2,
            upload_mb: 512,
            keep_days: 7,
            job_hours: 24,
            fold: None,
            fold_gpus: 1,
            name: None,
            web_dir: None,
            origins: Vec::new(),
        }
    }
}

/// The flags after `serve`, each with its value.
fn read_flags(args: &[String]) -> Result<Typed, String> {
    let mut typed = Typed::default();
    let mut i = 0;
    while i < args.len() {
        let flag = args[i].as_str();
        let v = || value(args, flag, i);
        match flag {
            "--listen" => {
                typed.listen = v()?
                    .to_socket_addrs()
                    .map_err(|e| format!("--listen {}: {e}", v().unwrap_or("")))?
                    .next()
                    .ok_or("--listen: no address found.")?;
            }
            "--token-file" => typed.token_file = Some(PathBuf::from(v()?)),
            "--data-dir" => typed.data_dir = PathBuf::from(v()?),
            "--gpus" => typed.gpus = v()?.to_string(),
            "--cpu-jobs" => typed.cpu_jobs = number(v()?, flag)?,
            "--max-upload-mb" => typed.upload_mb = number(v()?, flag)?,
            "--keep-days" => typed.keep_days = number(v()?, flag)?,
            "--job-hours" => typed.job_hours = number(v()?, flag)?,
            "--fold-command" => {
                // ponytail: split at spaces, so a path with a space in it or a quoted word is not understood; put the command in a
                // script and give its path, or add shell-like quoting if that is ever not enough
                let words: Vec<String> = v()?.split_whitespace().map(str::to_string).collect();
                typed.fold = (!words.is_empty()).then_some(words);
            }
            "--fold-gpus" => typed.fold_gpus = number(v()?, flag)?,
            "--name" => typed.name = Some(super::one_line(v()?)),
            "--web-dir" => typed.web_dir = Some(PathBuf::from(v()?)),
            "--allow-origin" => typed.origins.push(
                super::site::parse_origin(v()?).map_err(|why| format!("--allow-origin: {why}"))?,
            ),
            other => return Err(format!("Unknown option {other}.")),
        }
        i += 2;
    }
    Ok(typed)
}

/// The token: from the file named by `--token-file`, else from the environment (`env`).
fn token_for(typed: &Typed, env: impl Fn(&str) -> Option<String>) -> Result<String, String> {
    match &typed.token_file {
        Some(path) => read_token_file(path),
        None => env("HELIXLOOM_TOKEN").map(|t| t.trim().to_string()).filter(|t| !t.is_empty()).ok_or(
            "There is no token. Make one with `helixloom token > token.txt`, then use --token-file token.txt (or set HELIXLOOM_TOKEN)."
                .to_string(),
        ),
    }
}

/// Refuse a number that is zero, or more than anyone needs (and small enough that no sum or shift of them can overflow).
fn check_numbers(t: &Typed) -> Result<(), String> {
    let asked = [
        ("--cpu-jobs", t.cpu_jobs as u64, MOST_CPU_JOBS),
        ("--max-upload-mb", t.upload_mb as u64, MOST_UPLOAD_MB),
        ("--keep-days", t.keep_days, MOST_KEEP_DAYS),
        ("--job-hours", t.job_hours, MOST_JOB_HOURS),
        ("--fold-gpus", t.fold_gpus as u64, MOST_FOLD_GPUS),
    ];
    if asked.iter().any(|&(_, n, _)| n == 0) {
        return Err(
            "--cpu-jobs, --max-upload-mb, --keep-days, --job-hours and --fold-gpus must be at least 1."
                .into(),
        );
    }
    match asked.iter().find(|&&(_, n, most)| n > most) {
        Some((flag, _, most)) => Err(format!("{flag} may be at most {most}.")),
        None => Ok(()),
    }
}

/// `path` as the file system really has it, as far as it exists (a folder that is not made yet is its nearest existing
/// folder, really, and the rest).
fn resolved(path: &Path) -> PathBuf {
    let mut rest = PathBuf::new();
    let mut head = path.to_path_buf();
    loop {
        if let Ok(real) = std::fs::canonicalize(&head) {
            return real.join(rest);
        }
        match (head.parent(), head.file_name()) {
            (Some(parent), Some(name)) => {
                rest = Path::new(name).join(rest);
                head = parent.to_path_buf();
            }
            _ => return path.to_path_buf(),
        }
    }
}

/// The page is handed to anyone, so the folder it is read from must hold the web build and nothing the server keeps: not a file
/// with the token in it, not the token's own file, not the folder of the jobs and their results.
fn check_web_dir(
    dir: &Path,
    site: &Site,
    token: &str,
    data_dir: &Path,
    token_file: Option<&Path>,
) -> Result<(), String> {
    let public = format!(
        "{} is handed to anyone who asks, so it must hold only the web build",
        dir.display()
    );
    if let Some(path) = site.file_containing(token.as_bytes()) {
        return Err(format!("{public}, but {path} has the token in it."));
    }
    let root = resolved(dir);
    if resolved(data_dir).starts_with(&root) {
        return Err(format!(
            "{public}, but the folder of the jobs (--data-dir) is inside it."
        ));
    }
    if token_file.is_some_and(|file| resolved(file).starts_with(&root)) {
        return Err(format!("{public}, but the token's file is inside it."));
    }
    Ok(())
}

/// The programs of the analyses that are on this computer, with their versions, for the status page.
fn found_tools() -> Vec<(String, String)> {
    tools::ALL
        .iter()
        .filter_map(|p| {
            tools::detect(p).map(|f| (p.label.to_string(), short_version(&f.version).to_string()))
        })
        .collect()
}

/// Read the options of `serve`; the token comes from `--token-file` or `HELIXLOOM_TOKEN` (`env`).
pub fn serve_options(
    args: &[String],
    env: impl Fn(&str) -> Option<String>,
) -> Result<ServeOptions, String> {
    let typed = read_flags(args)?;
    let token = token_for(&typed, env)?;
    check_numbers(&typed)?;
    // an absolute path, so that the folders the jobs are given (their temporary files) are right wherever a program runs
    let data_dir = std::path::absolute(&typed.data_dir).unwrap_or_else(|_| typed.data_dir.clone());
    let mut config = Config::new(token, data_dir);
    config.cpu_jobs = typed.cpu_jobs;
    config.max_upload = typed.upload_mb << 20;
    config.keep = std::time::Duration::from_secs(typed.keep_days * 24 * 3600);
    config.job_limit = std::time::Duration::from_secs(typed.job_hours * 3600);
    config.fold_command = typed.fold;
    config.fold_gpus = typed.fold_gpus;
    if let Some(n) = typed.name.filter(|n| !n.is_empty()) {
        config.name = n;
    }
    config.gpus = gpus::select(&typed.gpus)?;
    config.tools = found_tools();
    config.site = Arc::new(match &typed.web_dir {
        Some(dir) => {
            let site = Site::load(dir, &typed.origins)?;
            check_web_dir(
                dir,
                &site,
                &config.token,
                &config.data_dir,
                typed.token_file.as_deref(),
            )?;
            site
        }
        None => Site::with_origins(&typed.origins),
    });
    Ok(ServeOptions {
        config,
        listen: typed.listen,
    })
}

/// What a person who has just started the server needs to see; the first line is what a script waits for.
fn announce(server: &server::Server) {
    let config = &server.shared.config;
    println!("HelixLoom server listening on {}", server.addr);
    println!("Name: {}", config.name);
    if config.gpus.is_empty() {
        println!("GPUs: none (jobs that need a GPU will be refused)");
    }
    for g in &config.gpus {
        println!("GPU {}: {}", g.id, g.name);
    }
    println!("Data folder: {}", config.data_dir.display());
    if config.site.has_page() {
        // a server that listens everywhere has no one address to open: HOST is the computer's own
        let at = if server.addr.ip().is_unspecified() {
            format!(
                "HOST:{} (HOST: this computer's address)",
                server.addr.port()
            )
        } else {
            server.addr.to_string()
        };
        println!(
            "Web version: open http://{at}/ in a browser, then Workstation, and type the token."
        );
    }
    for origin in config.site.origins() {
        println!("Pages of {origin} may call this server.");
    }
    println!(
        "A job may run {} hours; results are kept {} days.",
        config.job_limit.as_secs() / 3600,
        config.keep.as_secs() / 86400
    );
    if !server.addr.ip().is_loopback() {
        println!("Note: this port is open to other computers and nothing sent through it is encrypted: use a private network or an SSH tunnel, and never the open internet.");
    }
    println!(
        "In the app, open Workstation, type this computer's address with port {} and the token, and press Connect.",
        server.addr.port()
    );
    {
        use std::io::Write;
        let _ = std::io::stdout().flush();
    }
}

fn serve(args: &[String]) -> i32 {
    if args.iter().any(|a| a == "--help" || a == "-h") {
        print!("{USAGE}");
        return 0;
    }
    // the token is read once and then taken out of the environment, so that none of the programs this server runs
    // (nvidia-smi, the tools it looks for, the jobs) is handed it
    let token_in_env = std::env::var("HELIXLOOM_TOKEN").ok();
    std::env::remove_var("HELIXLOOM_TOKEN");
    let env = |k: &str| {
        if k == "HELIXLOOM_TOKEN" {
            token_in_env.clone()
        } else {
            std::env::var(k).ok()
        }
    };
    let options = match serve_options(args, env) {
        Ok(o) => o,
        Err(why) => {
            eprintln!("{why}\n\nRun `helixloom serve --help` for the options.");
            return 2;
        }
    };
    let server = match server::start(options.config, options.listen, true) {
        Ok(s) => s,
        Err(why) => {
            eprintln!("{why}");
            return 1;
        }
    };
    announce(&server);
    loop {
        std::thread::park();
    }
}

/// The server and the token file named after `status`, or why they are not.
fn status_flags(args: &[String]) -> Result<(String, PathBuf), String> {
    let (mut server, mut token_file) = (None, None);
    let mut i = 0;
    while i < args.len() {
        match args[i].as_str() {
            "--server" => server = Some(value(args, "--server", i)?.to_string()),
            "--token-file" => token_file = Some(PathBuf::from(value(args, "--token-file", i)?)),
            other => return Err(format!("Unknown option {other}.")),
        }
        i += 2;
    }
    server
        .zip(token_file)
        .ok_or_else(|| "Usage: helixloom status --server HOST:PORT --token-file FILE".to_string())
}

/// What a server said about itself, as lines.
fn print_status(remote: &Remote, s: &super::client::ServerStatus) {
    println!(
        "{} ({}), HelixLoom {}, protocol {}",
        s.name,
        remote.address(),
        s.version,
        s.protocol
    );
    println!(
        "{} running, {} waiting; {} CPU jobs at once; results kept {} hours",
        s.running, s.queued, s.cpu_jobs, s.keep_hours
    );
    if s.gpus.is_empty() {
        println!("GPUs: none");
    }
    for g in &s.gpus {
        let state = g
            .job
            .as_deref()
            .map_or("free".to_string(), |j| format!("busy with {j}"));
        println!("GPU {}: {} ({state})", g.id, g.name);
    }
    for (name, version) in &s.tools {
        println!("{name} {version}");
    }
    println!(
        "Local folding program: {}",
        if s.fold { "yes" } else { "no" }
    );
}

fn status(args: &[String]) -> i32 {
    let (server, token_file) = match status_flags(args) {
        Ok(found) => found,
        Err(why) => {
            eprintln!("{why}");
            return 2;
        }
    };
    let outcome = read_token_file(&token_file)
        .and_then(|t| Remote::new(&server, &t))
        .and_then(|r| r.status().map(|s| (r, s)));
    match outcome {
        Ok((remote, s)) => {
            print_status(&remote, &s);
            0
        }
        Err(why) => {
            eprintln!("{why}");
            1
        }
    }
}

/// How often a job's process looks whether the server that started it is still there.
#[cfg(unix)]
const SERVER_WATCH: std::time::Duration = std::time::Duration::from_secs(2);

/// A job's process ends when the server that started it is gone (it was stopped with Ctrl-C, or it crashed): nobody
/// would collect the result, and a long job would hold its GPU for hours. The whole process group is stopped, so the
/// programs the job started end with it.
// ponytail: Unix only (std has no parent process on Windows); there a job whose server is gone runs on until it ends.
#[cfg(unix)]
fn end_with_the_server(server: Option<u32>) {
    // the server says which process it is; a job started by hand has to go by its parent
    let server = server.unwrap_or_else(std::os::unix::process::parent_id);
    std::thread::spawn(move || loop {
        std::thread::sleep(SERVER_WATCH);
        if std::os::unix::process::parent_id() != server {
            // the server is gone: stop this job's process group, which is this process and what it started
            super::jobs::kill_tree(std::process::id());
            std::process::exit(1);
        }
    });
}

/// The process of the server that started a job, from `--server-pid N` among the words before any `--`.
#[cfg(unix)]
fn server_pid(args: &[String]) -> Option<u32> {
    let own: Vec<&String> = args.iter().take_while(|a| a.as_str() != "--").collect();
    own.windows(2)
        .find(|w| w[0] == "--server-pid")
        .and_then(|w| w[1].parse().ok())
}

/// The process a job runs in: read the request in the job's folder, run it, write the result beside it.
fn run_job(args: &[String]) -> i32 {
    let Some(dir) = args.first().map(PathBuf::from) else {
        eprintln!("Usage: helixloom run-job FOLDER [--fold-seconds N] [--flow-cache DIR] [-- fold command words]");
        return 2;
    };
    let mut env = Env::default();
    let mut rest = args[1..].iter();
    while let Some(a) = rest.next() {
        match a.as_str() {
            "--fold-seconds" => {
                env.fold_seconds = rest.next().and_then(|s| s.parse().ok()).unwrap_or(0)
            }
            "--flow-cache" => env.flow_cache = rest.next().map(PathBuf::from),
            // read by `dispatch` before the job starts
            "--server-pid" => {
                rest.next();
            }
            "--" => {
                let words: Vec<String> = rest.by_ref().cloned().collect();
                env.fold_command = (!words.is_empty()).then_some(words);
            }
            _ => {}
        }
    }
    let outcome = std::panic::catch_unwind(|| run_in(&dir, &env));
    let why = match outcome {
        Ok(Ok(())) => return 0,
        Ok(Err(why)) => why,
        Err(_) => "The job crashed.".to_string(),
    };
    let _ = std::fs::write(dir.join("error.txt"), super::one_line(&why));
    eprintln!("{why}");
    1
}

fn run_in(dir: &Path, env: &Env) -> Result<(), String> {
    let bytes = std::fs::read(dir.join("request.bin"))
        .map_err(|e| format!("Could not read the request: {e}"))?;
    let container = Container::decode(&bytes)?;
    let kind = meta(container.require("job")?)
        .get("kind")
        .and_then(|k| Kind::parse(k))
        .ok_or("The job's kind is not one this program runs.")?;
    let response = execute(Request::decode(kind, &container)?, env)?;
    let result = dir.join("result.bin");
    let partial = dir.join("result.part");
    std::fs::write(&partial, response.encode().encode())
        .and_then(|()| std::fs::rename(&partial, &result))
        .map_err(|e| format!("Could not write the result: {e}"))
}

/// Run the command if `args` (the program's arguments, the first being its name) ask for one of the server's;
/// `helixloom gpu-check`: try the GPUs against the processor (`accel::check`).
fn gpu_check(args: &[String]) -> i32 {
    #[cfg(all(feature = "gpu", not(target_arch = "wasm32")))]
    {
        crate::accel::check::run(args)
    }
    #[cfg(not(all(feature = "gpu", not(target_arch = "wasm32"))))]
    {
        let _ = args;
        eprintln!("This program was built without GPU support (`--no-default-features`): build it with the default features, or with `--features cuda`.");
        2
    }
}

/// `helixloom icon FILE.png [SIZE]`: the program's mark as the picture of the app, for the script that makes the Mac app
/// (`scripts/bundle-macos.sh`).
fn icon_command(args: &[String]) -> i32 {
    let (least, most) = crate::ui::brand::ICON_SIZES;
    let (Some(path), size) = (args.first(), args.get(1)) else {
        eprintln!("Usage: helixloom icon FILE.png [SIZE]   (SIZE in pixels each way, {least} to {most}; 1024 when left out)");
        return 1;
    };
    let size = match size.map(|s| s.parse::<u32>()) {
        None => 1024,
        Some(Ok(size)) if (least..=most).contains(&size) => size,
        Some(_) => {
            eprintln!("The size of an icon is a whole number of pixels from {least} to {most}.");
            return 1;
        }
    };
    let written = crate::ui::headless::png_bytes(&crate::ui::brand::icon(size))
        .and_then(|bytes| std::fs::write(path, bytes).map_err(|e| e.to_string()));
    match written {
        Ok(()) => 0,
        Err(why) => {
            eprintln!("Could not write {path}: {why}");
            1
        }
    }
}

/// `helixloom checks [WORDS…]`: how each analysis was checked, as the page of the Help menu says it; the words leave the checks that
/// hold them all.
fn checks_report(args: &[String]) -> i32 {
    let words = args.join(" ");
    let found = crate::checks::matching(&words);
    if found.is_empty() {
        eprintln!("No analysis matches {words}.");
        return 1;
    }
    print!("{}", crate::checks::to_text(&found));
    0
}

/// What `helixloom tools install` was asked: the program, the package manager if one was named, and whether the yes was given.
struct InstallAsk {
    program: &'static tools::Program,
    with: Option<install::Manager>,
    yes: bool,
}

fn install_ask(args: &[String]) -> Result<InstallAsk, String> {
    let called: Vec<&str> = tools::ALL.iter().map(|p| p.label).collect();
    let (mut name, mut with, mut yes) = (None, None, false);
    let mut words = args.iter();
    while let Some(word) = words.next() {
        match word.as_str() {
            "--yes" | "-y" => yes = true,
            "--with" => {
                let manager = words
                    .next()
                    .ok_or("--with needs a package manager: brew, micromamba, mamba or conda.")?;
                with = Some(install::Manager::parse(manager).ok_or_else(|| {
                    format!("{manager} is no package manager HelixLoom knows: brew, micromamba, mamba or conda.")
                })?);
            }
            option if option.starts_with('-') => {
                return Err(format!(
                    "The option {option} is not known; `tools install` takes --with and --yes."
                ))
            }
            program if name.is_none() => name = Some(program),
            extra => {
                return Err(format!(
                    "{extra} is one argument too many: one program is installed at a time."
                ))
            }
        }
    }
    let name =
        name.ok_or_else(|| format!("Which program? HelixLoom calls: {}.", called.join(", ")))?;
    let program = tools::by_name(name).ok_or_else(|| {
        format!(
            "{name} is not a program HelixLoom calls; it calls: {}.",
            called.join(", ")
        )
    })?;
    Ok(InstallAsk { program, with, yes })
}

/// `helixloom tools install NAME`: install one of the programs HelixLoom calls with a package manager that is on this computer.
/// The command is printed first, and nothing runs without a yes.
fn tools_install(args: &[String]) -> i32 {
    let ask = match install_ask(args) {
        Ok(ask) => ask,
        Err(why) => {
            eprintln!("{why}");
            return 1;
        }
    };
    let program = ask.program;
    if let Some(found) = tools::detect(program) {
        println!(
            "{} is installed already: {}",
            program.label,
            found.path.display()
        );
        return 0;
    }
    let active = std::env::var("CONDA_PREFIX").ok();
    let plan = match install::plan(program, &install::managers(), ask.with, active.as_deref()) {
        Ok(plan) => plan,
        Err(why) => {
            eprintln!("{why}");
            return 1;
        }
    };
    println!("{}\n", plan.notice(program));
    if !ask.yes {
        use std::io::Write;
        print!("Run it? [y/N] ");
        let _ = std::io::stdout().flush();
        let mut answer = String::new();
        let _ = std::io::stdin().read_line(&mut answer);
        if !install::agreed(&answer) {
            println!("Nothing was installed.");
            return 1;
        }
    }
    if let Err(why) = install::run(&plan, |line| println!("{line}")) {
        eprintln!("{why}");
        return 1;
    }
    match tools::detect(program) {
        Some(found) => {
            println!("{} is installed: {}", program.label, found.path.display());
            0
        }
        None => {
            eprintln!(
                "{} ended well, but {} is still not found on the PATH or in Homebrew's folders.",
                plan.line(),
                program.label
            );
            1
        }
    }
}

/// `helixloom tools`: the programs HelixLoom can call, whether each is found, its licence and how to install it. Runs each program
/// that is found once, to ask its version; installs nothing.
fn tools_report(args: &[String]) -> i32 {
    if args.first().map(String::as_str) == Some("install") {
        return tools_install(&args[1..]);
    }
    if let Some(extra) = args.first() {
        eprintln!("`tools` takes no argument but `install NAME`, not {extra}.");
        return 1;
    }
    let found: Vec<_> = tools::ALL.iter().map(|p| (*p, tools::detect(p))).collect();
    print!("{}", tools::report(&found));
    0
}

/// `None` when the app should open its window instead.
pub fn dispatch(args: &[String]) -> Option<i32> {
    let rest = args.get(2..).unwrap_or(&[]);
    Some(match args.get(1)?.as_str() {
        "serve" => serve(rest),
        "token" => {
            println!("{}", random_token());
            0
        }
        "status" => status(rest),
        "astral" | "species-tree" => crate::seq::astral::cli::run(rest),
        "motif" | "motifs" => crate::seq::motif::cli::run(rest),
        "gff" | "gffread" => crate::seq::gff::cli::run(rest),
        "update" => crate::update::cli::run(rest),
        "kat" => crate::seq::kmer::cli::run(rest),
        "pipeline" => crate::pipeline::flow::cli::main(rest),
        "tools" => tools_report(rest),
        "checks" => checks_report(rest),
        "icon" => icon_command(rest),
        "gpu-check" => gpu_check(rest),
        "run-job" => {
            #[cfg(unix)]
            end_with_the_server(server_pid(rest));
            run_job(rest)
        }
        "help" | "--help" | "-h" => {
            print!("{USAGE}");
            0
        }
        _ => return None,
    })
}

#[cfg(test)]
mod tests;

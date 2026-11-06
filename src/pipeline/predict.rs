//! Online protein structure prediction: the sequence goes to the ESM Metagenomic Atlas web service, which runs
//! ESMFold (Lin et al., Science 2023) and answers with a PDB file whose temperature-factor column holds the
//! confidence (pLDDT). The request is made with the `curl` program, so no network code is linked into this one.

use std::io::Write;
use std::path::Path;
use std::process::{Command, Stdio};

/// The service's address for sequences in, PDB text out (api.esmatlas.com/foldSequence/v1/pdb/).
pub const ESM_ATLAS_URL: &str = "https://api.esmatlas.com/foldSequence/v1/pdb/";
/// The longest sequence the service takes.
pub const MAX_LENGTH: usize = 400;
/// Shorter than this is a peptide, not a fold.
pub const MIN_LENGTH: usize = 10;
/// How long to wait for the answer before giving up.
const TIMEOUT_SECONDS: u32 = 180;
const AMINO_ACIDS: &[u8; 20] = b"ACDEFGHIKLMNPQRSTVWY";
/// How much of an unexpected answer is quoted in an error.
const QUOTE: usize = 160;

/// The amino acids in `text`: plain letters or FASTA (only the first record is taken). Spaces, digits and line
/// breaks are ignored; any other letter is an error, as the service reads only the twenty standard ones.
pub fn clean_sequence(text: &str) -> Result<Vec<u8>, String> {
    let trimmed = text.trim_start();
    let body: String = if trimmed.starts_with('>') {
        // FASTA: the lines after the first header, up to the next one
        trimmed
            .lines()
            .skip(1)
            .take_while(|l| !l.starts_with('>'))
            .collect()
    } else {
        trimmed.to_string()
    };
    let letters: Vec<u8> = body
        .bytes()
        .filter(u8::is_ascii_alphabetic)
        .map(|c| c.to_ascii_uppercase())
        .collect();
    if let Some(&bad) = letters.iter().find(|c| !AMINO_ACIDS.contains(c)) {
        return Err(format!(
            "'{}' is not one of the twenty amino acids ({}).",
            char::from(bad),
            String::from_utf8_lossy(AMINO_ACIDS)
        ));
    }
    match letters.len() {
        n if n < MIN_LENGTH => Err(format!("The sequence has {n} residues; at least {MIN_LENGTH} are needed.")),
        n if n > MAX_LENGTH => Err(format!(
            "The sequence has {n} residues; the service takes at most {MAX_LENGTH}. Use a domain or an ORF part."
        )),
        _ => Ok(letters),
    }
}

/// POST `sequence` to `url` with `curl` and return the PDB text of the answer.
pub fn predict(curl: &Path, url: &str, sequence: &[u8]) -> Result<String, String> {
    let mut child = Command::new(curl)
        .args(["--silent", "--show-error", "--max-time"])
        .arg(TIMEOUT_SECONDS.to_string())
        .args([
            "--request",
            "POST",
            "--data-binary",
            "@-",
            "--write-out",
            "\n%{http_code}",
        ])
        .arg(url)
        .stdin(Stdio::piped())
        .stdout(Stdio::piped())
        .stderr(Stdio::piped())
        .spawn()
        .map_err(|e| format!("Could not start curl: {e}"))?;
    if let Some(mut stdin) = child.stdin.take() {
        stdin
            .write_all(sequence)
            .map_err(|e| format!("Could not send the sequence to curl: {e}"))?;
    }
    let output = child
        .wait_with_output()
        .map_err(|e| format!("curl stopped unexpectedly: {e}"))?;
    if !output.status.success() {
        let why = String::from_utf8_lossy(&output.stderr);
        let why = why.trim();
        return Err(format!(
            "The request failed: {}",
            if why.is_empty() {
                "curl gave no reason"
            } else {
                why
            }
        ));
    }
    let answer = String::from_utf8_lossy(&output.stdout);
    let (body, status) = answer.rsplit_once('\n').unwrap_or((&answer, ""));
    if status.trim() != "200" {
        let quote: String = body.trim().chars().take(QUOTE).collect();
        return Err(format!("The service answered {} {quote}", status.trim())
            .trim_end()
            .to_string());
    }
    if !body.lines().any(|l| l.starts_with("ATOM")) {
        let quote: String = body.trim().chars().take(QUOTE).collect();
        return Err(format!("The service did not return a structure: {quote}"));
    }
    Ok(body.to_string())
}

#[cfg(test)]
mod tests;

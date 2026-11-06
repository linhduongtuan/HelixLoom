use super::*;
use crate::pipeline::tools;
use std::io::{Read, Write};
use std::net::TcpListener;

const PEPTIDE: &[u8] = b"MKTAYIAKQRQISFVKSHFSRQLEERLGLIEVQAPILSRVGDGTQDNLSGAEKAVQVKVKALPDAQFEVV";

const FAKE_PDB: &str =
    "ATOM      1  N   MET A   1       1.000   2.000   3.000  1.00 91.50           N\nEND\n";

fn curl() -> Option<std::path::PathBuf> {
    tools::locate(&["curl"])
}

/// A one-shot web server on this computer: it answers the first request with `status` and `body`, and hands back
/// what the request carried.
fn serve(status: &'static str, body: &'static str) -> (String, std::thread::JoinHandle<String>) {
    let listener = TcpListener::bind("127.0.0.1:0").expect("a free port");
    let url = format!("http://{}/fold", listener.local_addr().expect("an address"));
    let handle = std::thread::spawn(move || {
        let (mut stream, _) = listener.accept().expect("a request");
        let mut request = Vec::new();
        let mut buffer = [0u8; 1024];
        loop {
            let n = stream.read(&mut buffer).expect("request bytes");
            request.extend_from_slice(&buffer[..n]);
            let text = String::from_utf8_lossy(&request).into_owned();
            if let Some((head, rest)) = text.split_once("\r\n\r\n") {
                let length: usize = head
                    .lines()
                    .find_map(|l| {
                        l.to_ascii_lowercase()
                            .strip_prefix("content-length:")
                            .map(|v| v.trim().parse().unwrap_or(0))
                    })
                    .unwrap_or(0);
                if rest.len() >= length || n == 0 {
                    break;
                }
            } else if n == 0 {
                break;
            }
        }
        let reply = format!(
            "HTTP/1.1 {status}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{body}",
            body.len()
        );
        stream.write_all(reply.as_bytes()).expect("the reply");
        String::from_utf8_lossy(&request).into_owned()
    });
    (url, handle)
}

#[test]
fn a_sequence_is_cleaned_to_the_twenty_amino_acids() {
    assert_eq!(
        clean_sequence("mktayiakqr qisfvk\n 12 shf").unwrap(),
        b"MKTAYIAKQRQISFVKSHF"
    );
    // FASTA: the first record only
    let fasta = ">first protein\nMKTAYIAKQR\nQISFVKSHFS\n>second\nGGGGGGGGGGGG\n";
    assert_eq!(clean_sequence(fasta).unwrap(), b"MKTAYIAKQRQISFVKSHFS");
}

#[test]
fn sequences_that_cannot_be_predicted_are_refused_in_plain_words() {
    assert!(clean_sequence("MKTAYIAKQRQXISFV")
        .unwrap_err()
        .contains("'X' is not one of the twenty"));
    assert!(clean_sequence("MKT").unwrap_err().contains("at least 10"));
    assert!(clean_sequence("").unwrap_err().contains("0 residues"));
    let long = "A".repeat(MAX_LENGTH + 1);
    assert!(clean_sequence(&long).unwrap_err().contains("at most 400"));
    assert!(clean_sequence(&"A".repeat(MAX_LENGTH)).is_ok());
    // DNA letters are mostly amino acid letters: it is the person's job to translate first
    assert!(clean_sequence("ACGTACGTACGTACGT").is_ok());
}

#[test]
fn the_sequence_goes_out_as_the_body_and_the_structure_comes_back() {
    let Some(curl) = curl() else { return };
    let (url, server) = serve("200 OK", FAKE_PDB);
    let pdb = predict(&curl, &url, PEPTIDE).expect("a structure");
    assert_eq!(pdb, FAKE_PDB.trim_end_matches('\n').to_string() + "\n");
    let request = server.join().expect("the server");
    assert!(request.starts_with("POST /fold"), "{request}");
    assert!(
        request.ends_with(&String::from_utf8_lossy(PEPTIDE).into_owned()),
        "the body is the sequence"
    );
}

#[test]
fn a_refusal_a_plain_answer_and_no_server_each_give_a_reason() {
    let Some(curl) = curl() else { return };
    let (url, server) = serve("400 Bad Request", "sequence too long");
    let why = predict(&curl, &url, PEPTIDE).unwrap_err();
    server.join().expect("the server");
    assert!(
        why.contains("400") && why.contains("sequence too long"),
        "{why}"
    );

    let (url, server) = serve("200 OK", "<html>maintenance</html>");
    let why = predict(&curl, &url, PEPTIDE).unwrap_err();
    server.join().expect("the server");
    assert!(
        why.contains("did not return a structure") && why.contains("maintenance"),
        "{why}"
    );

    // nothing listens: take a free port and close it
    let closed = {
        let listener = TcpListener::bind("127.0.0.1:0").expect("a free port");
        format!("http://{}/fold", listener.local_addr().expect("an address"))
    };
    let why = predict(&curl, &closed, PEPTIDE).unwrap_err();
    assert!(why.starts_with("The request failed"), "{why}");
}

#[test]
fn a_missing_curl_is_reported() {
    let why = predict(
        std::path::Path::new("/nonexistent/curl"),
        "http://127.0.0.1:1/",
        PEPTIDE,
    )
    .unwrap_err();
    assert!(why.contains("Could not start curl"), "{why}");
}

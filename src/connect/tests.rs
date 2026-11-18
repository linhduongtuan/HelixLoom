use super::curl::{wait_needed, Curl};
use super::*;
use std::cell::RefCell;
use std::sync::atomic::{AtomicU32, Ordering};
use std::time::Duration;
use web_time::Instant;

static NEXT: AtomicU32 = AtomicU32::new(0);

struct Folder(PathBuf);

impl Folder {
    fn new() -> Folder {
        let path = std::env::temp_dir().join(format!(
            "helixloom-connect-test-{}-{}",
            std::process::id(),
            NEXT.fetch_add(1, Ordering::Relaxed)
        ));
        std::fs::create_dir_all(&path).expect("a folder");
        Folder(path)
    }
}

impl Drop for Folder {
    fn drop(&mut self) {
        let _ = std::fs::remove_dir_all(&self.0);
    }
}

/// A fetcher that answers from a table and remembers what it was asked.
struct Canned {
    replies: Vec<(String, Result<Reply, String>)>,
    asked: RefCell<Vec<String>>,
}

impl Canned {
    fn new(replies: Vec<(String, Result<Reply, String>)>) -> Canned {
        Canned {
            replies,
            asked: RefCell::new(Vec::new()),
        }
    }

    fn with(url: &str, status: u16, body: &str) -> Canned {
        Canned::new(vec![(
            url.to_string(),
            Ok(Reply {
                status,
                body: body.to_string(),
            }),
        )])
    }

    fn asked(&self) -> Vec<String> {
        self.asked.borrow().clone()
    }
}

impl Fetcher for Canned {
    fn get(&self, url: &str) -> Result<Reply, String> {
        self.asked.borrow_mut().push(url.to_string());
        self.replies
            .iter()
            .find(|(u, _)| u == url)
            .map(|(_, r)| r.clone())
            .unwrap_or(Ok(Reply {
                status: 404,
                body: String::new(),
            }))
    }
}

const GENBANK: &str = "LOCUS       AB123456                  24 bp    DNA     linear   BCT 01-JAN-2020\nDEFINITION  A test.\nACCESSION   AB123456\nFEATURES             Location/Qualifiers\n     CDS             1..12\n                     /gene=\"x\"\nORIGIN\n        1 atgaaatttg gccctaaaat ttaa\n//\n";
const FASTA: &str = ">sp|P69905|HBA_HUMAN Hemoglobin subunit alpha\nMVLSPADKTNVKAAWGKVGAHAGEYGAEALERMFLSFPTTKTYFPHF\n";
const PDB: &str = "HEADER    OXYGEN TRANSPORT\nATOM      1  N   VAL A   1      6.204  16.869   4.854  1.00 49.05           N\nEND\n";

fn request(source: Source, id: &str) -> Request {
    Request::check(source, id).expect("an accession")
}

// --- accessions ---------------------------------------------------------------------------------------------------------------

#[test]
fn accessions_of_each_database_are_taken_in_the_form_that_database_gives() {
    for (source, ids) in [
        (
            Source::NcbiNucleotide,
            &[
                "NC_000913.3",
                "U00096.3",
                "J01415",
                "CP000001.1",
                "NM_000518.5",
            ][..],
        ),
        (
            Source::NcbiProtein,
            &["NP_414542.1", "AAA12345.1", "P69905"][..],
        ),
        (
            Source::UniProt,
            &["P69905", "Q8N158", "O00533", "A2BC19", "A0A024R161"][..],
        ),
        (Source::AlphaFold, &["P69905", "A0A024R161"][..]),
        (Source::Pdb, &["1A3N", "4HHB", "7XYZ"][..]),
    ] {
        for id in ids {
            assert!(Request::check(source, id).is_ok(), "{source:?} {id}");
        }
    }
}

#[test]
fn what_is_not_an_accession_is_refused_with_what_one_looks_like() {
    for (source, ids) in [
        (
            Source::NcbiNucleotide,
            &[
                "",
                "NC_",
                "123456",
                "NC_000913.",
                "NC_000913.1234",
                "ABCDEFG123456",
                "NC 000913",
                "NC_00091",
            ][..],
        ),
        (
            Source::UniProt,
            &["P6990", "P699055", "1P69905", "PPPPPP", "", "A0A024R16"][..],
        ),
        (Source::Pdb, &["A3N1", "1A3", "1A3NN", "1A_N", ""][..]),
    ] {
        for id in ids {
            // "NC_00091" has five digits, the least; the others are all too short, too long or the wrong shape
            let r = Request::check(source, id);
            if *id == "NC_00091" {
                assert!(r.is_ok());
            } else {
                let why = r.unwrap_err();
                assert!(
                    why.contains("is not an accession of") && why.contains("look like"),
                    "{id}: {why}"
                );
            }
        }
    }
}

#[test]
fn text_that_could_add_to_a_request_is_never_an_accession() {
    for hostile in [
        "NC_000913.3&db=gene",
        "NC_000913.3\nHost: evil",
        "P69905/../../etc",
        "1A3N;ls",
        "P69905 --output /tmp/x",
        "$(id)",
        "NC_٠٠٠٠٠٠",
        "http://evil/NC_000913.3",
        "P69905%0d%0a",
    ] {
        for source in Source::ALL {
            assert!(
                Request::check(source, hostile).is_err(),
                "{source:?} {hostile:?}"
            );
        }
    }
}

#[test]
fn an_accession_is_trimmed_and_put_in_the_databases_letters() {
    let r = request(Source::UniProt, "  p69905\n");
    assert_eq!(r.id, "P69905");
    assert_eq!(
        request(Source::NcbiNucleotide, "nc_000913.3").id,
        "NC_000913.3"
    );
}

#[test]
fn each_request_goes_to_the_address_of_its_database() {
    let url = |s, id: &str| request(s, id).url();
    assert_eq!(
        url(Source::NcbiNucleotide, "NC_000913.3"),
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000913.3&rettype=gb&retmode=text&tool=HelixLoom"
    );
    assert_eq!(
        url(Source::NcbiProtein, "NP_414542.1"),
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=NP_414542.1&rettype=fasta&retmode=text&tool=HelixLoom"
    );
    assert_eq!(
        url(Source::UniProt, "P69905"),
        "https://rest.uniprot.org/uniprotkb/P69905.fasta"
    );
    assert_eq!(
        url(Source::Pdb, "1a3n"),
        "https://files.rcsb.org/download/1A3N.pdb"
    );
    assert_eq!(
        url(Source::AlphaFold, "P69905"),
        "https://alphafold.ebi.ac.uk/api/prediction/P69905"
    );
}

#[test]
fn every_database_has_a_licence_line_a_format_and_an_example_that_passes() {
    for s in Source::ALL {
        assert!(s.licence().contains("http"), "{s:?}");
        assert!(!s.label().is_empty());
        let first = s.example().split(" or ").next().unwrap();
        assert!(Request::check(s, first).is_ok(), "{s:?}: {first}");
    }
}

// --- fetching and keeping -------------------------------------------------------------------------------------------------------

#[test]
fn a_record_is_fetched_with_where_it_came_from_and_under_what_licence() {
    let r = request(Source::NcbiNucleotide, "AB123456.1");
    let net = Canned::with(&r.url(), 200, GENBANK);
    let got = fetch(&r, &net, None, false).unwrap();
    assert_eq!(got.body, GENBANK);
    assert!(!got.from_cache);
    let p = &got.provenance;
    assert_eq!(
        (p.source.as_str(), p.id.as_str()),
        ("NCBI nucleotide", "AB123456.1")
    );
    assert_eq!(p.url, r.url());
    assert_eq!(p.licence, Source::NcbiNucleotide.licence());
    assert_eq!(p.digest, digest(GENBANK.as_bytes()));
    assert_eq!(p.fetched.len(), 10);
    assert!(p.fetched.starts_with("20"), "{}", p.fetched);
    assert_eq!(net.asked(), [r.url()]);
}

#[test]
fn a_record_that_was_kept_is_not_asked_for_again_unless_a_fresh_copy_is_wanted() {
    let f = Folder::new();
    let store = Store::open(&f.0.join("kept")).unwrap();
    let r = request(Source::UniProt, "P69905");
    let net = Canned::with(&r.url(), 200, FASTA);
    let first = fetch(&r, &net, Some(&store), false).unwrap();
    assert!(!first.from_cache);
    let second = fetch(&r, &net, Some(&store), false).unwrap();
    assert!(second.from_cache);
    assert_eq!(second.body, FASTA);
    assert_eq!(
        second.provenance, first.provenance,
        "the day it was first fetched"
    );
    assert_eq!(
        net.asked().len(),
        1,
        "nothing was asked of the database the second time"
    );
    let third = fetch(&r, &net, Some(&store), true).unwrap();
    assert!(!third.from_cache);
    assert_eq!(net.asked().len(), 2);
}

#[test]
fn a_kept_file_that_was_changed_or_lost_its_provenance_is_fetched_again() {
    let f = Folder::new();
    let store = Store::open(&f.0).unwrap();
    let r = request(Source::UniProt, "P69905");
    let net = Canned::with(&r.url(), 200, FASTA);
    fetch(&r, &net, Some(&store), false).unwrap();
    let body = f.0.join("uniprot-P69905.body");
    std::fs::write(&body, ">changed\nAAAA\n").unwrap();
    assert!(store.get(&r).is_none(), "the digest no longer matches");
    let again = fetch(&r, &net, Some(&store), false).unwrap();
    assert!(!again.from_cache);
    assert_eq!(again.body, FASTA);
    std::fs::remove_file(f.0.join("uniprot-P69905.provenance")).unwrap();
    assert!(store.get(&r).is_none());
}

#[test]
fn an_answer_that_is_not_the_record_is_said_so_in_words_for_the_person() {
    let r = request(Source::NcbiNucleotide, "ZZ999999.1");
    let say = |status: u16, body: &str| {
        fetch(&r, &Canned::with(&r.url(), status, body), None, false).unwrap_err()
    };
    assert_eq!(
        say(404, ""),
        "NCBI nucleotide has no record called ZZ999999.1."
    );
    assert_eq!(
        say(400, "Bad request"),
        "NCBI nucleotide has no record called ZZ999999.1."
    );
    assert!(say(429, "").contains("too many requests"));
    assert!(say(503, "").contains("not answering just now (it said 503)"));
    assert!(say(418, "").contains("answered 418"));
    let page = say(200, "<html><body>Error: no such id</body></html>");
    assert!(
        page.contains("answered, but not with ZZ999999.1's record: <html>"),
        "{page}"
    );
    let unreachable = Canned::new(vec![(r.url(), Err("could not resolve host".into()))]);
    let why = fetch(&r, &unreachable, None, false).unwrap_err();
    assert_eq!(
        why,
        "Could not reach NCBI nucleotide: could not resolve host"
    );
}

#[test]
fn each_format_is_checked_for_what_a_record_of_it_holds() {
    assert!(check_body(Format::GenBank, "\n LOCUS x"));
    assert!(!check_body(Format::GenBank, ">x\nAC"));
    assert!(check_body(Format::Fasta, ">x\nAC"));
    assert!(!check_body(Format::Fasta, "LOCUS"));
    assert!(check_body(Format::Pdb, "HEADER\nHETATM    1 ZN\n"));
    assert!(!check_body(Format::Pdb, "HEADER only\nEND\n"));
}

#[test]
fn the_alphafold_structure_is_found_through_its_entry_and_only_on_its_own_server() {
    let r = request(Source::AlphaFold, "P69905");
    let file = "https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_v4.pdb";
    let entry =
        format!("[{{\"entryId\":\"AF-P69905-F1\",\"pdbUrl\": \"{file}\",\"cifUrl\":\"x\"}}]");
    let net = Canned::new(vec![
        (
            r.url(),
            Ok(Reply {
                status: 200,
                body: entry,
            }),
        ),
        (
            file.to_string(),
            Ok(Reply {
                status: 200,
                body: PDB.to_string(),
            }),
        ),
    ]);
    let got = fetch(&r, &net, None, false).unwrap();
    assert_eq!(net.asked(), [r.url(), file.to_string()]);
    assert_eq!(
        got.provenance.url, file,
        "the file's address, not the entry's"
    );
    assert_eq!(got.body, PDB);

    for bad in [
        "[{\"pdbUrl\":\"https://evil.example/files/x.pdb\"}]",
        "[{\"pdbUrl\":\"https://alphafold.ebi.ac.uk/files/../../x.pdb\"}]",
        "[{\"pdbUrl\":\"https://alphafold.ebi.ac.uk/files/a b.pdb\"}]",
        "[{\"pdbUrl\":\"http://alphafold.ebi.ac.uk/files/x.pdb\"}]",
        "[{\"cifUrl\":\"x\"}]",
        "[]",
    ] {
        let net = Canned::with(&r.url(), 200, bad);
        let why = fetch(&r, &net, None, false).unwrap_err();
        assert!(why.contains("names no structure file"), "{bad}: {why}");
        assert_eq!(net.asked().len(), 1, "nothing else was asked: {bad}");
    }
}

#[test]
fn when_the_alphafold_entry_service_refuses_the_files_are_tried_by_their_names_newest_first() {
    let r = request(Source::AlphaFold, "P69905");
    let file = |v: &str| format!("https://alphafold.ebi.ac.uk/files/AF-P69905-F1-model_{v}.pdb");
    let refuse = (
        r.url(),
        Ok(Reply {
            status: 403,
            body: "403 Forbidden".into(),
        }),
    );
    // the newest version is there
    let net = Canned::new(vec![
        refuse.clone(),
        (
            file("v6"),
            Ok(Reply {
                status: 200,
                body: PDB.to_string(),
            }),
        ),
    ]);
    let got = fetch(&r, &net, None, false).unwrap();
    assert_eq!(net.asked(), [r.url(), file("v6")]);
    assert_eq!(
        got.provenance.url,
        file("v6"),
        "the address that gave the file"
    );
    // an older one
    let net = Canned::new(vec![
        refuse.clone(),
        (
            file("v5"),
            Ok(Reply {
                status: 200,
                body: PDB.to_string(),
            }),
        ),
    ]);
    assert_eq!(
        fetch(&r, &net, None, false).unwrap().provenance.url,
        file("v5")
    );
    assert_eq!(net.asked(), [r.url(), file("v6"), file("v5")]);
    // none: the refusal is what is said
    let net = Canned::new(vec![refuse]);
    let why = fetch(&r, &net, None, false).unwrap_err();
    assert!(why.contains("refused the request (it said 403)"), "{why}");
    assert_eq!(net.asked().len(), 4);
    // a file of the right name that is not a structure is no record
    let net = Canned::new(vec![
        (
            r.url(),
            Ok(Reply {
                status: 500,
                body: String::new(),
            }),
        ),
        (
            file("v6"),
            Ok(Reply {
                status: 200,
                body: "<html>maintenance</html>".into(),
            }),
        ),
    ]);
    assert!(fetch(&r, &net, None, false).is_err());
}

#[test]
fn when_the_alphafold_entry_service_says_there_is_no_entry_no_file_is_looked_for() {
    let r = request(Source::AlphaFold, "P69905");
    let net = Canned::with(&r.url(), 404, "");
    let why = fetch(&r, &net, None, false).unwrap_err();
    assert_eq!(why, "AlphaFold DB has no record called P69905.");
    assert_eq!(net.asked(), [r.url()]);
}

#[test]
fn a_provenance_is_written_and_read_back_and_a_foreign_file_is_not_one() {
    let p = Provenance {
        source: "UniProt".into(),
        id: "P69905".into(),
        url: "https://rest.uniprot.org/uniprotkb/P69905.fasta".into(),
        fetched: "2026-11-11".into(),
        licence: "CC BY 4.0; see https://www.uniprot.org/help/license".into(),
        digest: "abc".into(),
    };
    assert_eq!(Provenance::from_text(&p.to_text()), Some(p.clone()));
    assert_eq!(Provenance::from_text("helixloom-provenance 2\n"), None);
    assert_eq!(
        Provenance::from_text("helixloom-provenance 1\nsource: x\n"),
        None
    );
    let line = p.comment();
    assert!(
        line.starts_with("Fetched from UniProt P69905 on 2026-11-11"),
        "{line}"
    );
    assert!(line.contains("Licence: CC BY 4.0"), "{line}");
}

#[test]
fn a_fetched_genbank_file_becomes_a_document_that_keeps_where_it_came_from_when_saved() {
    let r = request(Source::NcbiNucleotide, "AB123456.1");
    let got = fetch(&r, &Canned::with(&r.url(), 200, GENBANK), None, false).unwrap();
    let docs = got.documents().unwrap();
    assert_eq!(docs.len(), 1);
    assert_eq!(docs[0].len(), 24);
    assert_eq!(docs[0].features.len(), 1);
    assert!(
        docs[0]
            .comment
            .contains("Fetched from NCBI nucleotide AB123456.1"),
        "{}",
        docs[0].comment
    );
    let saved = crate::seq::io::to_genbank(&docs[0]);
    let again = crate::seq::io::parse_any(&saved, "x").unwrap().remove(0);
    assert!(
        again
            .comment
            .contains("Fetched from NCBI nucleotide AB123456.1"),
        "{saved}"
    );
    assert!(again.comment.contains("Licence:"), "{}", again.comment);
}

#[test]
fn a_fasta_record_becomes_a_document_and_a_structure_does_not() {
    let r = request(Source::UniProt, "P69905");
    let got = fetch(&r, &Canned::with(&r.url(), 200, FASTA), None, false).unwrap();
    let docs = got.documents().unwrap();
    assert_eq!(docs[0].seq.len(), 47);
    assert!(docs[0].comment.contains("UniProt P69905"));
    let r = request(Source::Pdb, "1A3N");
    let structure = fetch(&r, &Canned::with(&r.url(), 200, PDB), None, false).unwrap();
    assert!(structure
        .documents()
        .unwrap_err()
        .contains("no sequence document"));
}

// --- the date and the pace --------------------------------------------------------------------------------------------------------

#[test]
fn dates_are_counted_from_the_epoch_in_the_gregorian_calendar() {
    assert_eq!(utc_date(0), "1970-01-01");
    assert_eq!(utc_date(86_399), "1970-01-01");
    assert_eq!(utc_date(86_400), "1970-01-02");
    assert_eq!(
        utc_date(951_782_400),
        "2000-02-29",
        "a leap day of a year divisible by 400"
    );
    assert_eq!(utc_date(1_700_000_000), "2023-11-14");
    assert_eq!(
        utc_date(4_107_542_400),
        "2100-03-01",
        "2100 is not a leap year"
    );
    assert_eq!(utc_date(1_793_500_000), "2026-11-01");
}

#[test]
fn requests_to_ncbi_are_spaced_by_the_time_since_the_last() {
    let gap = Duration::from_millis(350);
    let now = Instant::now();
    assert_eq!(wait_needed(None, now, gap), Duration::ZERO);
    assert_eq!(wait_needed(Some(now), now, gap), gap);
    assert_eq!(
        wait_needed(Some(now), now + Duration::from_millis(100), gap),
        Duration::from_millis(250)
    );
    assert_eq!(
        wait_needed(Some(now), now + Duration::from_secs(1), gap),
        Duration::ZERO
    );
}

// --- curl, against a server on this computer ---------------------------------------------------------------------------------------

/// A server on this computer that answers `/ok` with 200 and anything else with 404, for `n` requests.
fn serve(n: usize) -> (String, std::thread::JoinHandle<()>) {
    use std::io::{Read, Write};
    let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
    let base = format!("http://{}", listener.local_addr().unwrap());
    let handle = std::thread::spawn(move || {
        for _ in 0..n {
            let (mut socket, _) = listener.accept().unwrap();
            let mut buffer = [0u8; 2048];
            let read = socket.read(&mut buffer).unwrap_or(0);
            let head = String::from_utf8_lossy(&buffer[..read]).into_owned();
            let (status, body) = if head.starts_with("GET /ok ") {
                ("200 OK", ">sp|X\nMKV\n")
            } else {
                ("404 Not Found", "no such thing\n")
            };
            let reply = format!(
                "HTTP/1.1 {status}\r\nContent-Length: {}\r\nConnection: close\r\n\r\n{body}",
                body.len()
            );
            let _ = socket.write_all(reply.as_bytes());
        }
    });
    (base, handle)
}

fn curl_here() -> Option<PathBuf> {
    crate::pipeline::tools::detect(&crate::pipeline::tools::CURL).map(|c| c.path)
}

#[test]
fn curl_brings_the_status_and_the_body_and_a_missing_record_is_a_status_not_a_failure() {
    let Some(program) = curl_here() else { return };
    let (base, server) = serve(2);
    let curl = Curl::for_tests(program);
    let ok = curl.get(&format!("{base}/ok")).unwrap();
    assert_eq!(
        ok,
        Reply {
            status: 200,
            body: ">sp|X\nMKV\n".into()
        }
    );
    let missing = curl.get(&format!("{base}/nothing")).unwrap();
    assert_eq!(missing.status, 404);
    server.join().unwrap();
}

#[test]
fn the_program_never_asks_for_a_plain_http_address_and_a_dead_address_is_an_error() {
    let Some(program) = curl_here() else { return };
    let strict = Curl::new(program.clone());
    let why = strict.get("http://127.0.0.1:9/ok").unwrap_err();
    assert!(
        why.to_lowercase().contains("protocol") || why.contains("not supported"),
        "{why}"
    );
    let lenient = Curl::for_tests(program);
    assert!(
        lenient.get("http://127.0.0.1:9/ok").is_err(),
        "nothing listens there"
    );
}

#[test]
fn refseq_genomes_and_contigs_and_uniprot_isoforms_are_accessions() {
    for id in ["NZ_CP009072.1", "NZ_JABCDE010000001.1", "NZ_LR134356.1"] {
        assert!(Request::check(Source::NcbiNucleotide, id).is_ok(), "{id}");
    }
    assert_eq!(request(Source::UniProt, "p69905-2").id, "P69905-2");
    assert_eq!(
        request(Source::UniProt, "P69905-2").url(),
        "https://rest.uniprot.org/uniprotkb/P69905-2.fasta"
    );
    for bad in ["P69905-", "P69905-x", "P69905-1234", "P69905--2", "-2"] {
        assert!(Request::check(Source::UniProt, bad).is_err(), "{bad}");
    }
    assert!(
        Request::check(Source::AlphaFold, "P69905-2").is_err(),
        "the AlphaFold DB has the canonical sequences"
    );
}

#[test]
fn a_pdb_entry_without_a_file_in_pdb_format_is_said_so_and_not_called_unknown() {
    let r = request(Source::Pdb, "7XYZ");
    let why = fetch(&r, &Canned::with(&r.url(), 404, ""), None, false).unwrap_err();
    assert!(
        why.contains("no file in PDB format for 7XYZ") && why.contains("mmCIF"),
        "{why}"
    );
}

#[test]
fn a_record_that_cannot_be_kept_is_still_the_record() {
    let f = Folder::new();
    let dir = f.0.join("gone");
    let store = Store::open(&dir).unwrap();
    std::fs::remove_dir_all(&dir).unwrap();
    let r = request(Source::UniProt, "P69905");
    let got = fetch(&r, &Canned::with(&r.url(), 200, FASTA), Some(&store), false).unwrap();
    assert_eq!(got.body, FASTA);
    assert!(
        got.warning
            .as_deref()
            .unwrap_or("")
            .contains("Could not keep the file"),
        "{:?}",
        got.warning
    );
    assert!(!got.from_cache);
}

/// A server on this computer that answers every request with a redirect to `to`, for `n` requests.
fn redirect(n: usize, to: &'static str) -> (String, std::thread::JoinHandle<()>) {
    use std::io::{Read, Write};
    let listener = std::net::TcpListener::bind("127.0.0.1:0").unwrap();
    let base = format!("http://{}", listener.local_addr().unwrap());
    let own = base.clone();
    let handle = std::thread::spawn(move || {
        for _ in 0..n {
            let Ok((mut socket, _)) = listener.accept() else {
                return;
            };
            let mut buffer = [0u8; 2048];
            let _ = socket.read(&mut buffer);
            let target = if to.is_empty() {
                format!("{own}/again")
            } else {
                to.to_string()
            };
            let reply = format!("HTTP/1.1 302 Found\r\nLocation: {target}\r\nContent-Length: 0\r\nConnection: close\r\n\r\n");
            let _ = socket.write_all(reply.as_bytes());
        }
    });
    (base, handle)
}

#[test]
fn curl_follows_at_most_five_redirects_and_never_to_a_file_or_another_protocol() {
    let Some(program) = curl_here() else { return };
    let curl = Curl::for_tests(program);
    // a loop of redirects ends after five, not fifty
    let (base, server) = redirect(8, "");
    let why = curl.get(&format!("{base}/start")).unwrap_err();
    assert!(why.to_lowercase().contains("redirect"), "{why}");
    drop(server);
    // a redirect to a file on this computer is refused
    let (base, server) = redirect(1, "file:///etc/hosts");
    let why = curl.get(&format!("{base}/start")).unwrap_err();
    assert!(
        why.to_lowercase().contains("protocol") || why.contains("not supported"),
        "{why}"
    );
    server.join().unwrap();
}

// --- the real databases (asked only when told to: `cargo test -- --ignored live`) ----------------------------------------------

/// One record of each database and a complex, with `curl`, through the same code the window uses. It asks the real services, so it
/// does not run with the others. The records are kept in `$HELIXLOOM_LIVE_DIR` (a folder of the temporary files when it is not
/// set), and what came is printed.
#[test]
#[ignore = "asks the real databases: run with `cargo test --bin helixloom -- --ignored live --nocapture`"]
fn live_one_record_of_each_database_and_a_complex() {
    let program = curl_here().expect("curl is installed");
    let dir = std::env::var_os("HELIXLOOM_LIVE_DIR")
        .map(PathBuf::from)
        .unwrap_or_else(|| std::env::temp_dir().join("helixloom-live"));
    let store = Store::open(&dir).unwrap();
    let curl = Curl::new(program);
    // `HELIXLOOM_LIVE_SKIP=ncbi` leaves a database out (say, one that has refused this computer)
    let skip = std::env::var("HELIXLOOM_LIVE_SKIP")
        .unwrap_or_default()
        .to_lowercase();
    let mut failed: Vec<String> = Vec::new();
    for (source, id, what) in [
        (Source::NcbiNucleotide, "L09137.2", "LOCUS"),
        (Source::UniProt, "P69905", ">sp|P69905|HBA_HUMAN"),
        (Source::AlphaFold, "P69905", "ATOM"),
        (Source::Pdb, "1UBQ", "ATOM"),
        (Source::Pdb, "1VFB", "ATOM"),
    ] {
        if skip.contains(
            &source
                .label()
                .to_lowercase()
                .split(' ')
                .next()
                .unwrap_or("")
                .to_string(),
        ) {
            println!("{:<16} {id}: left out", source.label());
            continue;
        }
        let r = request(source, id);
        let got = match fetch(&r, &curl, Some(&store), true) {
            Ok(got) => got,
            Err(why) => {
                println!("{:<16} {id}: FAILED: {why}", source.label());
                failed.push(format!("{} {id}", source.label()));
                continue;
            }
        };
        println!(
            "{:<16} {:<9} {:>7} bytes  from {}\n    fetched {}, licence {}",
            source.label(),
            id,
            got.body.len(),
            got.provenance.url,
            got.provenance.fetched,
            got.provenance.licence
        );
        assert!(
            got.body.contains(what),
            "{} {id} does not hold {what}",
            source.label()
        );
        assert!(!got.from_cache);
        assert!(got.warning.is_none(), "{:?}", got.warning);
        // a second ask is answered from the folder
        let again = fetch(&r, &curl, Some(&store), false).unwrap();
        assert!(again.from_cache && again.body == got.body);
    }
    assert!(failed.is_empty(), "these were not fetched: {failed:?}");
    println!("kept in {}", dir.display());
}

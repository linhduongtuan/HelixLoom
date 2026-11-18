use super::*;
use crate::pipeline::tools::work_folder;
use std::fs;
use std::path::PathBuf;

/// A folder with what a `trunk build` leaves: a page, a script, a program, a style in a sub-folder, and things that must
/// never be served.
fn built_page(label: &str) -> PathBuf {
    let dir = work_folder(label).expect("a folder");
    fs::write(dir.join("index.html"), "<html>HelixLoom</html>").expect("index");
    fs::write(dir.join("helixloom-1a2b.js"), "export default 1;").expect("js");
    fs::write(dir.join("helixloom-1a2b_bg.wasm"), b"\0asm\x01\0\0\0").expect("wasm");
    fs::create_dir_all(dir.join("assets")).expect("sub-folder");
    fs::write(dir.join("assets/style.css"), "body{}").expect("css");
    fs::write(dir.join(".secret"), "hidden").expect("hidden file");
    fs::write(dir.join("has space.txt"), "odd name").expect("odd name");
    dir
}

#[test]
fn a_site_serves_the_files_of_its_folder_and_nothing_else() {
    let dir = built_page("site-files");
    // a link out of the folder, and a file outside it that the link points at
    #[cfg(unix)]
    {
        let outside = work_folder("site-outside").expect("a folder");
        fs::write(outside.join("private.txt"), "not for the web").expect("outside file");
        std::os::unix::fs::symlink(outside.join("private.txt"), dir.join("link.txt"))
            .expect("a link");
        std::os::unix::fs::symlink(&outside, dir.join("linked-folder")).expect("a link");
    }
    let site = Site::load(&dir, &[]).expect("a site");
    assert_eq!(
        site.paths(),
        [
            "/",
            "/assets/style.css",
            "/helixloom-1a2b.js",
            "/helixloom-1a2b_bg.wasm",
            "/index.html",
        ]
    );
}

#[test]
fn the_root_is_the_index_page() {
    let site = Site::load(&built_page("site-root"), &[]).expect("a site");
    let head = head("GET", "/", &[]);
    let (reply, _) = site.answer(&head).expect("the page");
    assert_eq!(reply.status, 200);
    assert_eq!(&*reply.body, b"<html>HelixLoom</html>");
    assert!(reply.content_type.starts_with("text/html"));
}

#[test]
fn a_folder_that_is_not_a_built_page_is_refused() {
    let empty = work_folder("site-empty").expect("a folder");
    let why = Site::load(&empty, &[]).err().expect("refused");
    assert!(why.contains("index.html"), "{why}");
    let missing = empty.join("not-there");
    assert!(Site::load(&missing, &[]).is_err());
    let file = empty.join("a-file");
    fs::write(&file, "x").expect("a file");
    assert!(Site::load(&file, &[]).is_err());
}

#[test]
fn a_folder_with_too_many_files_is_refused() {
    let dir = built_page("site-many");
    for i in 0..super::MAX_FILES {
        fs::write(dir.join(format!("f{i}.txt")), "x").expect("a file");
    }
    let why = Site::load(&dir, &[]).err().expect("refused");
    assert!(why.contains("more than"), "{why}");
}

#[test]
fn a_browser_is_told_the_type_that_the_extension_says() {
    for (name, wanted) in [
        ("index.html", "text/html; charset=utf-8"),
        ("a.js", "text/javascript; charset=utf-8"),
        ("a.mjs", "text/javascript; charset=utf-8"),
        ("a.css", "text/css; charset=utf-8"),
        // a browser refuses to compile a program that is not sent with this exact type
        ("a_bg.wasm", "application/wasm"),
        ("m.json", "application/json"),
        // what lets a browser install the page: its own type, which the standard names for it
        ("manifest.webmanifest", "application/manifest+json"),
        // a service worker is a script, and a browser refuses one that is sent as anything else
        ("sw.js", "text/javascript; charset=utf-8"),
        ("i.png", "image/png"),
        ("i.svg", "image/svg+xml"),
        ("f.ico", "image/x-icon"),
        ("A.WASM", "application/wasm"),
        ("whatever.bin", "application/octet-stream"),
        ("no-extension", "application/octet-stream"),
    ] {
        assert_eq!(content_type(name), wanted, "{name}");
    }
}

#[test]
fn only_the_origins_that_were_named_may_call_across_origins() {
    let site = Site::with_origins(&[
        "https://lab.example.org".into(),
        "http://10.0.0.5:8080".into(),
    ]);
    assert_eq!(
        site.allowed_origin(Some("https://lab.example.org")),
        Some("https://lab.example.org")
    );
    assert_eq!(
        site.allowed_origin(Some("http://10.0.0.5:8080")),
        Some("http://10.0.0.5:8080")
    );
    for refused in [
        None,
        Some(""),
        Some("*"),
        Some("null"),
        Some("https://evil.example.org"),
        // the same host on another scheme or port, a longer host that starts the same, and another case
        Some("http://lab.example.org"),
        Some("https://lab.example.org:8443"),
        Some("https://lab.example.org.evil.net"),
        Some("https://LAB.example.org"),
    ] {
        assert_eq!(site.allowed_origin(refused), None, "{refused:?}");
    }
    // nothing named: nothing is allowed, not even everything
    assert_eq!(
        Site::default().allowed_origin(Some("https://lab.example.org")),
        None
    );
}

#[test]
fn an_allowed_origin_is_checked_and_written_in_one_spelling() {
    for (typed, kept) in [
        ("https://lab.example.org", "https://lab.example.org"),
        ("  HTTPS://Lab.Example.org  ", "https://lab.example.org"),
        ("http://10.0.0.5:8080", "http://10.0.0.5:8080"),
        // a browser leaves the default port out of the Origin it sends, so an origin spelt with it would never match
        ("https://lab.example.org:443", "https://lab.example.org"),
        ("http://lab.example.org:80", "http://lab.example.org"),
        // another scheme's default port is just a port
        ("https://lab.example.org:80", "https://lab.example.org:80"),
        ("http://lab.example.org:443", "http://lab.example.org:443"),
        ("http://localhost", "http://localhost"),
        ("http://[::1]:8080", "http://[::1]:8080"),
    ] {
        assert_eq!(parse_origin(typed).as_deref(), Ok(kept), "{typed}");
    }
    for bad in [
        "",
        "*",
        "null",
        "lab.example.org",
        "ftp://lab.example.org",
        "https://",
        "https://lab.example.org/",
        "https://lab.example.org/path",
        "https://user@lab.example.org",
        "https://lab example.org",
        "https://lab.example.org:notaport",
        "https://lab.example.org:99999",
        "https://lab.example.org?x=1",
        "https://lab.example.org\r\nX-Evil: 1",
    ] {
        assert!(parse_origin(bad).is_err(), "{bad:?}");
    }
}

#[test]
fn a_page_file_is_sent_with_its_tag_and_not_sent_again_to_a_browser_that_has_it() {
    let site = Site::load(&built_page("site-etag"), &[]).expect("a site");
    let (first, headers) = site
        .answer(&head("GET", "/helixloom-1a2b_bg.wasm", &[]))
        .expect("the file");
    assert_eq!(
        (first.status, first.content_type),
        (200, "application/wasm")
    );
    let etag = header(&headers, "ETag").expect("a tag");
    assert!(etag.starts_with('"') && etag.ends_with('"'), "{etag}");
    // the same bytes, the same tag; other bytes, another tag
    let (_, again) = site
        .answer(&head("GET", "/helixloom-1a2b_bg.wasm", &[]))
        .expect("the file");
    assert_eq!(header(&again, "ETag"), Some(etag));
    let (_, other) = site
        .answer(&head("GET", "/helixloom-1a2b.js", &[]))
        .expect("the file");
    assert_ne!(header(&other, "ETag"), Some(etag));
    // a browser that has it is told so, with no body; one that has an old one is sent the file
    let had = [("If-None-Match", etag)];
    let (same, headers) = site
        .answer(&head("GET", "/helixloom-1a2b_bg.wasm", &had))
        .expect("an answer");
    assert_eq!(same.status, 304);
    assert!(same.body.is_empty());
    assert_eq!(header(&headers, "ETag"), Some(etag));
    let listed = format!("W/\"old\", {etag}");
    let (same, _) = site
        .answer(&head(
            "GET",
            "/helixloom-1a2b_bg.wasm",
            &[("If-None-Match", listed.as_str())],
        ))
        .expect("an answer");
    assert_eq!(same.status, 304);
    let (old, _) = site
        .answer(&head(
            "GET",
            "/helixloom-1a2b_bg.wasm",
            &[("If-None-Match", "\"old\"")],
        ))
        .expect("an answer");
    assert_eq!(old.status, 200);
}

#[test]
fn only_a_read_of_a_file_that_is_in_the_page_is_answered_without_a_token() {
    let site = Site::load(
        &built_page("site-gets"),
        &["https://lab.example.org".into()],
    )
    .expect("a site");
    // not a file of the page, or not a read, or an address that only looks like one
    for (method, target) in [
        ("GET", "/nope"),
        ("GET", "/v1/status"),
        ("GET", "/.secret"),
        // (a name with a space cannot be asked for as it is: the request line would not parse)
        ("GET", "/has%20space.txt"),
        ("GET", "/has+space.txt"),
        ("GET", "/../Cargo.toml"),
        ("GET", "/assets/../index.html"),
        ("GET", "/%2e%2e/Cargo.toml"),
        ("GET", "//index.html"),
        ("GET", "/index.html/"),
        ("GET", "/INDEX.HTML"),
        ("POST", "/index.html"),
        ("DELETE", "/index.html"),
        ("PUT", "/index.html"),
    ] {
        assert!(
            site.answer(&head(method, target, &[])).is_none(),
            "{method} {target}"
        );
    }
    // a query is not part of the name
    assert!(site.answer(&head("GET", "/index.html?v=3", &[])).is_some());
}

#[test]
fn a_preflight_from_an_allowed_origin_is_answered_with_what_the_page_may_send() {
    let site = Site::with_origins(&["https://lab.example.org".into()]);
    let ask = [
        ("Origin", "https://lab.example.org"),
        ("Access-Control-Request-Method", "POST"),
        (
            "Access-Control-Request-Headers",
            "authorization,content-type",
        ),
    ];
    let (reply, headers) = site
        .answer(&head("OPTIONS", "/v1/jobs", &ask))
        .expect("an answer");
    assert_eq!(reply.status, 204);
    assert!(reply.body.is_empty());
    let methods = header(&headers, "Access-Control-Allow-Methods").expect("methods");
    for m in ["GET", "POST", "DELETE"] {
        assert!(methods.contains(m), "{methods}");
    }
    let allowed = header(&headers, "Access-Control-Allow-Headers").expect("headers");
    assert!(
        allowed.contains("authorization") && allowed.contains("content-type"),
        "{allowed}"
    );
    assert!(header(&headers, "Access-Control-Max-Age").is_some());
    // a browser that is about to reach a private address asks for this as well
    assert_eq!(
        header(&headers, "Access-Control-Allow-Private-Network"),
        None
    );
    let private = [
        ("Origin", "https://lab.example.org"),
        ("Access-Control-Request-Method", "GET"),
        ("Access-Control-Request-Private-Network", "true"),
    ];
    let (_, headers) = site
        .answer(&head("OPTIONS", "/v1/status", &private))
        .expect("an answer");
    assert_eq!(
        header(&headers, "Access-Control-Allow-Private-Network"),
        Some("true")
    );
    // from an origin that is not allowed, or with no question about a method, there is no answer to give
    for ask in [
        vec![
            ("Origin", "https://evil.example.org"),
            ("Access-Control-Request-Method", "POST"),
        ],
        vec![("Access-Control-Request-Method", "POST")],
        vec![("Origin", "https://lab.example.org")],
    ] {
        assert!(
            site.answer(&head("OPTIONS", "/v1/jobs", &ask)).is_none(),
            "{ask:?}"
        );
    }
}

#[test]
fn the_cors_headers_name_the_one_origin_and_never_everyone() {
    let headers = cors_headers("https://lab.example.org");
    assert_eq!(
        header(&headers, "Access-Control-Allow-Origin"),
        Some("https://lab.example.org")
    );
    // the answer depends on who asked, so a cache must not hand it to somebody else
    assert_eq!(header(&headers, "Vary"), Some("Origin"));
    assert!(headers.iter().all(|(_, v)| v != "*"));
}

/// A request's header as the server reads it.
fn head(method: &str, target: &str, headers: &[(&str, &str)]) -> Head {
    let mut text = format!("{method} {target} HTTP/1.1\r\n");
    for (n, v) in headers {
        text.push_str(&format!("{n}: {v}\r\n"));
    }
    text.push_str("\r\n");
    crate::remote::http::read_head(&mut std::io::Cursor::new(text.into_bytes()))
        .expect("a request header")
}

fn header<'a>(headers: &'a Headers, name: &str) -> Option<&'a str> {
    headers
        .iter()
        .find(|(n, _)| n.eq_ignore_ascii_case(name))
        .map(|(_, v)| v.as_str())
}

#[test]
fn a_folder_called_v1_would_shadow_the_server_and_is_refused() {
    let dir = built_page("site-v1");
    fs::create_dir_all(dir.join("v1")).expect("a folder");
    fs::write(dir.join("v1/status"), "fake").expect("a file");
    let why = Site::load(&dir, &[]).err().expect("refused");
    assert!(why.contains("v1") && why.contains("token"), "{why}");
}

#[test]
fn a_file_that_is_too_large_to_keep_in_memory_is_refused_without_being_read() {
    let dir = built_page("site-large");
    // a sparse file: it is as long as it says, and takes no room to make
    let big = fs::File::create(dir.join("big.wasm")).expect("a file");
    big.set_len(super::MAX_FILE + 1).expect("a length");
    let why = Site::load(&dir, &[]).err().expect("refused");
    assert!(why.contains("too large"), "{why}");
}

#[test]
fn only_so_many_folders_down_are_looked_at() {
    let dir = built_page("site-depth");
    let deep = |levels: usize| {
        let mut path = dir.clone();
        for level in 0..levels {
            path.push(format!("d{level}"));
        }
        fs::create_dir_all(&path).expect("folders");
        fs::write(path.join("f.txt"), "x").expect("a file");
    };
    deep(super::MAX_DEPTH);
    deep(super::MAX_DEPTH + 1);
    let paths = Site::load(&dir, &[]).expect("a site").paths();
    let at = |levels: usize| {
        let folders: Vec<String> = (0..levels).map(|l| format!("d{l}")).collect();
        format!("/{}/f.txt", folders.join("/"))
    };
    assert!(paths.contains(&at(super::MAX_DEPTH)), "{paths:?}");
    assert!(!paths.contains(&at(super::MAX_DEPTH + 1)), "{paths:?}");
}

#[test]
fn a_file_of_the_page_is_sent_to_only_so_many_at_once_and_the_places_are_given_back() {
    let site = Site::default();
    let held: Vec<_> = (0..super::MAX_PAGE_SENDS)
        .map(|_| site.take_send_place().expect("a place"))
        .collect();
    assert!(
        site.take_send_place().is_none(),
        "one more than there are places"
    );
    // a place is given back when its sending is over, however it ends
    drop(held);
    let again: Vec<_> = (0..super::MAX_PAGE_SENDS)
        .map(|_| site.take_send_place().expect("a place again"))
        .collect();
    assert_eq!(again.len(), super::MAX_PAGE_SENDS);
    assert!(site.take_send_place().is_none());
}

#[test]
fn the_page_is_given_a_minute_and_not_the_five_that_a_result_gets() {
    use crate::remote::http::MIN_UPLOAD_RATE;
    assert_eq!(page_time(0), std::time::Duration::from_secs(15));
    // the 13 MB program of the page: a link of 2 Mbit/s (256 kB a second) gets it in time
    let page = 13 << 20;
    assert!(
        page_time(page) >= std::time::Duration::from_secs(60),
        "{:?}",
        page_time(page)
    );
    assert!(
        page_time(page) <= std::time::Duration::from_secs(70),
        "{:?}",
        page_time(page)
    );
    // at the floor of a result (50 kB a second) it would be five minutes
    let result_floor = std::time::Duration::from_secs(60 + (page / MIN_UPLOAD_RATE) as u64);
    assert!(page_time(page) * 4 < result_floor);
}

#[test]
fn a_file_of_the_page_that_holds_a_secret_is_found_by_its_path() {
    let dir = built_page("site-secret");
    fs::write(dir.join("notes.txt"), "token: sec-ret-1234567890").expect("a file");
    fs::create_dir_all(dir.join("sub")).expect("a folder");
    fs::write(dir.join("sub/more.js"), "var t='sec-ret-1234567890'").expect("a file");
    let site = Site::load(&dir, &[]).expect("a site");
    // the first by path, and never the alias of the index
    assert_eq!(
        site.file_containing(b"sec-ret-1234567890").as_deref(),
        Some("/notes.txt")
    );
    assert_eq!(
        site.file_containing(b"HelixLoom").as_deref(),
        Some("/index.html")
    );
    assert_eq!(site.file_containing(b"not-in-the-page"), None);
    assert_eq!(site.file_containing(b""), None);
}

#[test]
fn a_file_of_the_page_may_not_be_framed_and_says_nothing_of_where_it_is_linked_from() {
    let site = Site::load(&built_page("site-frame"), &[]).expect("a site");
    let (_, headers) = site
        .answer(&head("GET", "/index.html", &[]))
        .expect("the page");
    assert_eq!(header(&headers, "X-Frame-Options"), Some("DENY"));
    assert_eq!(header(&headers, "Referrer-Policy"), Some("no-referrer"));
}

#[test]
fn an_answer_depends_on_who_asked_whenever_any_origin_was_named() {
    let site = Site::with_origins(&["https://lab.example.org".into()]);
    let named = site.origin_headers(Some("https://lab.example.org"));
    assert_eq!(
        header(&named, "Access-Control-Allow-Origin"),
        Some("https://lab.example.org")
    );
    assert_eq!(header(&named, "Vary"), Some("Origin"));
    // another origin, or none: no permission, but a cache in front must still not hand this answer to the named one
    for asked in [Some("https://evil.example.org"), None] {
        let other = site.origin_headers(asked);
        assert_eq!(
            header(&other, "Access-Control-Allow-Origin"),
            None,
            "{asked:?}"
        );
        assert_eq!(header(&other, "Vary"), Some("Origin"), "{asked:?}");
    }
    // nothing named: nothing to vary on
    assert!(Site::default()
        .origin_headers(Some("https://lab.example.org"))
        .is_empty());
}

/// Bytes that pass for the gzip of `file`: the two bytes a gzip file begins with, something, and at the end the checksum and
/// the length of the file (the server checks no more: it has nothing to unpack with).
fn packed_copy(file: &[u8]) -> Vec<u8> {
    let mut bytes = vec![0x1f, 0x8b, 8, 0, 0, 0, 0, 0, 2, 255];
    bytes.extend_from_slice(b"packed");
    bytes.extend_from_slice(&crc32(file).to_le_bytes());
    bytes.extend_from_slice(&(file.len() as u32).to_le_bytes());
    bytes
}

/// A built page whose program, page and style have a packed copy beside them, and whose script has one left over from
/// another build: of other bytes, though as many.
fn packed_page(label: &str) -> PathBuf {
    let dir = built_page(label);
    for name in ["helixloom-1a2b_bg.wasm", "index.html", "assets/style.css"] {
        let file = fs::read(dir.join(name)).expect("the file");
        fs::write(dir.join(format!("{name}.gz")), packed_copy(&file)).expect("a copy");
    }
    fs::write(
        dir.join("helixloom-1a2b.js.gz"),
        packed_copy(b"export default 2;"),
    )
    .expect("a stale copy");
    dir
}

#[test]
fn the_checksum_is_the_one_a_gzip_file_ends_with() {
    // the answers every table of CRC-32 gives
    assert_eq!(crc32(b""), 0);
    assert_eq!(crc32(b"123456789"), 0xcbf4_3926);
    assert_eq!(crc32(b"a"), 0xe8b7_be43);
    assert_eq!(
        crc32(b"The quick brown fox jumps over the lazy dog"),
        0x414f_a339
    );
}

#[test]
fn a_browser_that_takes_gzip_is_sent_the_packed_copy_and_any_other_the_file_itself() {
    let dir = packed_page("site-packed");
    let site = Site::load(&dir, &[]).expect("a site");
    let wasm = "/helixloom-1a2b_bg.wasm";
    let plain_bytes = fs::read(dir.join("helixloom-1a2b_bg.wasm")).expect("the program");
    let packed_bytes = packed_copy(&plain_bytes);

    // no word about what it takes: the file itself, and a note that the answer depends on that word
    let (plain, headers) = site.answer(&head("GET", wasm, &[])).expect("the file");
    assert_eq!(
        (plain.status, plain.content_type),
        (200, "application/wasm")
    );
    assert_eq!(&plain.body[..], &plain_bytes[..]);
    assert_eq!(header(&headers, "Content-Encoding"), None);
    assert_eq!(header(&headers, "Vary"), Some("Accept-Encoding"));
    assert_eq!(header(&headers, "ETag"), Some(etag(&plain_bytes).as_str()));

    // it takes gzip: the copy, said to be gzip, with the type of the file it is a copy of and the tag of its own bytes
    let takes = [("Accept-Encoding", "br, gzip;q=0.8, deflate")];
    let (packed, headers) = site.answer(&head("GET", wasm, &takes)).expect("the copy");
    assert_eq!(
        (packed.status, packed.content_type),
        (200, "application/wasm")
    );
    assert_eq!(&packed.body[..], &packed_bytes[..]);
    assert_eq!(header(&headers, "Content-Encoding"), Some("gzip"));
    assert_eq!(header(&headers, "Vary"), Some("Accept-Encoding"));
    let packed_tag = etag(&packed_bytes);
    assert_eq!(header(&headers, "ETag"), Some(packed_tag.as_str()));
    assert_ne!(packed_tag, etag(&plain_bytes));

    // a browser that has the copy is told so; the tag of the file itself is not the copy's
    let had = [
        ("Accept-Encoding", "gzip"),
        ("If-None-Match", packed_tag.as_str()),
    ];
    let (same, headers) = site.answer(&head("GET", wasm, &had)).expect("an answer");
    assert_eq!(same.status, 304);
    assert!(same.body.is_empty());
    assert_eq!(header(&headers, "Content-Encoding"), Some("gzip"));
    let plain_tag = etag(&plain_bytes);
    let other = [
        ("Accept-Encoding", "gzip"),
        ("If-None-Match", plain_tag.as_str()),
    ];
    let (sent, _) = site.answer(&head("GET", wasm, &other)).expect("an answer");
    assert_eq!(sent.status, 200);

    // the page itself, by its name and as the root, and a file of a folder in the folder
    for (path, kind) in [
        ("/", "text/html; charset=utf-8"),
        ("/index.html", "text/html; charset=utf-8"),
        ("/assets/style.css", "text/css; charset=utf-8"),
    ] {
        let (file, headers) = site.answer(&head("GET", path, &takes)).expect("the file");
        assert_eq!(file.content_type, kind, "{path}");
        assert_eq!(header(&headers, "Content-Encoding"), Some("gzip"), "{path}");
        assert!(file.body.starts_with(&[0x1f, 0x8b]), "{path}");
    }
    // a copy that is sent in the place of its file is not a file of the page by its own name
    assert!(site
        .answer(&head("GET", "/index.html.gz", &takes))
        .is_none());
    assert!(!site.paths().iter().any(|p| p == "/assets/style.css.gz"));
}

#[test]
fn a_copy_that_is_not_the_files_and_a_file_without_a_copy_are_sent_as_they_are() {
    let dir = packed_page("site-stale-copy");
    // the page of another build beside this one's: as long, since the names in a page are as long in every build
    let page = fs::read(dir.join("index.html")).expect("the page");
    let mut older = page.clone();
    older[7] ^= 1;
    fs::write(dir.join("index.html.gz"), packed_copy(&older)).expect("a stale copy");
    fs::remove_file(dir.join("assets/style.css.gz")).expect("removed");
    let site = Site::load(&dir, &[]).expect("a site");
    let takes = [("Accept-Encoding", "gzip")];
    // neither the script's copy nor the page's is of the bytes beside it: the files are sent, and nothing is said of gzip
    for (path, bytes) in [
        ("/helixloom-1a2b.js", &b"export default 1;"[..]),
        ("/index.html", &page[..]),
        ("/", &page[..]),
        ("/assets/style.css", &b"body{}"[..]),
    ] {
        let (file, headers) = site.answer(&head("GET", path, &takes)).expect("the file");
        assert_eq!(&file.body[..], bytes, "{path}");
        assert_eq!(header(&headers, "Content-Encoding"), None, "{path}");
        assert_eq!(header(&headers, "Vary"), None, "{path}");
    }
    // a packed file that is no copy of what is beside it is a file of the folder, sent as what it is
    let (copy, headers) = site
        .answer(&head("GET", "/index.html.gz", &takes))
        .expect("the file by its name");
    assert_eq!(copy.content_type, "application/gzip");
    assert_eq!(header(&headers, "Content-Encoding"), None);
}

#[test]
fn a_packed_copy_begins_as_gzip_does_and_ends_with_the_checksum_and_the_length_of_the_file() {
    let file = b"forty bytes of a file, or thereabouts...";
    assert!(is_packed_copy(&packed_copy(file), file));
    // other bytes of the same length; a byte more; another beginning; too short to say anything; nothing
    let mut other = *file;
    other[3] ^= 0x20;
    assert!(!is_packed_copy(&packed_copy(file), &other));
    assert!(!is_packed_copy(&packed_copy(file), &file[..file.len() - 1]));
    let mut wrong = packed_copy(file);
    wrong[1] = 0x8c;
    assert!(!is_packed_copy(&wrong, file));
    assert!(!is_packed_copy(&[0x1f, 0x8b, 1, 2, 3, 4, 5], file));
    assert!(!is_packed_copy(&[], b""));
    // the length and the checksum each have to be right: one of them swapped for another file's is not enough
    let mut length_only = packed_copy(file);
    let tail = length_only.len() - 8;
    length_only[tail..tail + 4].copy_from_slice(&crc32(&other).to_le_bytes());
    assert!(!is_packed_copy(&length_only, file));
    let mut checksum_only = packed_copy(file);
    let end = checksum_only.len() - 4;
    checksum_only[end..].copy_from_slice(&(file.len() as u32 + 1).to_le_bytes());
    assert!(!is_packed_copy(&checksum_only, file));
    // an empty file's copy says that nothing was packed, and a copy may be longer than its file
    assert!(is_packed_copy(&packed_copy(b""), b""));
    assert!(packed_copy(b"ab").len() > 2 && is_packed_copy(&packed_copy(b"ab"), b"ab"));
}

#[test]
fn a_request_takes_gzip_when_it_lists_it_and_not_with_a_weight_of_nothing() {
    let takes = |listed: &str| takes_gzip(&head("GET", "/", &[("Accept-Encoding", listed)]));
    assert!(takes("gzip"));
    assert!(takes("GZIP"));
    assert!(takes("deflate, gzip, br"));
    assert!(takes("br;q=1.0, gzip;q=0.5"));
    assert!(takes("gzip ; q=0.001"));
    assert!(!takes("gzip;q=0"));
    assert!(!takes("gzip;q=0.0, br"));
    assert!(!takes("br, deflate"));
    assert!(!takes("identity"));
    assert!(!takes("x-gzip-like"));
    assert!(!takes(""));
    assert!(!takes_gzip(&head("GET", "/", &[])));
}

#[test]
fn a_packed_file_for_a_page_of_another_origin_says_that_it_depends_on_both() {
    let dir = packed_page("site-packed-origin");
    let site = Site::load(&dir, &["https://lab.example.org".to_string()]).expect("a site");
    let asked = head(
        "GET",
        "/helixloom-1a2b_bg.wasm",
        &[
            ("Accept-Encoding", "gzip"),
            ("Origin", "https://lab.example.org"),
        ],
    );
    let (_, mut headers) = site.answer(&asked).expect("the copy");
    // the server adds what it says to every answer of a named origin beside what the file's answer says
    headers.extend(site.origin_headers(asked.header("origin")));
    let varies: Vec<&str> = headers
        .iter()
        .filter(|(name, _)| name.eq_ignore_ascii_case("Vary"))
        .map(|(_, value)| value.as_str())
        .collect();
    assert!(
        varies.contains(&"Accept-Encoding") && varies.contains(&"Origin"),
        "{varies:?}"
    );
}

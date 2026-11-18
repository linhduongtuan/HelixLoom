#!/bin/sh
# Make the web version of HelixLoom: a folder of a few files that any static file server can host (docs/WEB.md), and that
# `helixloom serve --web-dir FOLDER` serves beside the jobs: the page, its script, the program, a style sheet, and what
# lets a browser install the page and open it without a network (a manifest, two icons, a service worker). Run it from
# anywhere:
#
#     scripts/build-web.sh [FOLDER]        (FOLDER: where the files go; dist of the repository when left out)
#
# It does what `trunk build --release` does, without trunk: builds the program for the browser, lets wasm-bindgen write the
# script that loads it, and names the files after what is in them, so that a browser never keeps an old program beside a
# new page. It needs the Rust target wasm32-unknown-unknown and wasm-bindgen in the version Cargo.lock names (on the PATH,
# in cargo's bin folder, or the one trunk keeps in its cache). Nothing is downloaded: when something is missing the script
# says what, and stops.
# The program is not made smaller with wasm-opt, which trunk runs (packed, that gains little): a packed copy (.gz) is put
# beside each file instead, for the server to send.
set -eu

# the folder is the caller's to name: a path that is not whole is from where the script was called
OUT=${1:-}
case "$OUT" in
    "") ;;
    /*) ;;
    *) OUT="$(pwd)/$OUT" ;;
esac
# so is the folder of files that stand in for a build (below)
PREBUILT=${HELIXLOOM_WEB_PREBUILT:-}
case "$PREBUILT" in
    "" | /*) ;;
    *) PREBUILT="$(pwd)/$PREBUILT" ;;
esac
cd "$(dirname "$0")/.."
OUT=${OUT:-$(pwd)/dist}

# What this script wrote into the folder before is taken away from it, and some of that has the names of the script's own
# sources: the folder must not be one of the two the script reads from.
ROOT=$(pwd -P)
if [ -d "$OUT" ]; then
    case "$(cd "$OUT" && pwd -P)" in
        "$ROOT" | "$ROOT/web")
            echo "The web version cannot be put in $OUT: the script's own sources are there. Name another folder." >&2
            exit 1
            ;;
    esac
fi

WORK=$(mktemp -d)
trap 'rm -rf "$WORK"' EXIT
# The tests of this script (tests/web.rs) give it the two files wasm-bindgen writes, as small stand-ins, in the folder
# HELIXLOOM_WEB_PREBUILT: nothing is built then, and what is tested is how the page is put together from them. It says
# so, since a page made that way by accident (the setting left in a shell) would be a page of old files.
if [ -n "$PREBUILT" ]; then
    echo "HELIXLOOM_WEB_PREBUILT is set: nothing is built, the page is put together from the files in $PREBUILT." >&2
    cp "$PREBUILT/helixloom.js" "$PREBUILT/helixloom_bg.wasm" "$WORK/"
else
    # a Rust from Homebrew has no wasm32 standard library, whatever `rustup target list` says: rustup's own toolchain goes first
    if command -v rustup >/dev/null 2>&1; then
        CARGO=$(rustup which cargo 2>/dev/null || true)
        if [ -n "$CARGO" ]; then
            PATH="$(dirname "$CARGO"):$PATH"
            export PATH
        fi
        if ! rustup target list --installed 2>/dev/null | grep -q '^wasm32-unknown-unknown$'; then
            echo "The Rust target for the browser is not installed: rustup target add wasm32-unknown-unknown" >&2
            exit 1
        fi
    fi

    # wasm-bindgen must be the very version the program was built against
    WANT=$(awk '/^name = "wasm-bindgen"$/ { getline; gsub(/version = |"/, ""); print; exit }' Cargo.lock)
    BINDGEN=""
    for candidate in \
        "$(command -v wasm-bindgen 2>/dev/null || true)" \
        "${CARGO_HOME:-$HOME/.cargo}/bin/wasm-bindgen" \
        "$HOME/Library/Caches/dev.trunkrs.trunk/wasm-bindgen-$WANT/wasm-bindgen" \
        "${XDG_CACHE_HOME:-$HOME/.cache}/trunk/wasm-bindgen-$WANT/wasm-bindgen"; do
        if [ -n "$candidate" ] && [ -x "$candidate" ] && [ "$("$candidate" --version 2>/dev/null)" = "wasm-bindgen $WANT" ]; then
            BINDGEN=$candidate
            break
        fi
    done
    if [ -z "$BINDGEN" ]; then
        echo "wasm-bindgen $WANT was not found (on the PATH, in cargo's bin folder or in trunk's cache, which the system may clean away): cargo install --locked wasm-bindgen-cli --version $WANT" >&2
        exit 1
    fi

    # the folders of the computer that builds it stay out of the file (they are in the messages of a fault otherwise)
    RUSTFLAGS="${RUSTFLAGS:-} --remap-path-prefix=$HOME=~"
    export RUSTFLAGS
    cargo build --release --target wasm32-unknown-unknown --no-default-features --bin helixloom

    "$BINDGEN" --target web --no-typescript --out-dir "$WORK" --out-name helixloom \
        target/wasm32-unknown-unknown/release/helixloom.wasm
fi

# each file is named after what is in it (a checksum; the program that makes one has two names)
if command -v shasum >/dev/null 2>&1; then
    sum() { shasum -a 256 "$1"; }
elif command -v sha256sum >/dev/null 2>&1; then
    sum() { sha256sum "$1"; }
else
    echo "Neither shasum nor sha256sum was found: the files cannot be named after what is in them." >&2
    exit 1
fi
tag() {
    name=$(sum "$1" | cut -c1-16)
    # a name of sixteen letters, or nothing was summed: files with no name of their own would never be seen to change
    case "$name" in
        ????????????????) printf '%s' "$name" ;;
        *)
            echo "The checksum of $1 could not be made." >&2
            exit 1
            ;;
    esac
}
WASM=helixloom-$(tag "$WORK/helixloom_bg.wasm")_bg.wasm
CSS=style-$(tag web/style.css).css
# the script that loads the program names it: by its new name
sed "s/helixloom_bg\.wasm/$WASM/g" "$WORK/helixloom.js" >"$WORK/named.js"
JS=helixloom-$(tag "$WORK/named.js").js

# The page is put together in a folder of its own and goes into the one that was named only when all of it is there: a
# file that is missing here stops the script before anything of the page that was there before is touched.
SITE="$WORK/site"
mkdir "$SITE"
INSTALL="manifest.webmanifest icon-192.png icon-512.png"
cp "$WORK/helixloom_bg.wasm" "$SITE/$WASM"
cp "$WORK/named.js" "$SITE/$JS"
cp web/style.css "$SITE/$CSS"
# what a browser needs to install the page: its name and colours, and its picture in two sizes (made by the program
# itself: `helixloom icon web/icon-192.png 192`). They keep their names (the page names the manifest, and the manifest the
# icons): a server has to let the browser ask whether they are still the same, as `helixloom serve` does for every file.
for file in $INSTALL; do
    cp "web/$file" "$SITE/$file"
done
cat >"$SITE/index.html" <<PAGE
<!DOCTYPE html>
<html lang="en">
  <head>
    <meta charset="utf-8" />
    <meta name="viewport" content="width=device-width, initial-scale=1" />
    <meta name="color-scheme" content="light dark" />
    <meta name="theme-color" content="#5046e5" />
    <title>HelixLoom</title>
    <link rel="manifest" href="./manifest.webmanifest" />
    <link rel="icon" type="image/png" href="./icon-192.png" />
    <link rel="apple-touch-icon" href="./icon-192.png" />
    <link rel="stylesheet" href="./$CSS" />
  </head>
  <body>
    <canvas id="helixloom_canvas"></canvas>
    <div id="loading">Loading HelixLoom…</div>
    <script type="module">
      import init from "./$JS";
      init({ module_or_path: "./$WASM" });
    </script>
    <script>
      // The page is kept for use without a network where the browser allows a page that: over https, or on this
      // computer itself. Where it does not, or will not, the page works as it always did.
      if ("serviceWorker" in navigator) {
        window.addEventListener("load", () => {
          navigator.serviceWorker.register("./sw.js").catch((why) => {
            console.warn("HelixLoom is not kept for use without a network:", why);
          });
        });
      }
    </script>
  </body>
</html>
PAGE

# The service worker keeps the page and what the page needs, by name, and is a new version whenever one of them is
# another file: its version is a checksum of the names (which are checksums of the files that have one) and of the files
# that keep their names.
{
    printf '%s %s %s\n' "$JS" "$WASM" "$CSS"
    for file in $INSTALL index.html; do
        sum "$SITE/$file" | cut -c1-64
    done
    sum web/sw.js | cut -c1-64
} >"$WORK/build.txt"
VERSION=$(tag "$WORK/build.txt")
KEPT="\"./\", \"./$JS\", \"./$WASM\", \"./$CSS\""
for file in $INSTALL; do
    KEPT="$KEPT, \"./$file\""
done
sed -e "s|__VERSION__|$VERSION|" -e "s|__FILES__|$KEPT|" web/sw.js >"$SITE/sw.js"

# a packed copy beside each file: `helixloom serve --web-dir` sends it to a browser that takes gzip (so does nginx with
# `gzip_static on`), and the program is a third of its size that way. Without gzip the files are served as they are.
if command -v gzip >/dev/null 2>&1; then
    for file in "$SITE/$WASM" "$SITE/$JS" "$SITE/$CSS" "$SITE/index.html" "$SITE/sw.js" "$SITE/manifest.webmanifest"; do
        gzip -9 -n -c "$file" >"$file.gz"
    done
else
    echo "gzip was not found: no packed copies were made, and the files are served as they are." >&2
fi

# Only what this script wrote before is taken away from the folder: the files it names after a checksum (sixteen letters,
# so a helixloom-notes.js of the person's own stays), and those that keep their names. Then the new ones go in.
mkdir -p "$OUT"
H='????????????????'
for old in "$OUT"/helixloom-$H.js "$OUT"/helixloom-${H}_bg.wasm "$OUT"/style-$H.css; do
    rm -f "$old" "$old.gz"
done
for file in index.html sw.js $INSTALL; do
    rm -f "$OUT/$file" "$OUT/$file.gz"
done
cp "$SITE"/* "$OUT"/

echo "The web version is in $OUT:"
ls -l "$OUT"
echo "Try it: (cd \"$OUT\" && python3 -m http.server 8080 --bind 127.0.0.1), then open http://localhost:8080/"

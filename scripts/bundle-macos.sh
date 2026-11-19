#!/bin/sh
# Make HelixLoom.app, and a disk image of it to hand on, from a release build. Run it on a Mac, from anywhere:
#
#     scripts/bundle-macos.sh [FOLDER]        (FOLDER: where the app goes; target/bundle of the repository when left out)
#
# The app is signed "ad hoc": it starts on the Mac that made it, and on another Mac the first start needs a right-click and
# Open, because it is not signed with a developer identity and not notarized by Apple (both need an Apple Developer account;
# DEPLOY.md says what to do with one). BUNDLE_ID names the app to the system (org.helixloom.app when not set).
#
# The app declares no kinds of file as its own: a file that is double-clicked in the Finder is handed to an app in a way the
# program does not take yet, so "Open With" is left out rather than promised. Files are opened from the app, dropped on its
# window, or named on the command line.
set -eu

# the folder for the app is the caller's to name: a path that is not whole is from where the script was called
OUT=${1:-}
case "$OUT" in
    "") ;;
    /*) ;;
    *) OUT="$(pwd)/$OUT" ;;
esac
cd "$(dirname "$0")/.."

VERSION=$(sed -n 's/^version = "\(.*\)"/\1/p' Cargo.toml | head -1)
BUNDLE_ID=${BUNDLE_ID:-org.helixloom.app}
# the name goes into a file of markup as it is: letters, digits, dots and dashes only
case "$BUNDLE_ID" in
    *[!A-Za-z0-9.-]* | "")
        echo "BUNDLE_ID may hold letters, digits, dots and dashes only: $BUNDLE_ID" >&2
        exit 1
        ;;
esac
# where cargo builds, when another place than `target` was named for it
TARGET=${CARGO_TARGET_DIR:-target}
OUT=${OUT:-$TARGET/bundle}
APP="$OUT/HelixLoom.app"

# CARGO_FLAGS adds options to the build (a build without the GPU engine: CARGO_FLAGS=--no-default-features)
# shellcheck disable=SC2086
cargo build --release ${CARGO_FLAGS:-}
BIN="$TARGET/release/helixloom"

rm -rf "$APP"
mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
cp "$BIN" "$APP/Contents/MacOS/helixloom"

# the icon is the program's own mark, drawn by the program at each size an icon holds
ICONSET="$OUT/HelixLoom.iconset"
rm -rf "$ICONSET"
mkdir -p "$ICONSET"
# the pictures are for making the icon only: gone when the script ends, however it ends
trap 'rm -rf "$ICONSET"' EXIT
for size in 16 32 128 256 512; do
    "$BIN" icon "$ICONSET/icon_${size}x${size}.png" "$size"
    "$BIN" icon "$ICONSET/icon_${size}x${size}@2x.png" "$((size * 2))"
done
iconutil -c icns "$ICONSET" -o "$APP/Contents/Resources/HelixLoom.icns"

cat > "$APP/Contents/Info.plist" <<PLIST
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0">
<dict>
    <key>CFBundleName</key>
    <string>HelixLoom</string>
    <key>CFBundleDisplayName</key>
    <string>HelixLoom</string>
    <key>CFBundleIdentifier</key>
    <string>$BUNDLE_ID</string>
    <key>CFBundleExecutable</key>
    <string>helixloom</string>
    <key>CFBundleIconFile</key>
    <string>HelixLoom</string>
    <key>CFBundlePackageType</key>
    <string>APPL</string>
    <key>CFBundleShortVersionString</key>
    <string>$VERSION</string>
    <key>CFBundleVersion</key>
    <string>$VERSION</string>
    <key>CFBundleInfoDictionaryVersion</key>
    <string>6.0</string>
    <key>LSMinimumSystemVersion</key>
    <string>11.0</string>
    <key>NSHighResolutionCapable</key>
    <true/>
    <key>NSHumanReadableCopyright</key>
    <string>MIT licence. The programs it can call and the model files are other projects' own.</string>
</dict>
</plist>
PLIST
plutil -lint "$APP/Contents/Info.plist" > /dev/null

# ad hoc: no identity, so that the app starts here; another Mac asks before the first start
codesign --force --sign - "$APP"
codesign --verify "$APP"

DMG="$OUT/HelixLoom-$VERSION.dmg"
hdiutil create -quiet -ov -volname "HelixLoom $VERSION" -srcfolder "$APP" -format UDZO "$DMG"

echo "Made $APP"
echo "Made $DMG"

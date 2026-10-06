#!/bin/bash
# Builds "Proteomics Analyzer.app" (a double-clickable launcher with an icon) next to this
# repository's launch.sh. The app points at this folder, so rebuild it if you move the folder.
#   macos/build_app.sh [destination-folder]      (default: the repository folder)
set -eu
HERE="$(cd "$(dirname "$0")/.." && pwd)"
DEST="${1:-$HERE}"
APP="$DEST/Proteomics Analyzer.app"
VERSION="$(sed -n 's/^__version__ = "\(.*\)"/\1/p' "$HERE/proteomics_analyzer/__init__.py")"
rm -rf "$APP"; mkdir -p "$APP/Contents/MacOS" "$APP/Contents/Resources"
cat > "$APP/Contents/MacOS/ProteomicsAnalyzer" <<SH
#!/bin/bash
export MODE=gui
exec /bin/bash "$HERE/launch.sh"
SH
chmod +x "$APP/Contents/MacOS/ProteomicsAnalyzer"
cat > "$APP/Contents/Info.plist" <<PL
<?xml version="1.0" encoding="UTF-8"?>
<!DOCTYPE plist PUBLIC "-//Apple//DTD PLIST 1.0//EN" "http://www.apple.com/DTDs/PropertyList-1.0.dtd">
<plist version="1.0"><dict>
  <key>CFBundleName</key><string>Proteomics Analyzer</string>
  <key>CFBundleDisplayName</key><string>Proteomics Analyzer</string>
  <key>CFBundleIdentifier</key><string>io.github.proteomics-analyzer</string>
  <key>CFBundleVersion</key><string>$VERSION</string>
  <key>CFBundleShortVersionString</key><string>$VERSION</string>
  <key>CFBundlePackageType</key><string>APPL</string>
  <key>CFBundleExecutable</key><string>ProteomicsAnalyzer</string>
  <key>CFBundleIconFile</key><string>AppIcon</string>
  <key>LSMinimumSystemVersion</key><string>12.0</string>
  <key>NSHighResolutionCapable</key><true/>
</dict></plist>
PL
PY="$(command -v python3)"
"$PY" "$HERE/macos/make_icon.py" "$APP/Contents/Resources/AppIcon.icns"
touch "$APP"
echo "Built: $APP"

#!/usr/bin/env bash
# Adds HuntingTool to the Linux application menu (current user only).
set -e
DIR="$(cd "$(dirname "$0")/.." && pwd)"
mkdir -p "$HOME/.local/share/applications"
cat > "$HOME/.local/share/applications/huntingtool.desktop" <<DESK
[Desktop Entry]
Type=Application
Name=HuntingTool
Comment=HTT CAG-repeat genotyping from amplicon long reads
Exec=bash "$DIR/launchers/HuntingTool.sh"
Terminal=false
Categories=Science;Biology;
DESK
chmod +x "$DIR/launchers/HuntingTool.sh"
echo "HuntingTool added to the application menu."

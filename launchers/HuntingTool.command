#!/usr/bin/env bash
# macOS: double-click in Finder to open the HuntingTool window.
cd "$(dirname "$0")/.." || exit 1
if command -v huntingtool-gui >/dev/null 2>&1; then exec huntingtool-gui; fi
exec python3 huntingtool.py gui

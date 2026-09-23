#!/bin/bash

# ==============================================================================
# render_memo.sh — render a Quarto .md memo to PDF, matching Rubén's style
#
# The memo YAML should set `format: html` (Quarto's Bootstrap-based HTML theme
# gives us Rubén's look: sans-serif body, section rules, clean tables,
# monospace inline code). This script renders that HTML then uses headless
# Chrome to save it as PDF, no LaTeX/MiKTeX needed.
#
# Usage:
#   bash scripts/render_memo.sh output/FBX_reply_memo.md
#
# Produces:
#   output/FBX_reply_memo.html    (embedded resources; can be shared as-is)
#   output/FBX_reply_memo.pdf
# ==============================================================================

set -euo pipefail

if [ "$#" -ne 1 ]; then
  echo "usage: $0 path/to/memo.md" >&2
  exit 1
fi

memo="$1"
if [ ! -f "$memo" ]; then
  echo "ERROR: $memo not found" >&2
  exit 1
fi

dir=$(dirname "$memo")
base=$(basename "$memo" .md)
html="$dir/$base.html"
pdf="$dir/$base.pdf"

CHROME="/c/Program Files/Google/Chrome/Application/chrome.exe"
if [ ! -f "$CHROME" ]; then
  echo "ERROR: chrome.exe not at $CHROME" >&2
  echo "Edit render_memo.sh to point at your Chrome install." >&2
  exit 1
fi

echo "Rendering $memo -> HTML..."
( cd "$dir" && quarto render "$(basename "$memo")" --to html )

echo "Converting HTML -> PDF (headless Chrome)..."
# Chrome on Windows needs Windows-style absolute paths for both --print-to-pdf
# and the input file:/// URL. pwd -W in Git Bash gives us those.
abs_dir=$(cd "$dir" && pwd -W)
abs_html="$abs_dir/$base.html"
abs_pdf="$abs_dir/$base.pdf"
"$CHROME" --headless --disable-gpu --no-pdf-header-footer \
  --print-to-pdf="$abs_pdf" "file:///$abs_html" 2>&1 | tail -2

ls -la "$html" "$pdf"

#!/usr/bin/env bash
# Convert a Markdown file to PDF using pandoc.
#
# Usage:
#   md_to_pdf.sh input.md [output.pdf]
#
# If output path is omitted, the PDF is placed alongside the input file.
# Requires: pandoc (brew install pandoc)

set -euo pipefail

if [[ $# -lt 1 ]]; then
    echo "Usage: $(basename "$0") input.md [output.pdf]" >&2
    exit 1
fi

INPUT="$1"
OUTPUT="${2:-${INPUT%.md}.pdf}"

pandoc \
    --pdf-engine=xelatex \
    -V geometry:margin=2cm \
    -V fontsize=11pt \
    --metadata title="$(basename "${INPUT%.md}")" \
    -o "$OUTPUT" \
    "$INPUT"

echo "Written: $OUTPUT"

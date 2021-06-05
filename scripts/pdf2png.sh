#!/usr/bin/env bash
# $1: pdf file

set -euo pipefail

pdf="$1"
num_pages=$(pdfinfo "$pdf" | \grep "^Pages" | awk '{print $2}')
echo "${num_pages} pages"

out_base="${pdf%.*}"

for i in $(seq 1 "$num_pages"); do
    echo "$i"
    inkscape --export-type="png" --export-filename="${out_base}${i}.png" --pdf-page="$i" --export-dpi=600 "$pdf"
done

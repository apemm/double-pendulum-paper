#!/usr/bin/env bash
# Build every PDF in this folder with a fixed clock, so that the PDF metadata
# carries no time zone (the anonymous submission must not hint at where it was
# made) and rebuilding does not change the files. Run from this folder.
#
#   bash build.sh
#
# manuscript.tex is the source. The wrappers select the variants:
#   manuscript.tex             anonymous submission (default build)
#   manuscript_named.tex       with authors and acknowledgments
#   manuscript_long.tex        anonymous, with the material cut for length
#   manuscript_twocolumn.tex   two-column, only to measure length
#   manuscript_long_twocolumn.tex
set -e
export SOURCE_DATE_EPOCH=0
export FORCE_SOURCE_DATE=1
rm -f *Notes.bib *.aux *.bbl *.blg
for f in manuscript manuscript_named manuscript_long manuscript_twocolumn manuscript_long_twocolumn; do
  pdflatex -interaction=nonstopmode "$f" > /dev/null
  bibtex "$f" > "$f.bib.buildlog" 2>&1 || true
  pdflatex -interaction=nonstopmode "$f" > /dev/null
  pdflatex -interaction=nonstopmode "$f" > "$f.buildlog"
  echo "$f: $(grep -E 'Output written' "$f.buildlog")"
  grep -E '^!|Reference.*undefined|Citation.*undefined|multiply defined' "$f.buildlog" || true
done
pdflatex -interaction=nonstopmode cover_letter > cover_letter.buildlog
echo "cover_letter: $(grep -E 'Output written' cover_letter.buildlog)"

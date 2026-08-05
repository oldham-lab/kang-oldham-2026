#!/usr/bin/env bash
# Optimize SVG files in-place using scour.
# Usage: ./optimize_svgs.sh [svg_dir] [jobs]

SVG_DIR="${1:-/home/gugene/ShinyApps/copacabana/www/dCoPA_snapshots}"
JOBS="${2:-4}"

if [[ ! -d "$SVG_DIR" ]]; then
  echo "Directory not found: $SVG_DIR" >&2
  exit 1
fi

mapfile -t SVGS < <(find "$SVG_DIR" -name "*.svg" | sort)
TOTAL=${#SVGS[@]}

if [[ $TOTAL -eq 0 ]]; then
  echo "No SVG files found in $SVG_DIR" >&2
  exit 1
fi

echo "Found $TOTAL SVG files — running scour with $JOBS parallel jobs"
echo "Target: $SVG_DIR"
echo ""

RESULTS_DIR=$(mktemp -d)

FIX_CSS_PY='
import re, sys
src, dst = sys.argv[1], sys.argv[2]
with open(src) as f:
    svg = f.read()
shape_re = re.compile(r"<(line|polyline|polygon|path|rect|circle)(\s[^>]*)?(/>|>)", re.DOTALL)
def fix_element(m):
    tag, attrs, close = m.group(1), m.group(2) or "", m.group(3)
    if not re.search(r"\bfill\s*=", attrs):
        attrs = " fill=\"none\"" + attrs
    if not re.search(r"\bstroke\s*=", attrs):
        attrs = " stroke=\"#000000\"" + attrs
    return f"<{tag}{attrs}{close}"
svg = shape_re.sub(fix_element, svg)
svg = re.sub(r"\s*fill:\s*none;\s*\n?", " ", svg)
svg = re.sub(r"\s*stroke:\s*#000000;\s*\n?", " ", svg)
with open(dst, "w") as f:
    f.write(svg)
'

optimize_one() {
  local f="$1"
  local results_dir="$2"
  local tmp tmp2
  tmp=$(mktemp --suffix=.svg)
  tmp2=$(mktemp --suffix=.svg)
  local orig_size
  orig_size=$(wc -c < "$f")

  scour \
    --enable-viewboxing \
    --enable-id-stripping \
    --enable-comment-stripping \
    --indent=none \
    --set-precision=4 \
    -i "$f" -o "$tmp" > /dev/null 2>&1

  # Fix CSS specificity: scour converts inline styles to presentation attrs,
  # but the .svglite CSS class overrides those. Explicitly add fill/stroke attrs
  # to elements relying on class defaults, then remove the conflicting CSS rules.
  python3 - "$tmp" "$tmp2" <<< "$FIX_CSS_PY" 2>/dev/null && mv "$tmp2" "$tmp"

  local new_size
  new_size=$(wc -c < "$tmp")

  if [[ $new_size -gt 0 && $new_size -lt $orig_size ]]; then
    mv "$tmp" "$f"
  else
    rm -f "$tmp"
    new_size=$orig_size
  fi

  # write result to a per-file temp file to avoid interleaved output
  echo "$orig_size $new_size" > "$results_dir/$(basename "$f").result"
}

export FIX_CSS_PY
export -f optimize_one

printf '%s\n' "${SVGS[@]}" | xargs -P "$JOBS" -I{} bash -c 'optimize_one "$@"' _ {} "$RESULTS_DIR"

echo "Done."
echo ""

ORIG_TOTAL=0
NEW_TOTAL=0
COUNT=0
while IFS=' ' read -r orig new; do
  ORIG_TOTAL=$(( ORIG_TOTAL + orig ))
  NEW_TOTAL=$(( NEW_TOTAL + new ))
  COUNT=$(( COUNT + 1 ))
done < <(cat "$RESULTS_DIR"/*.result 2>/dev/null)

rm -rf "$RESULTS_DIR"

if [[ $ORIG_TOTAL -gt 0 ]]; then
  awk -v orig="$ORIG_TOTAL" -v new="$NEW_TOTAL" -v count="$COUNT" '
    BEGIN {
      saved = orig - new
      pct   = saved * 100 / orig
      printf "  Files processed: %d\n", count
      printf "  Before: %.1f MB\n", orig / 1048576
      printf "  After:  %.1f MB\n", new  / 1048576
      printf "  Saved:  %.1f MB (%.0f%%)\n", saved / 1048576, pct
    }'
else
  echo "  No results collected."
fi

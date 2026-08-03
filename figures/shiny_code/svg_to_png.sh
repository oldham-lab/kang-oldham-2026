#!/usr/bin/env bash
# Convert dCoPA snapshot SVGs to PNGs.
# Reads svg_map CSVs to copy only app-referenced files, then converts in-place.
# Usage: ./svg_to_png.sh [src_dir] [dst_dir] [dpi] [jobs]

SRC_DIR="${1:-${REPO_DIR:-/home/gugene/code/git/kang-oldham-2026}/figures/fig_6/dCoPA_snapshots}"
DST_DIR="${2:-/home/gugene/ShinyApps/CoPA/www/dCoPA_snapshots}"
SVG_MAP_DIR="${SVG_MAP_DIR:-/home/gugene/ShinyApps/CoPA/www/dcopa/svg_map}"
DPI="${3:-144}"
JOBS="${4:-4}"

if [[ ! -d "$SRC_DIR" ]]; then
  echo "Source directory not found: $SRC_DIR" >&2; exit 1
fi
if [[ ! -d "$SVG_MAP_DIR" ]]; then
  echo "svg_map dir not found: $SVG_MAP_DIR" >&2; exit 1
fi

# Build tab-separated list of src→dst SVG paths from svg_map CSVs
COPY_LIST=$(mktemp)

python3 - "$SVG_MAP_DIR" "$SRC_DIR" "$DST_DIR" "$COPY_LIST" << 'PYEOF'
import csv, sys, os

svg_map_dir, src_dir, dst_dir, copy_list = sys.argv[1], sys.argv[2], sys.argv[3], sys.argv[4]

csv_to_subdir = {
    "ad_dfc.csv":          "AllADVsCon_DFC",
    "ad_dfc_admods.csv":   "AllADVsCon_DFC_ROSMAP",
    "ad_mtg.csv":          "AllADVsCon_MTG",
    "ad_mtg_admods.csv":   "AllADVsCon_MTG_ROSMAP",
    "scz_dfc.csv":         "SCZvsCon_DFC_bulkmegaset",
    "scz_dfc_sczmods.csv": "SCZvsCon_DFC_bulkmegaset_Brainseq",
}

total = 0
missing = 0
with open(copy_list, "w") as out:
    for csv_name, subdir in csv_to_subdir.items():
        csv_path = os.path.join(svg_map_dir, csv_name)
        with open(csv_path) as f:
            mods = sorted(set(row["mod"] for row in csv.DictReader(f)), key=int)
        for mod in mods:
            src = os.path.join(src_dir, subdir, f"{mod}.svg")
            dst = os.path.join(dst_dir, subdir, f"{mod}.svg")
            if not os.path.exists(src):
                print(f"WARNING: missing source {src}", file=sys.stderr)
                missing += 1
                continue
            out.write(f"{src}\t{dst}\n")
            total += 1

print(f"  {total} files to copy ({missing} missing in source)")
PYEOF

TOTAL=$(wc -l < "$COPY_LIST")
echo "Copying $TOTAL SVG files to $DST_DIR ..."

while IFS=$'\t' read -r src dst; do
  mkdir -p "$(dirname "$dst")"
  cp "$src" "$dst"
done < "$COPY_LIST"

echo "Done copying."
echo ""

# Convert all copied SVGs to PNG in parallel
mapfile -t SVGS < <(awk -F'\t' '{print $2}' "$COPY_LIST")
rm -f "$COPY_LIST"

echo "Converting $TOTAL SVGs to PNG at ${DPI} DPI with $JOBS parallel jobs ..."
echo ""

RESULTS_DIR=$(mktemp -d)

convert_one() {
  local f="$1"
  local dpi="$2"
  local results_dir="$3"
  local png="${f%.svg}.png"
  local svg_size
  svg_size=$(wc -c < "$f")

  rsvg-convert --dpi-x "$dpi" --dpi-y "$dpi" -o "$png" "$f" 2>/dev/null

  if [[ $? -eq 0 && -s "$png" ]]; then
    rm -f "$f"
    local png_size
    png_size=$(wc -c < "$png")
    echo "$svg_size $png_size ok" > "$results_dir/$(basename "$f").result"
  else
    rm -f "$png"
    echo "$svg_size 0 fail" > "$results_dir/$(basename "$f").result"
    echo "FAILED: $f" >&2
  fi
}

export -f convert_one

printf '%s\n' "${SVGS[@]}" | xargs -P "$JOBS" -I{} bash -c 'convert_one "$@"' _ {} "$DPI" "$RESULTS_DIR"

echo "Done."
echo ""

SVG_TOTAL=0
PNG_TOTAL=0
OK=0
FAIL=0
while IFS=' ' read -r svg_size png_size status; do
  SVG_TOTAL=$(( SVG_TOTAL + svg_size ))
  PNG_TOTAL=$(( PNG_TOTAL + png_size ))
  if [[ $status == ok ]]; then (( OK++ )); else (( FAIL++ )); fi
done < <(cat "$RESULTS_DIR"/*.result 2>/dev/null)

rm -rf "$RESULTS_DIR"

awk -v svg="$SVG_TOTAL" -v png="$PNG_TOTAL" -v ok="$OK" -v fail="$FAIL" '
  BEGIN {
    printf "  Converted: %d  Failed: %d\n", ok, fail
    printf "  SVG total: %.1f MB\n", svg / 1048576
    printf "  PNG total: %.1f MB\n", png / 1048576
    printf "  Saved:     %.1f MB (%.0f%%)\n", (svg-png)/1048576, (svg-png)*100/svg
  }'

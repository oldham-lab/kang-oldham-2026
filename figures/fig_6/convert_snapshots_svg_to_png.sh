#!/usr/bin/env bash
#
# Recreate the dCoPA snapshot PNGs from their source SVGs.
#
# For every PNG that already exists under TARGET_DIR, the matching SVG is
# located under SOURCE_DIR (same relative path, .svg instead of .png) and
# re-rendered, overwriting the PNG in place. No new files or subdirectories
# are created — this reproduces the existing 471-file PNG set exactly.
#
# Rendering: 144 DPI, which matches the originals (2x the SVG's point size,
# e.g. 432x864 pt -> 864x1728 px). This scales correctly for every plot
# regardless of its facet count / aspect ratio.
#
# Requires: rsvg-convert (librsvg).

set -euo pipefail

SOURCE_DIR="${REPO_DIR:-/home/gugene/code/git/kang-oldham-2026}/figures/fig_6/dCoPA_snapshots"
TARGET_DIR="${SHINYAPP_DIR:-/home/gugene/ShinyApps/copacabana}/www/dCoPA_snapshots"
DPI=144

command -v rsvg-convert >/dev/null 2>&1 || {
    echo "ERROR: rsvg-convert not found (install librsvg2-bin)." >&2
    exit 1
}

converted=0
missing=0

# Drive off the existing PNGs so only the current set is recreated.
while IFS= read -r -d '' png; do
    rel="${png#"$TARGET_DIR"/}"          # e.g. AllADVsCon_DFC/11.png
    svg="$SOURCE_DIR/${rel%.png}.svg"    # matching source SVG

    if [[ -f "$svg" ]]; then
        rsvg-convert -d "$DPI" -p "$DPI" "$svg" -o "$png"
        converted=$((converted + 1))
    else
        echo "WARNING: no source SVG for $rel (expected $svg)" >&2
        missing=$((missing + 1))
    fi
done < <(find "$TARGET_DIR" -type f -name '*.png' -print0)

echo "Done. Converted: $converted   Missing source SVG: $missing"

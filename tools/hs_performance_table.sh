#!/usr/bin/env bash
#
# Regenerate the Performance table of the user guide (web/index.html,
# section "Performance"): runs the Hock-Schittkowski test suite
# (test/test_hs_suite.f90, 305 problems) once for each globalization
# configuration below (line search / merit function / penalty, or trust region), and prints the table rows (HTML by
# default, or Markdown with --markdown), ready to paste into the guide.
#
# usage (from the repository root):
#
#    pixi run tools/hs_performance_table.sh [--markdown]
#
# Each run takes about a second (release build). The per-run Markdown
# reports are left in a temporary directory, printed at the end.

set -euo pipefail

FPM=${FPM:-fpm}
format=html
if [[ "${1:-}" == "--markdown" ]]; then
    format=markdown
fi

# label (HTML) | label (Markdown) | harness options  (the first row is the default)
rows=(
  "<strong>filter</strong> (default, with interpolation)|**filter** (default, with interpolation)|"
  "filter, without interpolation|filter, without interpolation|--no-interpolate"
  "funnel (with interpolation)|funnel (with interpolation)|--linesearch=funnel"
  "Armijo / \\( \\ell_1 \\) / multipliers|Armijo / ℓ1 / multipliers|--linesearch=armijo --no-interpolate"
  "Armijo / \\( \\ell_1 \\) / model (Byrd&ndash;Nocedal)|Armijo / ℓ1 / model (Byrd-Nocedal)|--linesearch=armijo --penalty=model --no-interpolate"
  "Armijo / augmented Lagrangian / multipliers|Armijo / augmented Lagrangian / multipliers|--linesearch=armijo --merit=al --no-interpolate"
  "Armijo / augmented Lagrangian / multipliers, interpolation + non-monotone (10)|Armijo / augmented Lagrangian / multipliers, interpolation + non-monotone (10)|--linesearch=armijo --merit=al --nonmonotone=10"
  "Armijo / augmented Lagrangian / model (GMSW)|Armijo / augmented Lagrangian / model (GMSW)|--linesearch=armijo --merit=al --penalty=model --no-interpolate"
  "watchdog / \\( \\ell_1 \\) / model|watchdog / ℓ1 / model|--linesearch=watchdog --penalty=model --no-interpolate"
  "trust region / filter|trust region / filter|--trust-region"
  "trust region / funnel|trust region / funnel|--trust-region --linesearch=funnel"
)
# (and the footnote's "Armijo / l1 / multipliers, with interpolation" figures)
extra="--linesearch=armijo"

outdir=$(mktemp -d)

# thousands separators (1234567 -> 1,234,567)
commas() { awk -v n="$1" 'BEGIN { s = ""; while (length(n) > 3) { s = "," substr(n, length(n)-2) s; n = substr(n, 1, length(n)-3) }; print n s }'; }

# run the suite with the given options; prints "solved local failed nf"
run() {
    local report=$1; shift
    local line
    line=$($FPM test test_hs_suite --profile release -- "$report" "$@" 2>&1 | grep -a '^summary:' || true)
    if [[ -z "$line" ]]; then
        echo "error: no summary from test_hs_suite $*" >&2
        exit 1
    fi
    echo "$line" | sed -E 's/summary: solved=([0-9]+) local=([0-9]+) failed=([0-9]+) nf=([0-9]+) ng=([0-9]+)/\1 \2 \3 \4/'
}

# build first, so that the build output doesn't mix with the table
if ! build_log=$($FPM build --profile release --tests 2>&1); then
    echo "$build_log" >&2
    exit 1
fi

if [[ $format == markdown ]]; then
    echo "| globalization | solved | local | failed | \`fc\` calls (solved) |"
    echo "|---|--:|--:|--:|--:|"
fi

i=0
for row in "${rows[@]}"; do
    IFS='|' read -r html_label md_label opts <<< "$row"
    # shellcheck disable=SC2086
    read -r solved local failed nf <<< "$(run "$outdir/row$i.md" $opts)"
    nf=$(commas "$nf")
    if [[ $format == markdown ]]; then
        if [[ $i == 0 ]]; then
            echo "| $md_label | **$solved** | $local | **$failed** | **$nf** |"
        else
            echo "| $md_label | $solved | $local | $failed | $nf |"
        fi
    else
        if [[ $i == 0 ]]; then
            echo "<tr class=\"best\"><td>$html_label</td><td class=\"num\"><strong>$solved</strong></td><td class=\"num\">$local</td><td class=\"num\"><strong>$failed</strong></td><td class=\"num\"><strong>$nf</strong></td></tr>"
        else
            echo "<tr><td>$html_label</td><td class=\"num\">$solved</td><td class=\"num\">$local</td><td class=\"num\">$failed</td><td class=\"num\">$nf</td></tr>"
        fi
    fi
    i=$((i+1))
done

# shellcheck disable=SC2086
read -r solved local failed nf <<< "$(run "$outdir/extra.md" $extra)"
echo
echo "footnote: Armijo / l1 / multipliers with interpolation: solved=$solved local=$local failed=$failed nf=$(commas "$nf")"
echo "reports: $outdir"

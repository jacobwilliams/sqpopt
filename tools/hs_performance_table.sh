#!/usr/bin/env bash
#
# Regenerate the Hock-Schittkowski table of the Performance page
# (web/performance.html, "By configuration"): runs the Hock-Schittkowski test suite
# (test/test_hs_suite.f90, 305 problems) once for each globalization
# configuration below (line search / merit function / penalty, or trust region), and prints the table rows (HTML by
# default, or Markdown with --markdown), ready to paste into the guide. The
# default configuration's run also regenerates the data of the interactive
# results page (web/js/hs_results_data.js, for the same page's results of every problem), and the
# SLSQP comparison's data is regenerated too (test/test_hs_slsqp.f90,
# web/js/hs_slsqp_data.js).
#
# usage (from the repository root):
#
#    pixi run tools/hs_performance_table.sh [--markdown] [--mumps]
#
# The rows of the options that factor matrices (options%inertia_control and
# options%direct_qp) use the default sparse solver, QDLDL. With --mumps, the
# library is built with MUMPS too (the HAS_MUMPS preprocessor directive, which
# needs the sequential MUMPS library of the pixi environment), and the table
# gets the same rows with options%linear_solver = MUMPS, for comparison. The
# other rows don't depend on the sparse solver.
#
# Each run takes about a second (release build). The per-run Markdown
# reports are left in a temporary directory, printed at the end.

set -euo pipefail

FPM=${FPM:-fpm}
format=html
mumps=0
for arg in "$@"; do
    case "$arg" in
        --markdown) format=markdown ;;
        --mumps)    mumps=1 ;;
        *) echo "usage: $0 [--markdown] [--mumps]" >&2; exit 2 ;;
    esac
done

# extra fpm options, for the build with MUMPS
fpm_flags=()
if [[ $mumps == 1 ]]; then
    fpm_flags=(--flag "-DHAS_MUMPS -I${CONDA_PREFIX:?run this in the pixi environment}/include" --link-flag "-ldmumps_seq")
fi

# label (HTML) | label (Markdown) | harness options  (the first row is the default)
rows=(
  "<strong>filter</strong> (default, with interpolation)|**filter** (default, with interpolation)|--web-data=web/js/hs_results_data.js"
  "filter, without interpolation|filter, without interpolation|--no-interpolate"
  "funnel (with interpolation)|funnel (with interpolation)|--linesearch=funnel"
  "filter, exact Hessian (finite differences of the gradients)|filter, exact Hessian (finite differences of the gradients)|--hessian=exact"
  "Armijo / \\( \\ell_1 \\) / multipliers|Armijo / ℓ1 / multipliers|--linesearch=armijo --no-interpolate"
  "Armijo / \\( \\ell_1 \\) / model (Byrd&ndash;Nocedal)|Armijo / ℓ1 / model (Byrd-Nocedal)|--linesearch=armijo --penalty=model --no-interpolate"
  "Armijo / augmented Lagrangian / multipliers|Armijo / augmented Lagrangian / multipliers|--linesearch=armijo --merit=al --no-interpolate"
  "Armijo / augmented Lagrangian / multipliers, interpolation + non-monotone (10)|Armijo / augmented Lagrangian / multipliers, interpolation + non-monotone (10)|--linesearch=armijo --merit=al --nonmonotone=10"
  "Armijo / augmented Lagrangian / model (GMSW)|Armijo / augmented Lagrangian / model (GMSW)|--linesearch=armijo --merit=al --penalty=model --no-interpolate"
  "watchdog / \\( \\ell_1 \\) / model|watchdog / ℓ1 / model|--linesearch=watchdog --penalty=model --no-interpolate"
  "trust region / filter|trust region / filter|--trust-region"
  "trust region / funnel|trust region / funnel|--trust-region --linesearch=funnel"
)
# the rows of the options that factor matrices (after the exact-Hessian row),
# with QDLDL, and with MUMPS too if it is built in
factored=(
  "exact Hessian with inertia control|--hessian=exact --inertia"
  "exact Hessian with inertia control and direct QP|--hessian=exact --inertia --direct"
  "L-SR1 with inertia control|--hessian=sr1 --inertia"
  "L-BFGS with direct QP|--direct"
)
with_factored=()
for row in "${rows[@]}"; do
    with_factored+=("$row")
    if [[ "$row" == *"|--hessian=exact" ]]; then
        for f in "${factored[@]}"; do
            IFS='|' read -r label opts <<< "$f"
            with_factored+=("filter, $label|filter, $label|$opts")
            if [[ $mumps == 1 ]]; then
                with_factored+=("filter, $label (MUMPS)|filter, $label (MUMPS)|$opts --linear-solver=mumps")
            fi
        done
    fi
done
rows=("${with_factored[@]}")
# (and the footnote's "Armijo / l1 / multipliers, with interpolation" figures)
extra="--linesearch=armijo"

outdir=$(mktemp -d)

# thousands separators (1234567 -> 1,234,567)
commas() { awk -v n="$1" 'BEGIN { s = ""; while (length(n) > 3) { s = "," substr(n, length(n)-2) s; n = substr(n, 1, length(n)-3) }; print n s }'; }

# run the suite with the given options; prints "solved local failed nf"
run() {
    local report=$1; shift
    local line
    line=$($FPM test test_hs_suite --profile release ${fpm_flags[@]+"${fpm_flags[@]}"} -- "$report" "$@" 2>&1 \
           | grep -a '^summary:' || true)
    if [[ -z "$line" ]]; then
        echo "error: no summary from test_hs_suite $*" >&2
        exit 1
    fi
    echo "$line" | sed -E 's/summary: solved=([0-9]+) local=([0-9]+) failed=([0-9]+) nf=([0-9]+) ng=([0-9]+)/\1 \2 \3 \4/'
}

# build first, so that the build output doesn't mix with the table
if ! build_log=$($FPM build --profile release --tests ${fpm_flags[@]+"${fpm_flags[@]}"} 2>&1); then
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

# the SLSQP comparison of the results page
slsqp=$($FPM test test_hs_slsqp --profile release ${fpm_flags[@]+"${fpm_flags[@]}"} -- --web-data=web/js/hs_slsqp_data.js 2>&1 \
        | grep -a '^summary:' || true)
if [[ -z "$slsqp" ]]; then
    echo "error: no summary from test_hs_slsqp" >&2
    exit 1
fi
echo "SLSQP (results page data): ${slsqp#summary: }"
echo "reports: $outdir"

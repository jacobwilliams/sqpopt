#!/usr/bin/env bash
#
# Run one Hock-Schittkowski problem with both SQPOPT and SLSQP, printing each
# solver's iterations, for debugging a difference between them (e.g. on the
# interactive results page, web/hs_results.html).
#
# usage (from the repository root):
#
#    pixi run tools/hs_compare.sh N [SQPOPT options...]
#
# N is the problem number (e.g. 220). Any further options are passed to the
# SQPOPT harness (test/test_hs_suite.f90; e.g. --linesearch=funnel,
# --hessian=exact, --print=1); its iteration log is printed with --print=2
# unless another --print=L is given. SLSQP is run by test/test_hs_slsqp.f90
# with its fixed settings.

set -euo pipefail

FPM=${FPM:-fpm}
if [[ $# -lt 1 || ! "$1" =~ ^[0-9]+$ ]]; then
    echo "usage: $0 N [SQPOPT options...]" >&2
    exit 1
fi
id=$1; shift
opts=("$@")
if ! printf '%s\n' "${opts[@]+"${opts[@]}"}" | grep -q '^--print='; then
    opts+=(--print=2)
fi

# build first, so that the build output doesn't mix with the logs
if ! build_log=$($FPM build --profile release --tests 2>&1); then
    echo "$build_log" >&2
    exit 1
fi

echo "==================== SQPOPT: TP$id ===================="
$FPM test test_hs_suite --profile release -- --problem="$id" "${opts[@]}" 2>&1 \
    | grep -av -e '^ *-\{5,\}' -e '^ test_hs_suite$' -e '^test_hs_suite: non-default' \
                -e '^problems:' -e '^solved:' -e '^local solutions:' -e '^failed:' -e '^with FD' \
                -e '^evaluations (solved' -e '^ *NLPQLP' -e '^time:' -e '^summary:'
echo
echo "==================== SLSQP: TP$id ===================="
$FPM test test_hs_slsqp --profile release -- --problem="$id" --print 2>&1 \
    | grep -av -e '^ *-\{5,\}' -e '^ test_hs_slsqp$' -e '^test_hs_slsqp PASSED' -e '^summary:'

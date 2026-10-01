#!/bin/bash
# Script to run tests with coverage and generate reports
#
# usage:  ./coverage.sh [--mumps]

set -e

COV_DIR="coverage"

# Find the correct gcov (matching gfortran version)
GCOV=$(dirname $(which gfortran))/gcov

echo "Using gcov: $GCOV"
echo ""
echo "================================================"
echo "Running tests with coverage..."
echo "================================================"

# Clear the execution counts of any previous run (so they don't accumulate)
find build -name "*.gcda" -delete 2>/dev/null || true

# Build and run tests with coverage flags. With --mumps, the library is built
# with MUMPS (the HAS_MUMPS preprocessor directive, see the README), so that
# the report covers the code that uses it (in a pixi environment, which has
# the sequential MUMPS library).
if [[ "${1:-}" == "--mumps" ]]; then
    fpm test --flag "--coverage -DHAS_MUMPS -I${CONDA_PREFIX:?run this in the pixi environment}/include" --link-flag "-ldmumps_seq"
else
    fpm test --flag "--coverage"
fi

echo ""
echo "================================================"
echo "Generating coverage report with lcov..."
echo "================================================"

# Create a fresh coverage directory and copy the build artifacts into it. Only
# from the build just made: older coverage builds (with a different build-flag
# hash, e.g. after an fpm.toml change) hold stale data for the same files.
LATEST_GCNO=$(ls -t build/gfortran_*/*/*.gcno | head -1)
BUILD_DIR=$(dirname "$(dirname "$LATEST_GCNO")")
echo "Using build directory: $BUILD_DIR"
rm -rf "$COV_DIR"
mkdir -p "$COV_DIR"
cp -r "$BUILD_DIR"/*/* "$COV_DIR/" 2>/dev/null || true
# don't show coverage of test/examples programs:
rm -rf "$COV_DIR"/test_*
rm -rf "$COV_DIR"/example_*
# we want to skip the files that start with 'build_dependencies' because they are not part of the source code
find "$COV_DIR" -name "build_dependencies*" -exec rm -rf {} +

# Capture initial coverage (baseline)
lcov --ignore-errors inconsistent --gcov-tool "$GCOV" --capture --initial --base-directory . --directory "$COV_DIR" --output-file "$COV_DIR/coverage.base" --quiet

# Capture test execution coverage
lcov --ignore-errors inconsistent --gcov-tool "$GCOV" --capture --base-directory . --directory "$COV_DIR" --output-file "$COV_DIR/coverage.capture" --quiet

# Combine baseline and test coverage
lcov --add-tracefile "$COV_DIR/coverage.base" --add-tracefile "$COV_DIR/coverage.capture" --output-file "$COV_DIR/coverage.info" --quiet

echo ""
echo "================================================"
echo "Coverage Summary"
echo "================================================"

# Display summary
lcov --summary "$COV_DIR/coverage.info"

echo ""
echo "Coverage data saved to: $COV_DIR/coverage.info"
echo ""
echo "generate HTML report:"

genhtml $COV_DIR/coverage.info --output-directory $COV_DIR/html --dark-mode

echo ""

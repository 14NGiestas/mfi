#!/usr/bin/env bash

set -ex

# Target directory to deploy mfi to
destdir="${DESTDIR:-mfi-fpm}"

# Get fypp preprocessor
fypp="${FYPP:-$(which fypp)}"

# Arguments for the fypp preprocessor
fyflags="${FYFLAGS:--DMAXRANK=4}"

# Number of parallel jobs for preprocessing
njob="$(nproc)"

# Additional files to include
include=(
  "fpm.toml"
  "LICENSE"
)

# Files to remove from collection
prune=(
)

mkdir -p "$destdir/src" "$destdir/test"

# Preprocess mfi sources
find src -maxdepth 1 -iname "*.fypp" \
  | cut -f1 -d. | xargs -P "$njob" -I{} "$fypp" "{}.fypp" "$destdir/{}.f90" $fyflags

# Collect mfi source files
find src -maxdepth 1 -iname "*.f90" -exec cp {} "$destdir/src/" \;
find src/tests -name "test_*.f90" -exec cp {} "$destdir/test/" \;
find src/tests -name "*.dat" -exec cp {} "$destdir/" \;

# Include additional files
cp "${include[@]}" "$destdir/"

# Source file workarounds for fpm
rm "${prune[@]}"

# List mfi-fpm package contents
ls -R "$destdir"

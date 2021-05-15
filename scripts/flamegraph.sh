#!/usr/bin/env bash
# Runs the given implementation on the given scene through perf and generates a flame graph. 
# This requires perf as well as the flame graph scripts in the PATH: 
# https://github.com/brendangregg/FlameGraph

set -euo pipefail

if [ $# -lt 2 ]; then
    echo "Usage :  $0 <impl> <scene> [<width> [<height>]]"
    exit 1
fi

this=$(realpath "$(dirname "$0")")

build_dir="$this/../build"

impl="$1"
scene="$2"

if [ $# -ge 3 ]; then
    width="$3"
else
    width=1920
fi

if [ $# -ge 4 ]; then
    height="$4"
else
    height=1080
fi

scene_name="$(basename "$scene" ".json")"
name="${impl}-${scene_name}-${width}x${height}"
data="$name.data"
svg="$name.svg"

"$this/build.sh" Release "-gdwarf"

perf record -F 200 --call-graph dwarf -o "$data" "$build_dir/main" "$impl" "$scene" "$width" "$height"
perf script --input="$data" | stackcollapse-perf.pl | flamegraph.pl --title "Flame Graph: $impl - $scene_name ${width}x${height}" > "$svg"

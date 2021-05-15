#!/usr/bin/env bash

set -euo pipefail

script_dir=$(realpath "$(dirname "$0")")
dat_dir="$script_dir/../plots"

# TODO 
impls=("ref" "opt0" "opt1" "opt3")

bench_types=("size" "all" "box" "sphere" "cone" "torus" "octahedron")

for impl in "${impls[@]}"; do
    fname_base="$dat_dir/bench-${impl}"

    for bench_type in "${bench_types[@]}"; do
        echo -e "\e[32;1mRunning benchmark type '${bench_type}' for implementation '${impl}'\e[0m" >&2
        "$script_dir/bench.py" "$impl" "$bench_type" | tee "${fname_base}-${bench_type}.dat"
    done
done

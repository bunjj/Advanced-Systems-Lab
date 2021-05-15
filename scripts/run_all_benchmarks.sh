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

        case "${bench_type}" in
            size)
                x_base=200
                x_end=4200
                x_step=200
                ;;
            all|box|sphere|cone|torus|octahedron)
                x_base=10
                x_end=80
                x_step=10
                ;;
            *)
                echo "Unknown bench_type '${bench_type}'" >&2
                exit 1
                ;;
        esac

        echo -e "\e[32;1mRunning benchmark type '${bench_type}' for implementation '${impl}'\e[0m" >&2
        "$script_dir/bench.py" "$impl" "$bench_type" "$x_base" "$x_end" "$x_step" | tee "${fname_base}-${bench_type}.dat"
    done
done

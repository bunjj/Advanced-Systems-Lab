#!/usr/bin/env bash
# Generate all plots for the given implementation
# $1 Measurement

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage: $0 <measurement>"
    echo -e "\tmeasurement\t Either flops, seconds, cycles, or perf"
    exit 1
fi

names=("Gflops" "seconds" "cycles" "flops/cycle")
indices=(2 3 4 5)

declare -A plot_types=([flops]=0 [seconds]=1 [cycles]=2 [perf]=3)

measurement="$1"
idx=${plot_types[$measurement]}

echo -e "\e[32;1m======== Generating Plots for ${names[$idx]} ========\e[0m"

bench_types=("all" "box" "cone" "octahedron" "size" "sphere" "torus")

base_dir=$(realpath "$(dirname "$0")")
plot_dir="${base_dir}/../plots"

cd "$plot_dir"

for bench in "${bench_types[@]}"; do
    basename="${measurement}-${bench}"
    echo "$basename"
    gnuplot -c size-time.plt "$bench" "${indices[$idx]}" "${names[$idx]}" "$basename"
done

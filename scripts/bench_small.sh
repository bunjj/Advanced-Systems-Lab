#!/usr/bin/env bash
# Runs all local benchmarks with the given implementation
# Uses the currently built executable
# $1 implementation

set -euo pipefail

if [ $# -lt 1 ]; then
    echo "Usage :  $0 <impl>"
    exit 1
fi

impl="$1"

base_dir=$(realpath "$(dirname "$0")")
build_dir="${base_dir}/../build"
scene_dir="${base_dir}/../scenes"

scene0="$scene_dir/scene0.json"
benchmark_small="$(find "${scene_dir}/benchmark_small" -name "*.json")"

all_scenes=$(echo "${scene0}"$'\n'"${benchmark_small}" | sort)

echo "$all_scenes" | while read -r scene; do
    echo -e "\e[32;1m======== $(basename "$scene") ========\e[0m"
    "${build_dir}/main" "$impl" "${scene}" 1920 1080
done

#!/usr/bin/env bash
# Runs all local benchmarks with the given implementation and validates it
# against an existing image for each scene.  Each scene must have a
# corresponding .ppm file in the same folder for validation
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
finish() { 
    if [ -e "$tmp" ]; then
        rm "$tmp"
    fi
}

trap finish EXIT

tmp=$(mktemp --suffix=.ppm)

echo "$all_scenes" | while read -r scene; do
    echo -e "\e[32;1m======== $(basename "$scene") ========\e[0m"
    if "${build_dir}/main" "$impl" "${scene}" 1920 1080 "$tmp" "${scene%.*}.ppm"; then
        echo -e "\e[32;1mPASSED\e[0m"
    else
        echo -e "\e[31;1mFAILED\e[0m"
    fi
done

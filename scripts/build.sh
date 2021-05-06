#!/usr/bin/env bash
# $1 Build Type
# $2 Compiler Flags

set -euo pipefail

if [ $# -lt 2 ] || [ -z "$1" ]; then
    echo "Usage :  $0 <BUILD_TYPE> <CXX_FLAGS> [<cmake flags>...]"
    exit 1
fi

build_dir=$(realpath "$(dirname "$0")/../build")

if [ -d "$build_dir" ]; then
    rm -Rf "$build_dir"
fi

mkdir "$build_dir"
cd "$build_dir"

build_type="$1"
cflags="$2"

cmake -DCMAKE_BUILD_TYPE="$build_type" -DCMAKE_CXX_FLAGS="$cflags" "${@:3}" ..
make -j

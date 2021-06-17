#!/usr/bin/env python3

import bench
from bench import script_dir
from bench import Config
from bench import get_range
from pathlib import Path
import sys
import json

plots_dir = script_dir / ".." / "plots"

out_all = plots_dir / "bench.json"
out_compiler_flags = plots_dir / "bench-compiler-flags.json"

impls = [
    "ref",
    "opt0",
    "opt1",
    "opt2",
    "opt3",
    "opt4",
    "opt5",
    "vec0",
    "vec1",
    "vec2",
    "vec3",
    "vec4"]

shapes = ["all", "box", "sphere", "cone", "torus", "octahedron"]
bench_types = ["size"] + shapes
flags = [
    "-O2",
    "-O3 -fno-tree-vectorize",
    "-O3",
    "-Ofast",
    "-Ofast -march=native"]

compilers = ["g++", "clang++"]

if __name__ == "__main__":

    assert(len(sys.argv) == 2)
    d = Path(sys.argv[1])
    if not d.exists():
        d.mkdir(parents = True)

    assert(d.is_dir())

    conf_all = Config(["g++"], ["-Ofast"], impls, bench_types, shapes, get_range)
    conf_compiler_flags = Config(compilers, flags, ["vec4"], bench_types, shapes, get_range)

    data_all = bench.main(conf_all, d)
    data_compiler_flags = bench.main(conf_compiler_flags, d)

    with open(str(out_all), "w") as outfile:
        json.dump(data_all, outfile, indent=4)

    with open(str(out_compiler_flags), "w") as outfile:
        json.dump(data_compiler_flags, outfile, indent=4)



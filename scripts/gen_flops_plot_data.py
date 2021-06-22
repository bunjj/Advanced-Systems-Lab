#!/usr/bin/env python3

import sys
import json
from itertools import groupby

from pathlib import Path
plot_dir = Path(__file__).parent / ".." / "plots"


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def print_header(name):
    eprint(f"\033[32;1m{name}\033[0m")


def print_line(file, line, quotes=False):

    if quotes:
        line = [f'"{l}"' for l in line]

    print("\t".join(line), file=file)


impls = ["ref", "opt1", "opt3", "opt5", "vec4"]
impl_names = ["Base", "Inlining", "Term", "MVMs", "Vect"]
bench_types = []
shapes = []

flop_types = [
    ("ADD", 1),
    ("MUL", 1),
    ("FMA", 2),
    ("DIV", 1),
    ("SQRT", 1),
    ("ABS", 1),
    ("CMP", 1),
    ("MAX", 1),
    ("POW", 30),
    ("TAN", 30),
]


def collect_categories(datapoints: list):
    global bench_types
    global shapes

    for dp in datapoints:
        if not dp["has_flops"]:
            continue

        bench_type = dp["type"]

        bench_types.append(bench_type)

        if bench_type != "size":
            shapes.append(bench_type)

    bench_types = sorted(list(set(bench_types)))
    shapes = sorted(list(set(shapes)))

    eprint(f"Loaded {len(datapoints)} datapoints")
    eprint(f"Loaded benchmark types: {', '.join(bench_types)}")


def gen_data_file(path: Path, filtered, bench_type):
    print_header(str(path.name))

    if bench_type == "size":
        x_field = "width"
    else:
        x_field = "num_shapes"

    max_x = max(filtered, key=lambda v: v[x_field])[x_field]

    filtered = [v for v in filtered if v[x_field] == max_x]

    # We should only have one entry per implementation
    assert(len(filtered) == len(impls))

    header = ["impl", "flops"] + list(map(lambda v: v[0], flop_types))
    lines = []

    file = open(path, 'w')

    for i, impl in enumerate(impls):
        d = [v for v in filtered if v["impl"] == impl]
        assert(len(d) == 1)

        line = [impl_names[i], d[0]["Flops"]]

        for t in flop_types:
            line.append(str(int(d[0][t[0]]) * t[1]))

        lines.append(line)

    print_line(file, header, True)

    for line in lines:
        print_line(file, line, False)


def main():
    if len(sys.argv) < 2:
        eprint(f"Usage: {sys.argv[0]} <json>")
        sys.exit(1)

    input_file = sys.argv[1]

    if input_file == "-":
        json_data = "\n".join(sys.stdin.read())
        datapoints = json.loads(json_data)
    else:
        datapoints = json.load(open(input_file))

    collect_categories(datapoints)

    datapoints = [v for v in datapoints if v["has_flops"] and v["impl"] in impls]

    for bench_type in bench_types:
        for bench_type in bench_types:
            data_file = plot_dir / f"flops-{bench_type}.dat"
            filtered = [v for v in datapoints if v["type"] == bench_type]
            gen_data_file(data_file, filtered, bench_type)


if __name__ == "__main__":
    main()

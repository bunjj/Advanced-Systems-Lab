#!/usr/bin/env python3

import sys
import json

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

impls = []
flags = []
bench_types = []
shapes = []


def collect_categories(datapoints: list):
    global impls
    global flags
    global bench_types
    global shapes

    for dp in datapoints:
        impls.append(dp["impl"])

        if not dp["has_flops"]:
            flags.append(dp["flags"])

        bench_type = dp["type"]

        bench_types.append(bench_type)

        if bench_type != "size":
            shapes.append(bench_type)

    impls = sorted(list(set(impls)))
    flags = sorted(list(set(flags)))
    bench_types = sorted(list(set(bench_types)))
    shapes = sorted(list(set(shapes)))

    eprint(f"Loaded {len(datapoints)} datapoints")
    eprint(f"Loaded implementations: {', '.join(impls)}")
    eprint(f"Loaded flags: {', '.join(flags)}")
    eprint(f"Loaded benchmark types: {', '.join(bench_types)}")


def gen_data_file(path : Path, filtered, bench_type):
    print_header(str(path.name))

    if bench_type == "size":
        x_field = "width"
        x_name = f"width"
    else:
        x_field = "num_shapes"
        x_name = f"number of shapes ({bench_type})"

    flops = {v[x_field]: v["Flops"] for v in filtered if v["has_flops"]}

    # Filter and sort so that only non-instrumented measurements are included,
    # all entries with the same x value appear next to each other and are
    # sorted by flags
    filtered = [v for v in filtered if not v["has_flops"]]
    filtered = sorted(filtered, key=lambda v: (v[x_field], v["flags"]))

    header = [x_name, "flops"]
    lines = []

    for flag in flags:
        header.extend(
            [f"{flag} - seconds", f"{flag} - cycles", f"{flag} - flops/cycle"])

    file = open(path, 'w')

    i = 0
    x_val_old = -1

    for dp in filtered:
        x_val = dp[x_field]

        flop = float(flops[x_val])

        if x_val != x_val_old:
            i = 0
            lines.append([])
            lines[-1].extend([str(x_val), str(flop)])

        flag = flags[i]
        seconds = float(dp["Microseconds"]) / 1e6
        cycles = float(dp["Cycles"])

        assert(dp["flags"] == flag)

        lines[-1].extend([str(seconds), str(cycles), str(flop / cycles)])

        i += 1

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

    for impl in impls:
        for bench_type in bench_types:
            data_file = plot_dir / f"{impl}-{bench_type}.dat"
            filtered = [
                v for v in datapoints if (
                    v["impl"] == impl and v["type"] == bench_type)]
            gen_data_file(data_file, filtered, bench_type)


if __name__ == "__main__":
    main()

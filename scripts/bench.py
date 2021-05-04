#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys

script_dir = Path(__file__).parent
build_exec = script_dir / "build.sh"
run_exec = (script_dir / ".." / "build" / "main").resolve()
input_scene = script_dir / ".." / "scenes" / "scene0.json"

fields = ["Width", "Height", "Flops", "Cycles", "Microseconds"]

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def run_benchmark():
    results = []

    widths = range(100, 2100, 100)

    for width in widths:
        eprint("{}x{}".format(width, width))
        p = subprocess.run([str(run_exec), str(input_scene), str(width), str(width)], check=True, capture_output=True)

        values = {}

        for line in p.stderr.decode().split('\n'):
            if ': ' in line:
                key, value = line.split(': ')

                if key in fields:
                    # eprint("{}={}".format(key, value))
                    values[key] = value
        results.append(values)

    return results

def build(flag, instrument):
    cmd = [str(build_exec), "Release", flag]

    if instrument:
        cmd.append("-DINSTRUMENT=ON")
    else:
        cmd.append("-DINSTRUMENT=OFF")

    eprint(cmd)

    p = subprocess.run(cmd, capture_output=True)

    if p.returncode != 0:
        eprint(p.stdout.decode())
        eprint(p.stderr.decode())
        raise subprocess.CalledProcessError(p.returncode, p.args, p.stdout, p.stderr)

def print_line(line, quotes=False):
    for l in line:
        if quotes:
            print('\t"{}"'.format(l), end='')
        else:
            print('\t{}'.format(l), end='')

    print('')

def main():
    flags = ["-O3 -fno-tree-vectorize", "-O3", "-Ofast"]

    eprint("Flops Counter")
    # Count flops first
    build("", True)

    flops = run_benchmark()

    for v in flops:
        assert("Width" in v)
        assert("Height" in v)
        assert("Flops" in v)

    results = []

    for flag in flags:
        eprint(flag)

        build(flag, False)
        values = run_benchmark()

        for v in values:
            assert("Width" in v)
            assert("Height" in v)
            assert("Cycles" in v)
            assert("Microseconds" in v)

        assert(len(values) == len(flops))

        results.append(values)

    header = ["pixels", "flops"]

    for flag in flags:
        header.extend([flag + " - seconds", flag + " - cycles", flag + " - flops/cycle"])

    print_line(header, True)

    for i in range(len(flops)):
        flop_val = flops[i]
        width = float(flop_val["Width"])
        height = float(flop_val["Height"])
        num_flops = float(flop_val["Flops"])

        line = [width * height, num_flops]

        for result in results:
            r = result[i]
            usec = float(r["Microseconds"])
            cycles = float(r["Cycles"])
            line.extend([usec / 1e6, cycles, num_flops / cycles])

        print_line(line, False)

if __name__=="__main__":
    main()   

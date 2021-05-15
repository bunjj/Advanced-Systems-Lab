#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys
from collections import OrderedDict
import tempfile

script_dir = Path(__file__).parent
build_exec = script_dir / "build.sh"
run_exec = (script_dir / ".." / "build" / "main").resolve()
input_scene = script_dir / ".." / "scenes" / "scene0.json"
shape_scene_creator = script_dir / ".." / "benchmark" / "create_scene.py"
all_shapes_scene_creator = script_dir / ".." / \
    "benchmark" / "create_scene_all_shapes.py"

fields = ["Width", "Height", "Flops", "Cycles", "Microseconds"]

x_base = 0
x_end = 0
x_step = 0


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)

def get_range():
    return range(x_base, x_end, x_step)

def run_single(impl, scene, width, height):
    p = subprocess.run([str(run_exec), impl, str(scene), str(
        width), str(height)], check=True, capture_output=True)

    values = {}

    for line in p.stderr.decode().split('\n'):
        if ': ' in line:
            key, value = line.split(': ')

            if key in fields:
                # eprint("{}={}".format(key, value))
                values[key] = value

    return values


def run_size_benchmark(impl, name, results: dict):
    """
    Runs the benchmark on the given implementation across multiple image sizes.

    results is a dictionary mapping from widths to a dictionary containing the measurements

    Results from this benchmark run are added to the dictionary for the widths under 'name'
    """

    for width in get_range():
        eprint("{}x{}".format(width, width))
        r = results.setdefault(width, {})
        r[name] = run_single(impl, input_scene, width, width)


def run_shape_benchmark(impl, temp_dir, shape, name, results):
    """
    Runs the benchmark on the given implementation on a scene with the given shape with increasing number of shapes

    results is a dictionary mapping from number of shapes to a dictionary containing the measurements

    Results from this benchmark run are added to the dictionary for the number of shapes under 'name'
    """

    # Fixed size
    width = 1920
    height = 1080

    for num_shapes in get_range():
        eprint(f"{num_shapes} x {shape} @ {width}x{height}")
        fname = f"{impl}-{shape}-{num_shapes}-{width}x{height}-{name}.json"
        scene_path = temp_dir / fname

        if shape == "all":
            subprocess.run([str(all_shapes_scene_creator), str(
                num_shapes // 5), str(scene_path)], check=True)
        else:
            subprocess.run([str(shape_scene_creator), shape, str(
                num_shapes), str(scene_path)], check=True)

        r = results.setdefault(num_shapes, {})
        r[name] = run_single(impl, scene_path, width, height)

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
        raise subprocess.CalledProcessError(
            p.returncode, p.args, p.stdout, p.stderr)


def print_line(line, quotes=False):
    for l in line:
        if quotes:
            print('\t"{}"'.format(l), end='')
        else:
            print('\t{}'.format(l), end='')

    print('')


def print_results(results, x_name, flags):
    header = [x_name, "flops"]

    for flag in flags:
        header.extend([flag +
                       " - seconds", flag +
                       " - cycles", flag +
                       " - flops/cycle"])

    print_line(header, True)

    for x_val, measurements in results.items():
        num_flops = int(measurements["Flops"]["Flops"])

        line = [x_val, num_flops]

        for flag in flags:
            r = measurements[flag]
            usec = float(r["Microseconds"])
            cycles = float(r["Cycles"])
            line.extend([usec / 1e6, cycles, num_flops / cycles])

        print_line(line, False)


def main(temp_dir):

    flags = [
        "-O2",
        "-O3 -fno-tree-vectorize",
        "-O3",
        "-Ofast",
        "-Ofast -march=native"]

    if len(sys.argv) < 6:
        eprint(f"Usage: {sys.argv[0]} <impl> <type> <x-base> <x-end> <x-step>")
        eprint("\t<type> must be one of the following: 'size', 'all', 'box', 'sphere', 'cone', 'torus', 'octahedron'")
        sys.exit(1)

    impl = sys.argv[1]
    benchmark_type = sys.argv[2]
    global x_base
    global x_end
    global x_step
    x_base = int(sys.argv[3])
    x_end = int(sys.argv[4])
    x_step = int(sys.argv[5])

    if benchmark_type == "size":
        flops_runner = lambda results: run_size_benchmark(
            impl, "Flops", results)
        benchmark_runner = lambda flag, results: run_size_benchmark(
            impl, flag, results)
        x_name = "width/height"
    else:
        flops_runner = lambda results: run_shape_benchmark(
            impl, temp_dir, benchmark_type, "Flops", results)
        benchmark_runner = lambda flag, results: run_shape_benchmark(
            impl, temp_dir, benchmark_type, flag, results)
        x_name = f"number of shapes ({benchmark_type})"

    eprint("Flops Counter")
    # Count flops first
    build("", True)

    results = OrderedDict()
    flops_runner(results)

    for flag in flags:
        build(flag, False)
        benchmark_runner(flag, results)

    print_results(results, x_name, flags)

if __name__ == "__main__":
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            main(Path(temp_dir))
    except subprocess.CalledProcessError as e:
        print(e.stdout.decode())
        print(f"\033[31;1m{e.stderr.decode()}\033[0m")
        raise e

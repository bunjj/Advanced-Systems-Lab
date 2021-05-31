#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys
import tempfile
import json

script_dir = Path(__file__).parent
build_exec = script_dir / "build.sh"
run_exec = (script_dir / ".." / "build" / "main").resolve()
input_scene = script_dir / ".." / "scenes" / "scene0.json"
shape_scene_creator = script_dir / ".." / "benchmark" / "create_scene.py"
all_shapes_scene_creator = script_dir / ".." / \
    "benchmark" / "create_scene_all_shapes.py"

# Fields to collect from the program output
fields = [
    "Flops",
    "Cycles",
    "Microseconds",
    "ADD",
    "MUL",
    "FMA",
    "DIV",
    "SQRT",
    "ABS",
    "CMP",
    "MAX",
    "POW",
    "TAN",
    "SPHERE",
    "SPHERE_R",
    "PLANE",
    "BOX",
    "BOX_R",
    "TORUS",
    "TORUS_R",
    "CONE",
    "CONE_R",
    "OCTA",
    "OCTA_R",
    "SPHERE_N",
    "PLANE_N",
    "BOX_N",
    "TORUS_N",
    "CONE_N",
    "OCTA_N"]


datapoints = []

shape_scenes = {}

impls = [
    "ref",
    "opt0",
    "opt1",
    "opt3",
    "opt4",
    "opt5",
    "vec1",
    "vec2",
    "vec3",
    "vec4"]

shapes = ["all", "box", "sphere", "cone", "torus", "octahedron"]
bench_types = ["size"] + shapes
# flags = [
#     "-O2",
#     "-O3 -fno-tree-vectorize",
#     "-O3",
#     "-Ofast",
#     "-Ofast -march=native"]

flags = ["-Ofast"]


def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def get_range(bench_type):
    if bench_type == "size":
        x_base = 200
        x_end = 4200
        x_step = 200
    elif bench_type in ["all", "box", "sphere", "cone", "torus", "octahedron"]:
        x_base = 10
        x_end = 110
        x_step = 10
    else:
        raise ValueError(bench_type)

    return range(x_base, x_end, x_step)


def gen_shape_scenes(temp_dir):

    for shape in shapes:
        for num_shapes in get_range(shape):
            fname = f"{shape}-{num_shapes}.json"
            scene_path = temp_dir / fname

            if shape == "all":
                subprocess.run([str(all_shapes_scene_creator), str(
                    num_shapes // 5), str(scene_path)], check=True)
            else:
                subprocess.run([str(shape_scene_creator), shape, str(
                    num_shapes), str(scene_path)], check=True)

            shape_scenes.setdefault(shape, {})[num_shapes] = scene_path


def build(flag, instrument):
    cmd = [str(build_exec), "Release", flag]

    if instrument:
        cmd.append("-DINSTRUMENT=ON")
    else:
        cmd.append("-DINSTRUMENT=OFF")

    eprint(" ".join(cmd))

    p = subprocess.run(cmd, capture_output=True)

    if p.returncode != 0:
        eprint(p.stdout.decode())
        eprint(f"\033[31;1m{e.stderr.decode()}\033[0m")
        raise subprocess.CalledProcessError(
            p.returncode, p.args, p.stdout, p.stderr)


def run_single(impl, scene, width, height, datapoint_base: dict):
    p = subprocess.run([str(run_exec), impl, str(scene), str(
        width), str(height)], check=True, capture_output=True)

    datapoint = datapoint_base.copy()
    datapoint["scene"] = str(scene)
    datapoint["width"] = width
    datapoint["height"] = height

    for line in p.stderr.decode().split('\n'):
        if ': ' in line:
            key, value = line.split(': ')
            key = key.strip()
            value = value.strip()

            if key in fields:
                datapoint[key] = value

    datapoints.append(datapoint)


def run_size_benchmark(impl, datapoint_base: dict):
    """
    Runs the benchmark on the given implementation across multiple image sizes.
    """

    for width in get_range("size"):
        eprint("{}x{}".format(width, width))
        datapoint_base2 = datapoint_base.copy()
        run_single(impl, input_scene, width, width, datapoint_base2)


def run_shape_benchmark(impl, shape, datapoint_base: dict):
    """
    Runs the benchmark on the given implementation on a scene with the given shape with increasing number of shapes
    """

    # Fixed size
    width = 1920
    height = 1080

    for num_shapes in get_range(shape):
        eprint(f"{num_shapes} x {shape} @ {width}x{height}")
        datapoint_base2 = datapoint_base.copy()
        datapoint_base2["num_shapes"] = num_shapes

        scene_path = shape_scenes[shape][num_shapes]

        run_single(impl, scene_path, width, height, datapoint_base2)


def run_with_flags(flag, do_instrument, datapoint_base: dict):
    build(flag, do_instrument)
    for impl in impls:
        for bench_type in bench_types:
            print_header(
                f"Running benchmark type '{bench_type}' for implementation '{impl}' with flags '{flag}'")

            datapoint_base2 = datapoint_base.copy()
            datapoint_base2.update({
                "flags": flag,
                "impl": impl,
                "type": bench_type})

            if bench_type == "size":
                run_size_benchmark(impl, datapoint_base2)
            else:
                run_shape_benchmark(impl, bench_type, datapoint_base2)


def print_header(name):
    eprint(f"\033[32;1m{name}\033[0m")


def main(temp_dir):

    print_header("Generate Shape Scenes")
    gen_shape_scenes(temp_dir)

    print_header("Flops Counter")
    # Count flops first
    run_with_flags("", True, {"has_flops": True})

    for flag in flags:
        run_with_flags(flag, False, {"has_flops": False})

    print(json.dumps(datapoints, indent=4))


if __name__ == "__main__":
    try:
        with tempfile.TemporaryDirectory() as temp_dir:
            main(Path(temp_dir))
    except subprocess.CalledProcessError as e:
        eprint(e.stdout.decode())
        eprint(f"\033[31;1m{e.stderr.decode()}\033[0m")
        raise e

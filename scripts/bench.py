#!/usr/bin/env python3

from pathlib import Path
import subprocess
import sys
import tempfile
import json

class Config:
    def __init__(self, compilers, flags, impls, bench_types, shapes, range_func):
        self._compilers = compilers
        self._flags = flags
        self._impls = impls
        self._bench_types = bench_types
        self._shapes = shapes
        self._range_func = range_func

    @property
    def compilers(self):
        return self._compilers

    @property
    def flags(self):
        return self._flags

    @property
    def impls(self):
        return self._impls

    @property
    def bench_types(self):
        return self._bench_types

    @property
    def shapes(self):
        return self._shapes

    @property
    def range_func(self):
        return self._range_func
        

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

def eprint(*args, **kwargs):
    print(*args, file=sys.stderr, **kwargs)


def get_range(bench_type):
    if bench_type == "size":
        x_base = 200
        x_end = 4200
        x_step = 200
    elif bench_type in ["all", "box", "sphere", "cone", "torus", "octahedron"]:
        x_base = 10
        x_end = 371
        x_step = 40
    else:
        raise ValueError(bench_type)

    return range(x_base, x_end, x_step)


def gen_shape_scenes(conf : Config, temp_dir):

    for shape in conf.shapes:
        for num_shapes in conf.range_func(shape):
            fname = f"{shape}-{num_shapes}.json"
            scene_path = temp_dir / fname

            if scene_path.exists():
                eprint(f"\033[34mReusing scene file {str(scene_path)}\033[0m")
            elif shape == "all":
                subprocess.run([str(all_shapes_scene_creator), str(
                    num_shapes // 5), str(scene_path)], check=True)
            else:
                subprocess.run([str(shape_scene_creator), shape, str(
                    num_shapes), str(scene_path)], check=True)

            shape_scenes.setdefault(shape, {})[num_shapes] = scene_path


def build(compiler, flag, instrument):
    cmd = [str(build_exec), "Release", flag]

    if instrument:
        cmd.append("-DINSTRUMENT=ON")
    else:
        cmd.append("-DINSTRUMENT=OFF")

    cmd.append(f"-DCMAKE_CXX_COMPILER={compiler}");

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


def run_size_benchmark(conf : Config, impl, datapoint_base: dict):
    """
    Runs the benchmark on the given implementation across multiple image sizes.
    """

    for width in conf.range_func("size"):
        eprint("{}x{}".format(width, width))
        datapoint_base2 = datapoint_base.copy()
        run_single(impl, input_scene, width, width, datapoint_base2)


def run_shape_benchmark(conf : Config, impl, shape, datapoint_base: dict):
    """
    Runs the benchmark on the given implementation on a scene with the given shape with increasing number of shapes
    """

    # Fixed size
    width = 1920
    height = 1080

    for num_shapes in conf.range_func(shape):
        eprint(f"{num_shapes} x {shape} @ {width}x{height}")
        datapoint_base2 = datapoint_base.copy()
        datapoint_base2["num_shapes"] = num_shapes

        scene_path = shape_scenes[shape][num_shapes]

        run_single(impl, scene_path, width, height, datapoint_base2)


def run_with_flags(conf : Config, compiler, flag, do_instrument, datapoint_base: dict):
    build(compiler, flag, do_instrument)
    for impl in conf.impls:
        for bench_type in conf.bench_types:
            print_header(
                f"Running benchmark type '{bench_type}' for implementation '{impl}' with flags '{flag}'")

            datapoint_base2 = datapoint_base.copy()
            datapoint_base2.update({
                "compiler": compiler,
                "flags": flag,
                "impl": impl,
                "type": bench_type})

            if bench_type == "size":
                run_size_benchmark(conf, impl, datapoint_base2)
            else:
                run_shape_benchmark(conf, impl, bench_type, datapoint_base2)


def print_header(name):
    eprint(f"\033[32;1m{name}\033[0m")


def main(conf, temp_dir):

    print_header("Generate Shape Scenes")
    gen_shape_scenes(conf, temp_dir)

    print_header("Flops Counter")
    # Count flops first
    run_with_flags(conf, conf.compilers[0], "", True, {"has_flops": True})

    for compiler in conf.compilers:
        for flag in conf.flags:
            run_with_flags(conf, compiler, flag, False, {"has_flops": False})

    return datapoints


if __name__ == "__main__":

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

    compilers = ["g++", "clang++"]

    conf = Config(compilers, flags, impls, bench_types, shapes, get_range)

    try:
        if len(sys.argv) == 2:
            d = Path(sys.argv[1])
            if not d.exists():
                d.mkdir(parents = True)

            assert(d.is_dir())

            data = main(conf, d)
        else:
            with tempfile.TemporaryDirectory() as temp_dir:
                data = main(conf, Path(temp_dir))

        print(json.dumps(data, indent=4))
    except subprocess.CalledProcessError as e:
        eprint(e.stdout.decode())
        eprint(f"\033[31;1m{e.stderr.decode()}\033[0m")
        raise e

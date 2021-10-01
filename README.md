# Sphere Tracing - ASL Project

*By Felix Sarnthein, Irfan Bunjaku, Lukas Heimes, and Patrick Ziegler*

Project Work for MSc Advanced Systems Lab course offered at ETH Zurich (263-0007-00L). For details on implementation, algorithm, optimizations and evaluation check out **03_report.pdf**.

## Compiling

This project uses the standard cmake build process:

```bash
mkdir build
cd build
cmake ..
make -j
```

Several options can be set for the `cmake` command:

* `-DINSTRUMENT=ON`: Turns on instrumentation to count flops
* `-DSANITIZE=ON`: Turns on the address and undefined behavior sanitizer for debugging.
* `-DCMAKE_BUILD_TYPE=Release`/`-DCMAKE_BUILD_TYPE=Debug` Use release mode or debug mode.

On macOS:
* `-DCMAKE_C_COMPILER=/usr/local/bin/gcc-10` and `-DCMAKE_CXX_COMPILER=/usr/local/bin/g++-10`: Actually use gcc instead of Apple clang. `gcc-10` can be installed using Homebrew.

For benchmarking 

```bash
cmake -DCMAKE_BUILD_TYPE=Release -DINSTRUMENT=OFF -DSANITIZE=OFF
```

should be used and for regular development

```bash
cmake -DCMAKE_BUILD_TYPE=Debug -DINSTRUMENT=OFF -DSANITIZE=OFF
```

The sanitizers should only be turned on for finding bugs and as a sanity check
(to see if it can detect undefined behavior) on the code because it
dramatically increases compile and runtime.

## Running

Compiling will produce an executable called `main` in the `build` directory.
It takes four mandatory arguments and two optional ones:

```
./main <impl> <input> <width> <height> [<output>] [<reference>]
```

The first one is the name of the implementation to run. Each `impl_*` folder in
`src` is one such implementation, the name for it is everything after `impl_`.

The second parameter is the path to a JSON file containing a scene definition
followed by the desired width and height of the image in pixels.

The fifth argument is optional and denotes the path where the resulting ppm (or pfm)
image file should be written (any existing file will be overwritten). If the
argument is not specified, no image will be written.  The last optional argument
is the path to the reference ppm image (used for output validation).

## Create Benchmark scenes
Running the Python scripts will create benchmark scenes with random objects.
create_scene.py creates n scenes of type shape

```
python create_scene.py <shape> <n> <output>
```

create_scene_all_shapes.py creates n shapes of each shape type.

```
python create_scene_all_shapes.py <n> <output>
```

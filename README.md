# Sphere Tracing - ASL Project

*By Felix Sarnthein, Irfan Bunjaku, Lukas Heimes, and Patrick Ziegler*

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
It takes two mandatory arguments and one optional:

```
./main <input> <output> [<reference>]
```

The first one is the path to a JSON file containing a scene definition and the second one is the path where the resulting ppm image file should be written to (any existing file will be overwritten). The optional argument is the path to the reference ppm image (used for output validation).

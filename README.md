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
It takes two arguments: 

```
./main <input> <output>
```

The first one is the path to a JSON file containing a scene definition and the second one is the path where the resulting ppm image file should be written too (Any existing file will be overwritten).

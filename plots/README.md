# Plots for our sphere tracing optimizations

Use `./scripts/bench_runner.py` to run all benchmarks across all implementations. This
generates a `bench.json` and a `bench-compiler-flags.json` file in the `plots` directory.

You can then use the `./scripts/gen_plot_data.py` script to generate gnuplot data
files from these results:

```bash
./scripts/gen_plot_data.py plots/bench.json
./scripts/gen_plot_data.py plots/bench-compiler-flags.json flags
```

You can then generate plots using `gnuplot`:

```bash
gnuplot -c size-time.plt <bench-type> <y-index> <y-label> <output-basename>
```

This will generate a single plots with one line per optimization.
* `<bench-type>` selects which benchmark should be plotted (`size`, `all`, `box`, ...)
* `<y-index>` is the 1-based index of the column in the .dat files to use (`i` is the 0-based index of the compiler flags):
    * `2`: Selects Gflops
    * `3 * i + 3`: Seconds
    * `3 * i + 4`: Cycles
    * `3 * i + 5`: flops/cycle
* `<y-label>` sets the label for the y axis
* `<output-basename>` Name of the output file without extension

The `plots` directory also contains a host of other `.plt` gnuplot files that can be used to generate other plots.
The individual files have usage instructions at the top

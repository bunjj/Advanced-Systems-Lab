# Plots for our sphere tracing optimizations

Use `scripts/bench.py` to run all benchmarks across all implementations. This
generates a `bench.json` file in the `plots` directory.

You can then use the `scripts/gen_plot_data.py` script to generate gnuplot data
files from these results.

You can then generate plots using `gnuplot`:

```bash
gnuplot -c size-time.plt <bench-type> <y-index> <y-label>
```

This will generate a single plots with one line per optimization.
* `<bench-type>` selects which benchmark should be plotted (`size`, `all`, `box`, ...)
* `<y-index>` is the 1-based index of the column in the .dat files to use (`i` is the 0-based index of the compiler flags):
    * `2`: Selects Gflops
    * `3 * i + 3`: Seconds
    * `3 * i + 4`: Cycles
    * `3 * i + 5`: flops/cycle
* `<y-label>` sets the label for the y axis
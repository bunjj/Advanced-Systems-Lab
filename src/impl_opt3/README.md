# Opt3

## Changes (w.r.t. opt1)

Added short circuit termination in the distance functions for all shapes to avoid unneccessary sqrts.
Further replaced 4x4 matrices by 3x3 matrices.


## Timing

| Benchmark  | opt1 [s] | opt3 [s] | Improvement |
|------------|------:|------:|---------------:|
|`scene0`    |  3.64 |  3.39 | -0.25  (+6.8%) |
|`all`       | 59.04 | 52.18 | -6.86 (+11.6%) |
|`boxes`     | 36.26 | 32.88 | -3.38  (+9.3%) |
|`cones`     | 33.76 | 29.88 | -3.88 (+11.4%) |
|`octahedra` | 36.34 | 35.71 | -0.63  (+1.7%) |
|`spheres`   | 32.81 | 26.66 | -6.15 (+18.7%) |
|`tori`      | 40.22 | 31.76 | -8.46 (+21.0%) |

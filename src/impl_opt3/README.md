# Opt3

## Changes (w.r.t. opt1)

Added short circuit termination in the distance functions for all shapes to avoid unneccessary sqrts.
Further replaced 4x4 matrices by 3x3 matrices.

Short circuit termination extended. Precompute a radius r for each shape type and compute distance function similar to spheres to determine
if the object is too far away. Makes computation of exact distance functions in many cases obsolete. Bajor improvement to previous adjustments
for complex shapes

## Timing

| Benchmark  | opt1 [s] | opt3 [s] |     Improvement |
|------------|---------:|---------:|----------------:|
|`scene0`    |     3.67 |     2.47 |  -1.2s (+32.7%) |
|`all`       |    55.27 |    31.20 | -24.0s (+43.5%) |
|`boxes`     |    31.98 |    20.22 | -11.8s (+36.8%) |
|`cones`     |    31.92 |    12.90 | -19.0s (+59.6%) |
|`octahedra` |    36.50 |    18.64 | -17.9s (+48.9%) |
|`spheres`   |    27.37 |    26.42 |  -1.0s  (+3.5%) |
|`tori`      |    38.55 |    24.43 | -14.1s (+36.5%) |


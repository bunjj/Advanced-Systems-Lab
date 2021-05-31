# vec4

## Changes (w.r.t. vec3)

Add dummy objects to scene if the number of objects for a shape is not divisible by 8 and the remainder is > 2.

## Timing

| Benchmark  | vec3 [s] | vec4 [s] | Improvement         |
|------------|---------:|---------:|--------------------:|
|`scene0`    |     2.03 |     2.05 | -                   |
|`all`       |     7.88 |     5.40 | -2.5s (31.5%, 1.5x) |
|`boxes`     |     3.41 |     3.42 | -                   |
|`cones`     |     3.91 |     3.23 | -0.7s (17.4%, 1.2x) |
|`octahedra` |     3.21 |     3.29 | -                   |
|`spheres`   |     5.29 |     5.03 | -0.3s  (4.9%, 1.1x) |
|`tori`      |     4.03 |     3.58 | -0.5s (11.2%, 1.1x) |

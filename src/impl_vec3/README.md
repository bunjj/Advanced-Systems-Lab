# Vec3

## Changes (w.r.t. vec2)

* Adapted all changes in opt4 and opt5 into the vectorized version.
* Vectorized the finding of the minimum distance from 8 distances
* Used aligned loads and stores
* Vectorized ray precomputation

## Timing

| Benchmark  | vec2 [s] | vec3 [s] | Improvement         |
|------------|---------:|---------:|--------------------:|
|`scene0`    |     2.52 |     1.97 | -0.6s (21.8%, 1.3x) |
|`all`       |    11.65 |     7.69 | -4.0s (34.0%, 1.5x) |
|`boxes`     |     4.95 |     3.25 | -1.7s (34.3%, 1.5x) |
|`cones`     |     5.79 |     3.79 | -2.0s (34.5%, 1.5x) |
|`octahedra` |     4.86 |     3.12 | -1.7s (35.8%, 1.6x) |
|`spheres`   |     6.82 |     5.23 | -1.6s (23.3%, 1.3x) |
|`tori`      |     5.62 |     3.90 | -1.7s (30.6%, 1.4x) |

| Benchmark  | opt5 [s] | vec3 [s] | Improvement          |
|------------|---------:|---------:|---------------------:|
|`scene0`    |     1.58 |     1.97 |  +0.4s (-24.7% 0.8x) |
|`all`       |    16.29 |     7.69 |  -8.6s  (52.8% 2.1x) |
|`boxes`     |     9.87 |     3.25 |  -6.6s  (67.1% 3.0x) |
|`cones`     |     6.83 |     3.79 |  -3.0s  (44.5% 1.8x) |
|`octahedra` |     9.08 |     3.12 |  -6.0s  (65.6% 2.9x) |
|`spheres`   |    19.86 |     5.23 | -14.6s  (73.7% 3.8x) |
|`tori`      |    12.25 |     3.90 |  -8.4s  (68.2% 3.1x) |

| Benchmark  |  ref [s] | vec3 [s] | Improvement |
|------------|---------:|---------:|------------:|
|`scene0`    |     7.49 |     1.97 |        3.8x |
|`all`       |   113.75 |     7.69 |       14.8x |
|`boxes`     |    62.63 |     3.25 |       19.3x |
|`cones`     |    49.17 |     3.79 |       13.0x |
|`octahedra` |    68.45 |     3.12 |       21.9x |
|`spheres`   |   128.54 |     5.23 |       24.6x |
|`tori`      |    88.66 |     3.90 |       22.7x |

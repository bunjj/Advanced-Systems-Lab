# opt4

## Changes (w.r.t. opt3)
1. Move the computation of 'World->Object transform' from distance functions into rendering loop. Reuse this newly created p_obj for distance and normal computation.
2. Introduce rays and accociated methods in prepraration for the ray precomputations.
3. Precompute 'World->Object transform' for the input ray in render loop. There are now many representations of the same ray in all object coordinates. This computation was done for every distance query and is now only done once per object. 


## Timing

| Benchmark  | opt3 [s] | opt4 [s] |    Improvement |
|------------|---------:|---------:|---------------:|
|`scene0`    |     2.47 |     1.80 |  -0.7s (37.2%) |
|`all`       |    29.53 |    18.95 | -10.6s (55.8%) |
|`boxes`     |    19.96 |    12.39 |  -7.6s (61.1%) |
|`cones`     |    12.80 |     8.06 |  -4.7s (58.8%) |
|`octahedra` |    18.34 |    10.97 |  -7.4s (67.1%) |
|`spheres`   |    23.98 |    24.24 |  +0.2s (-1.1%) |
|`tori`      |    24.25 |    14.57 |  -9.7s (66.4%) |


# opt4

## Changes (w.r.t. opt3)
1. Move the computation of 'World->Object transform' from distance functions into rendering loop. Reuse this newly created p_obj for distance and normal computation.
2. Introduce rays and accociated methods in prepraration for the ray precomputations.
3. Precompute 'World->Object transform' for the input ray in render loop. There are now many representations of the same ray in all object coordinates. This computation was done for every distance query and is now only done once per object. 


## Timing

| Benchmark  | opt3 [s] | opt4 [s] |     Improvement |
|------------|---------:|---------:|----------------:|
|`scene0`    |     2.63 |    2.19  |   -0.4s (16.7%) |
|`all`       |    34.26 |    24.18 |  -10.1s (29.4%) |
|`boxes`     |    23.27 |    15.66 |   -7.6s (32.7%) |
|`cones`     |    14.39 |     9.59 |   -4.8s (33.4%) |
|`octahedra` |    20.70 |    13.94 |   -6.8s (32.7%) |
|`spheres`   |    27.66 |    28.54 |   +0.9s (-3.2%) |
|`tori`      |    27.65 |    19.13 |   -8.5s (30.8%) |


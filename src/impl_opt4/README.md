# opt4

## Changes (w.r.t. opt3)
1. Move the computation of 'World->Object transform' from distance functions into rendering loop. Reuse this newly created p_obj for distance and normal computation.
2. Introduce rays and accociated methods in prepraration for the ray precomputations.
3. Precompute 'World->Object transform' for the input ray in render loop. There are now many representations of the same ray in all object coordinates. This computation was done for every distance query and is now only done once per object. 

I don't understand why there is no speed-up. When I tested it before the introduction of early termination, I recall an improvement of 20% and early termination shouldn't be save those matrix multiplications.

## Timing

| Benchmark  | opt3 [s] | opt4 [s] |     Improvement |
|------------|---------:|---------:|----------------:|
|`scene0`    |     2.71 |     2.88 |  +0.2s  (-6.3%) |
|`all`       |    34.36 |    35.59 |  +1.2s  (-3.6%) |
|`boxes`     |    24.93 |    22.83 |  -2.1s (+08.4%) |
|`cones`     |    14.65 |    15.94 |  +1.3s  (-8.8%) |
|`octahedra` |    21.44 |    27.07 |  +5.6s (-26.3%) |
|`spheres`   |    28.16 |    32.66 |  +5.5s (-16.0%) |
|`tori`      |    28.27 |    28.27 |  -0.0s   (0.0%) |


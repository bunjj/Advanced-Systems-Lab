# Vec2

## Changes (w.r.t. vec1)

Adds early termination by first computing the distance to an enclosing sphere. This is the vectorized version of (full) `opt3`.

Furthermore, the early termination has been simplified: The distance functions are no longer split into two functions. The single vectorized distance function indicates to the caller whether early termination was applied. While the caller still has to decide whether or not to process the results, the caller no longer needs to call the second function.

## Timing

(measured by Lukas)

| Benchmark  | vec1 [s] | vec2 [s] | Improvement |
|------------|------:|------:|---------------:|
|`scene0`    |  2.81 |  2.01 | 1.40x |
|`all`       | 15.26 | 10.76 | 1.42x |
|`boxes`     |  5.37 |  4.29 | 1.25x |
|`cones`     |  9.31 |  5.04 | 1.85x |
|`octahedra` |  7.15 |  4.28 | 1.67x |
|`spheres`   |  5.10 |  5.10 | - |
|`tori`      |  5.97 |  5.64 | 1.06x |

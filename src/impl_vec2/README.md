# Vec2

## Changes (w.r.t. vec1)

Adds early termination by first computing the distance to an enclosing sphere. This is the vectorized version of (full) `opt3`.

## Timing

(measured by Lukas)

| Benchmark  | vec1 [s] | vec2 [s] | Improvement |
|------------|------:|------:|---------------:|
|`scene0`    |  2.85 |  2.01 | 1.42x |
|`all`       | 15.84 |  9.85 | 1.61x |
|`boxes`     |  5.36 |  4.26 | 1.26x |
|`cones`     |  9.30 |  5.25 | 1.77x |
|`octahedra` |  9.12 |  4.15 | 2.20x |
|`spheres`   |  5.21 |  5.84 | 0.89x (-) |
|`tori`      |  6.28 |  5.34 | 1.18x |

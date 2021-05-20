# Vec1

## Changes (w.r.t. opt3)

Vectorized distance functions, i.e., operate on 8 shapes of the same type in parallel.
Note that this includes the early termination _to avoid sqrt computations_, but does not include trick with the surrounding sphere that was introduced later in `opt3`.


## Timing

(measured by Lukas)

| Benchmark  | opt3 [s] | vec1 [s] | Improvement |
|------------|------:|------:|---------------:|
|`scene0`    |  2.59 |  2.81 | 0.92x (-)
|`all`       | 38.77 | 17.30 | 2.24x
|`boxes`     | 26.06 |  5.26 | 4.95x
|`cones`     | 22.35 |  9.37 | 2.39x
|`octahedra` | 29.20 |  7.54 | 3.88x
|`spheres`   | 17.91 |  5.45 | 3.29x
|`tori`      | 23.57 |  7.08 | 3.33x

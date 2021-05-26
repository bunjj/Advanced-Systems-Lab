# opt5

## Changes (w.r.t. opt4)

* Optimizes distance and normal functions. It is not always worth it to
  precompute certain properties and store them in the shape. Likely because the
  compiler can better optimize when the computation is in the function.
* Remove all m44 from shape structs. Fewer data that needs to be moved around
* Smaller optimization in the sphere tracing: Precompute, reduce local variables, strength reduction.

---

Try to add accumulators for the different shapes
* Sphere: Does reduce pure-sphere scene execution time by 1 second, but the overhead makes all other scenes slower.
* Box: Only accumulators in sphere_trace are worth it, not in the shadow tracing
* Torus, Cone, Octahedron: Accumulators have no benefit.

In the end the accumulators decreased the overall performance except for the
pure-sphere scene, so I decided to not include them

## Timing

| Benchmark  | opt4 [s] | opt5 [s] |    Improvement |
|------------|---------:|---------:|---------------:|
|`scene0`    |     2.32 |     2.12 |  -0.2s  (9.4%) |
|`all`       |    23.09 |    22.39 |  -0.7s  (3.1%) |
|`boxes`     |    15.21 |    13.88 |  -1.3s  (9.6%) |
|`cones`     |    11.06 |     9.52 |  -1.5s (16.2%) |
|`octahedra` |    12.53 |    12.42 |  -0.1s  (0.9%) |
|`spheres`   |    27.52 |    27.50 |  -0.0s  (0.0%) |
|`tori`      |    17.48 |    17.48 |  -0.0s  (0.0%) |


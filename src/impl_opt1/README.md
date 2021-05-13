# Opt1

## Changes (w.r.t. reference)
- Data representation:
    - Get rid of generic shape struct, different shapes store only what they actually need.
    - Scene now has arrays for each shape type, i.e., a scene is a struct containing pointers to arrays of structs.
- Algorithm:
    - Restructure s.t. the type of shape is known statically, i.e., first process all spheres, then all boxes, etc.
- Other:
    - Move distance and normal functions into header file s.t. the compiler can inline them.

## Timings

TODO

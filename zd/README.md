# Pipeline

After compiling all 3 programs, to calculate a particular *even* d, there are 3 steps:

1. `./canonCount d - 1 > c{$d-1}.txt`.
2. `./SATCount d - 2 > s{$d-2}.txt`.
3. `./finalCount d s{$d-2}.txt c{$d-1}.txt`.

For example: with `d = 4`, the result is 2970.

## CanonCount

Given a dimension `d`, the program `canonCount` outputs all canonical choices of `2^(d - 1)` vertices from `2^d` as a `2^d` bit integer followed by the number of distinct automorphisms of that vertex choice.

## SATCount

Given a dimension `d`, the program `SATCount` outputs a sorted list of triples of the form `(v, e, count)` where:

* `v` is a `2^d` bit integer, each bit corresponding to the parity of the difference in in/out degree of the vertex
* `e` is a `2^d` bit integer, each bit corresponding to the direction of an edge extending from a particular vertex into space.
* `count` is the number of valid ways to pick edge directions in the d-box given the constraints on `v` and `e`.

## FinalCount

Using the results from `./canonCount d-1` and `./SATCount d-2`, `FinalCount` calculates the number of Eulerian Orientations for a d-box.

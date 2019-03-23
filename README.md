# GraphModularDecomposition.jl

This package implements [modular decomposition](https://en.wikipedia.org/wiki/Modular_decomposition)
of directed graphs using algorithms from the following three papers:

* McConnell & Montgolfier 2004: _"Linear-time modular decomposition of directed graphs"_
* Capelle, Habib & de Montgolfier 2002: _"Graph decompositions and factorizing permutations"_
* Habib, Paul & Viennot 1999: _"Partition refinement techniques: an interesting algorithmic tool kit"_

The implementation as it currently stands is not linear time for fairly silly reasons. Specifically,
I used a brute force algorithm for computing the overlap components of two strong module trees; a
linear time algorithm exists, as described in Dalhaus 1998, "Parallel algorithms for hierarchical
clustering and applications to split decomposition and parity graph recognition", but I didn't
bother to implement it. There are also parts of strong module tree construction from a factorizing
permutation that may be super-linear unnecessarily.

The strong module tree construction code uses an efficient and novel (to me at least) approach to
manipulating trees: instead of building tree data structures explicitly, since the tree operations
always keep the leaves in the same order, it works with counts of how many open and close
"parentheses" occur before and after each node, thereby implicitly maintaining a tree structure by
changing how nodes are parenthesized. This technique allows allocating linear storage once at the
start of the algorithm, rather than needing to continually allocate and free tree data structures.
The entire computation is done by incrementing and decrementing a few vectors of per-node "parens"
counters. Only at the end when the final tree structure is determined is an actual tree created.

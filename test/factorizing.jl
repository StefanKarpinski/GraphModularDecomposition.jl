# Regression tests for #5: symgraph_factorizing_permutation used to return
# permutations in which a strong module was not contiguous, which is the one
# thing a factorizing permutation is supposed to guarantee. The cause was
# add_pivot dropping half of a split part instead of carrying both halves
# forward, so some vertices were never used as pivots.

# graphs with a twin pair, and starting orders that used to tear the pair apart
const TORN_APART = [
    ([0 0 1 0 0 1 1 0 0 0
      0 0 1 0 0 1 1 0 0 0
      1 1 0 0 1 0 0 1 1 0
      0 0 0 0 0 1 1 1 1 0
      0 0 1 0 0 0 1 0 1 1
      1 1 0 1 0 0 1 1 0 0
      1 1 0 1 1 1 0 1 0 0
      0 0 1 1 0 1 1 0 0 0
      0 0 1 1 1 0 0 0 0 0
      0 0 0 0 1 0 0 0 0 0], [2, 4, 9, 6, 7, 5, 10, 3, 8, 1], [1, 2]),
    ([0 1 1 1 0 0 0 0 0 0
      1 0 1 1 1 0 0 1 1 0
      1 1 0 0 1 1 0 0 1 0
      1 1 0 0 1 0 0 0 1 0
      0 1 1 1 0 1 0 0 0 0
      0 0 1 0 1 0 1 1 0 1
      0 0 0 0 0 1 0 1 0 1
      0 1 0 0 0 1 1 0 0 0
      0 1 1 1 0 0 0 0 0 0
      0 0 0 0 0 1 1 0 0 0], [9, 8, 3, 5, 7, 10, 1, 4, 6, 2], [1, 9]),
]

@testset "factorizing permutations: regressions" begin
    for (G, order, M) in TORN_APART
        n = size(G, 1)
        @test G == G'
        @test is_module(G, M)
        p = symgraph_factorizing_permutation(G, copy(order))
        @test sort(p) == collect(1:n)
        @test diff(sort([findfirst(isequal(x), p) for x in M])) == [1]
        @test is_modular_permutation(G, p)
    end
end

# The bug only showed up on graphs carrying a module that the rest of the graph
# does not decompose neatly around: uniformly random graphs rarely have one at
# all, and recursively substituted graphs decompose too cleanly to trip it.
# Planting a twin pair in an otherwise random graph does reproduce it, but only
# for about one starting order in three thousand, so this checks a lot of them
# at once and asserts once. Against the code before #5 was fixed this seed
# produces 33 failures; a handful of orders would have found nothing.
@testset "factorizing permutations: planted twins" begin
    rng = MersenneTwister(0xf1e1d) # a private stream: see test/overlap.jl
    failures, checked, first_failure = 0, 0, nothing
    for _ = 1:400
        n = rand(rng, 9:11)
        G = zeros(Int, n, n)
        for i = 1:n-1, j = i+1:n
            G[i,j] = G[j,i] = rand(rng, 0:1)
        end
        for k = 3:n # vertices 1 and 2 are twins, so {1,2} is a module
            G[2,k] = G[k,2] = G[1,k]
        end
        modules = all_strong_modules(G)
        for _ = 1:250
            p = symgraph_factorizing_permutation(G, shuffle(rng, 1:n))
            checked += 1
            is_modular_permutation(G, p, modules=modules) && continue
            failures += 1
            first_failure === nothing && (first_failure = (copy(G), copy(p)))
        end
    end
    @test checked == 100_000
    failures == 0 ||
        @warn "non-modular permutation" G=first_failure[1] p=first_failure[2]
    @test failures == 0
end

# Regression for the digraph path. intersect_permutation builds a laminar family
# of sets and walks it to lay the vertices out. A strong module can fail to be
# one of those sets and still have to come out contiguous: below, {1,4} is a
# module, is not one of the sets, and stays together only because the leaves
# belonging to a node directly are not separated by that node's children.
@testset "factorizing permutations: digraph regressions" begin
    G = [0 0 1 0 0
         0 0 0 0 0
         0 1 0 0 0
         0 0 1 0 0
         0 1 1 0 0]
    @test is_module(G, [1, 4])
    p = digraph_factorizing_permutation(G)
    @test sort(p) == collect(1:5)
    @test diff(sort([findfirst(isequal(x), p) for x in (1, 4)])) == [1]
    @test is_modular_permutation(G, p)

    # That graph turned up about once in 50,000 random digraphs, so a sweep of a
    # few thousand would have proved nothing -- and did not, when it was tried.
    # This one reports two failures against the code the regression is for.
    rng = MersenneTwister(3) # a private stream: see test/overlap.jl
    failures, first_failure = 0, nothing
    for _ = 1:25_000
        n = rand(rng, 4:7)
        H = Int[i != j && rand(rng, Bool) for i = 1:n, j = 1:n]
        q = digraph_factorizing_permutation(H)
        sort(q) == collect(1:n) && is_modular_permutation(H, q) && continue
        failures += 1
        first_failure === nothing && (first_failure = (copy(H), copy(q)))
    end
    failures == 0 ||
        @warn "non-modular permutation" G=first_failure[1] p=first_failure[2]
    @test failures == 0
end

# LinearAlgebra.issymmetric can report a sparse matrix symmetric when it is not
# (JuliaSparse/SparseArrays.jl#748), and choosing the undirected algorithm for a
# digraph on that basis returns a permutation that is not factorizing at all.
# This is one of the matrices it gets wrong, and the exhaustive check below is
# what turned it up: it is 2 of the 4096 digraphs on four vertices, so nothing
# smaller than all of them was going to find it.
@testset "symmetry is tested without trusting issymmetric" begin
    G = [0 0 1 0
         0 0 0 1
         1 1 0 0
         0 1 0 0]
    @test G != transpose(G)
    @test !GraphModularDecomposition.is_symmetric(sparse(G))
    @test is_modular_permutation(G, graph_factorizing_permutation(sparse(G)))

    disagreements = 0
    pairs = [(i, j) for i = 1:4 for j = 1:4 if i != j]
    for bits = 0:(2^length(pairs) - 1)
        H = zeros(Int, 4, 4)
        for (b, (i, j)) in enumerate(pairs)
            (bits >> (b-1)) & 1 == 1 && (H[i,j] = 1)
        end
        repr(StrongModuleTree(H)) == repr(StrongModuleTree(sparse(H))) ||
            (disagreements += 1)
    end
    @test disagreements == 0
end

# Cutters are found two different ways, by scanning a dense matrix and by
# merging the relation lists of a sparse one, so the representation is now
# something the answer could depend on and must not.
@testset "cutters agree across representations" begin
    rng = MersenneTwister(0x5ca7) # a private stream: see test/overlap.jl
    for _ = 1:2000
        n = rand(rng, 2:14)
        G = zeros(Int, n, n)
        kind = rand(rng, 1:3)
        if kind == 1                                        # symmetric graph
            for i = 1:n-1, j = i+1:n
                G[i,j] = G[j,i] = rand(rng, 0:1)
            end
        elseif kind == 2                                    # digraph
            for i = 1:n, j = 1:n
                i != j && (G[i,j] = rand(rng, 0:1))
            end
        else                                                # 0/1/2 two-structure
            for i = 1:n-1, j = i+1:n
                G[i,j] = G[j,i] = rand(rng, 0:2)
            end
        end
        p = shuffle(rng, 1:n)
        @test repr(StrongModuleTree(G, copy(p))) ==
              repr(StrongModuleTree(sparse(G), copy(p)))
    end
end

# Regression for #13. sort_leaves! orders the children of a degenerate node of
# the H tree, and that order is what decides whether a module of G that is only
# a *union of siblings* comes out contiguous. Here {4,5,7} is such a module: it
# is a module of G, a weak module of the symmetric part, and a module of no tree
# the algorithm builds, so nothing but the child order keeps it together.
@testset "factorizing permutations: sibling unions (#13)" begin
    G = [0 0 0 0 0 1 0
         1 0 0 0 0 1 0
         0 0 0 0 0 0 0
         1 0 0 0 0 0 0
         1 0 0 0 0 0 0
         0 0 0 0 0 0 0
         1 0 0 0 1 0 0]
    @test is_module(G, [4, 5, 7])
    @test is_module(G, [5, 7])
    for rep in (G, sparse(G))
        p = graph_factorizing_permutation(rep)
        @test sort(p) == collect(1:7)
        @test diff(sort([findfirst(isequal(x), p) for x in (4,5,7)])) == [1, 1]
        @test is_modular_permutation(G, p)
    end

    # Breadth, but read the margin before trusting it: the defect turns up in
    # about one sparse digraph in forty thousand, and this seed finds exactly
    # one in twenty-five thousand. It is the graph above that is the real guard;
    # this is here to catch anything of the same shape that is not that graph.
    rng = MersenneTwister(7) # a private stream: see test/overlap.jl
    failures = 0
    for _ = 1:6_000
        n = rand(rng, 5:9)
        H = Int[i != j && rand(rng, Bool) for i = 1:n, j = 1:n] .* rand(rng, 0:1, n, n)
        q = graph_factorizing_permutation(H)
        sort(q) == collect(1:n) && is_modular_permutation(H, q) || (failures += 1)
    end
    @test failures == 0
end

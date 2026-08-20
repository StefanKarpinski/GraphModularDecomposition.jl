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

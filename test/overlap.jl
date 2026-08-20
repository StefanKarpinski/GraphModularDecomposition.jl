using GraphModularDecomposition.OverlapComponents
using GraphModularDecomposition.OverlapComponents: dahlhaus_max, lex_order, sortperm_desc

## brute force reference: build the overlap graph and flood fill it ##

function brute_overlap_classes(F::Vector{Vector{Int}})
    m = length(F)
    S = map(sort, F)
    adjacent = [[j for j = 1:m if j != i && overlap(S[i], S[j])] for i = 1:m]
    class, nclasses = zeros(Int, m), 0
    for i = 1:m
        class[i] == 0 || continue
        nclasses += 1
        stack = [i]
        while !isempty(stack)
            x = pop!(stack)
            class[x] == 0 || continue
            class[x] = nclasses
            append!(stack, adjacent[x])
        end
    end
    return class, nclasses
end

function brute_overlap_components(F::Vector{Vector{Int}})
    class, nclasses = brute_overlap_classes(F)
    U = [Int[] for _ = 1:nclasses]
    for (i, X) in enumerate(F)
        append!(U[class[i]], X)
    end
    return sort!(map(sort!∘unique!, U))
end

# Max(X) is the first set in LF order that overlaps X and is at least as big
function brute_max(F::Vector{Vector{Int}}, n::Int)
    sizes, S = map(length, F), map(sort, F)
    order = sortperm_desc(1:length(F), sizes, n)
    map(1:length(F)) do i
        j = findfirst(k -> sizes[k] ≥ sizes[i] && overlap(S[i], S[k]), order)
        j === nothing ? 0 : order[j]
    end
end

fast_max(F::Vector{Vector{Int}}, n::Int) =
    let sizes = map(length, F), order = sortperm_desc(1:length(F), sizes, n)
        P = lex_order(F, order, n)
        first(dahlhaus_max(F, order, sizes, P.elems, P.pos, n))
    end

# two set families induce the same partition of their common index set
same_classes(a::Vector{Int}, b::Vector{Int}) =
    Set(Set(findall(isequal(c), a)) for c in unique(a)) ==
    Set(Set(findall(isequal(c), b)) for c in unique(b))

## the worked example from Charbit et al. 2007, figures 1 and 2 ##

@testset "overlap classes: paper example" begin
    a, b, c, d, e, f, g, h, i, j, k, l = 1:12
    F = [[a,b], [b,c,d], [c,d,e,f,g,h], [d,e], [e,f,g,h], [f,g],
         [b,h], [j,k], [b,i,j,k,l], [k,l], [a,i]]
    X1, X2, X3, X4, X5, X6, X7, X8, X9, X10, X11 = 1:11

    order = sortperm_desc(1:11, map(length, F), 12)
    @test order == [X3, X9, X5, X2, X1, X4, X6, X7, X8, X10, X11] # figure 2 rows

    # figure 2 columns: a i l j k b c d h f g e -- f and g are indistinguishable
    P = lex_order(F, order, 12)
    @test P.elems[1:9] == [a, i, l, j, k, b, c, d, h]
    @test sort(P.elems[10:11]) == [f, g]
    @test P.elems[12] == e

    @test fast_max(F, 12) == brute_max(F, 12)
    @test fast_max(F, 12)[[X1, X8, X10, X11]] == [X9, X10, X8, X9] # stated in §3

    # X6 = {f,g} is nested inside X5 and X3 and meets nothing else; X8 and X10
    # overlap each other but are nested inside X9; everything else is one class
    class, nclasses = overlap_classes(F, 12)
    @test nclasses == 3
    @test same_classes(class, [1,1,1,1,1,2,1,3,1,3,1])
    @test sort(overlap_components(F, 12)) == [collect(1:12), [f,g], [j,k,l]]
    @test sort(overlap_components(F, 12)) == brute_overlap_components(F)
end

## randomized differential testing against the brute force reference ##

# these tests draw from a private stream rather than the global one, so that
# including this file does not shift what the randomized tests elsewhere see
const RNG = MersenneTwister(0xf1e1d)

rand_sets(n, m) = [sort!(shuffle(RNG, 1:n)[1:rand(RNG, 0:n)]) for _ = 1:m]
rand_small_sets(n, m) = [sort!(shuffle(RNG, 1:n)[1:rand(RNG, 1:min(3, n))]) for _ = 1:m]
rand_intervals(n, m) =
    [collect(UnitRange(minmax(rand(RNG, 1:n), rand(RNG, 1:n))...)) for _ = 1:m]

# a random hierarchy over 1:n: laminar, so no two of its sets ever overlap
function rand_laminar(n::Int)
    F = Vector{Int}[]
    function split(V)
        push!(F, sort(V))
        length(V) > 1 || return
        k = rand(RNG, 1:length(V)-1)
        split(V[1:k]); split(V[k+1:end])
    end
    split(shuffle(RNG, 1:n))
    return F
end

@testset "overlap classes: random families" begin
    families = Dict(
        "subsets"    => n -> rand_sets(n, rand(RNG, 1:12)),
        "small sets" => n -> rand_small_sets(n, rand(RNG, 1:20)),
        "intervals"  => n -> rand_intervals(n, rand(RNG, 1:15)),
        "laminar"    => n -> rand_laminar(n),
        # two hierarchies over one ground set: what modular decomposition feeds
        # to overlap_components, and the only case that reaches it in this package
        "two trees"  => n -> [rand_laminar(n); rand_laminar(n)],
    )
    for (name, generate) in families
        @testset "$name" begin
            for _ = 1:200
                n = rand(RNG, 1:14)
                F = generate(n)
                @test fast_max(F, n) == brute_max(F, n)
                @test same_classes(first(overlap_classes(F, n)),
                                   first(brute_overlap_classes(F)))
                @test sort(overlap_components(F, n)) == brute_overlap_components(F)
            end
        end
    end
end

@testset "overlap classes: degenerate families" begin
    @test overlap_components(Vector{Int}[], 0) == []
    @test overlap_components([Int[]], 0) == [[]]
    @test overlap_components([[1], [1], [1]], 1) == [[1], [1], [1]]
    @test overlap_components([[1,2,3], [1,2,3]], 3) == [[1,2,3], [1,2,3]]
    # {1,2} and {2,3} overlap; {1,2,3} overlaps neither, so it is a class of
    # its own -- two classes whose unions happen to coincide
    @test overlap_components([[1,2], [2,3], [1,2,3]], 3) == [[1,2,3], [1,2,3]]
    @test overlap_components([[1,2], [3,4]], 4) == [[1,2], [3,4]]
    # elements need not be dense in 1:n, nor start at 1
    @test overlap_components([[2,4], [4,6]]) == [[2,4,6]]
    @test overlap_components([[-3,0], [0,5]]) == [[-3,0,5]]
    @test overlap_components([["a","b"], ["b","c"], ["d"]]) == [["a","b","c"], ["d"]]
end

## the strong module tree entry point ##

@testset "overlap components of strong module trees" begin
    for _ = 1:200
        n = rand(RNG, 3:9)
        G = Int[i != j && rand(RNG, Bool) for i = 1:n, j = 1:n]
        Gs, Gd = G .| G', G .& G'
        s = StrongModuleTree(Gs, symgraph_factorizing_permutation(Gs))
        t = StrongModuleTree(Gd, symgraph_factorizing_permutation(Gd))
        M, N = strong_modules(s), strong_modules(t)
        got = sort!(GraphModularDecomposition.overlap_components(s, t, M, N))
        @test got == brute_overlap_components(M ∪ N)
    end
end

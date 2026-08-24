module GraphModularDecomposition

export StrongModuleTree, strong_modules,
    graph_factorizing_permutation,
    digraph_factorizing_permutation,
    symgraph_factorizing_permutation,
    tournament_factorizing_permutation

using LinearAlgebra
using LinearAlgebra: checksquare
using SparseArrays
using SparseArrays: nonzeros, nzrange, rowvals

include("OverlapComponents.jl")
using .OverlapComponents

include("StrongModuleTrees.jl")
using .StrongModuleTrees

const \ = setdiff

## iterating a vertex's neighbours ##

"""
    neighbors(G, x)

Iterate the vertices adjacent to `x`, reading column `x` of `G`.

This is the one graph operation the refinement below needs, and each matrix
representation can serve it in time proportional to the part of itself that it
has to touch: scanning a dense column is `O(n)`, which is that column's whole
size, while walking the stored entries of a sparse column is `O(deg x)`.
Supporting a new representation efficiently is one more method here.

The callers below are for undirected graphs, where column `x` and row `x` agree.
Column order is deliberate: it is what a `SparseMatrixCSC` stores contiguously,
and what a dense column-major `Matrix` reads without striding.
"""
neighbors(G::AbstractMatrix, x::Integer) =
    (y for y in axes(G, 1) if y != x && !iszero(G[y, x]))

function neighbors(G::SparseMatrixCSC, x::Integer)
    rows, vals = rowvals(G), nonzeros(G)
    (rows[k] for k in nzrange(G, x) if rows[k] != x && !iszero(vals[k]))
end

function symgraph_factorizing_permutation(
    G :: AbstractMatrix,
    V :: Vector{Int} = collect(1:checksquare(G)),
)
    P = [V]
    center::Int = 0
    pivots::Vector{Vector{Int}} = []
    modules::Vector{Vector{Int}} = []
    first_pivot = Dict{Vector{Int},Int}()

    # membership flags, so that a part can be split in one pass over it rather
    # than by intersecting it with the splitting set, and so that the pivot part
    # can be subtracted from a neighbourhood as the neighbourhood is built
    n = checksquare(G)
    in_splitter = falses(n)   # the set currently being split by
    in_pivot_part = falses(n) # the part the current pivot vertices came from
    in_V = falses(n)          # V may be a subset of the vertices of G
    for y in V
        in_V[y] = true
    end

    smaller_larger(A, B) = length(A) <= length(B) ? (A, B) : (B, A)

    function refine!(P, S, x)
        for y in S
            in_splitter[y] = true
        end
        i, between = 0, false
        while (i += 1) <= length(P)
            X = P[i]
            if center in X || x in X
                between = !between
                continue
            end
            # X ∩ S followed by X \ (X ∩ S) walks S once per part, so refining
            # the whole partition costs O(|P|·|S|); splitting X against the
            # membership flags instead costs O(|X|), and O(Σ|X|) over the
            # partition. element order within each half is the order within X,
            # exactly as intersect and setdiff give it
            k = count(y -> in_splitter[y], X)
            (k == 0 || k == length(X)) && continue
            Xₐ, Xₙ = Vector{Int}(undef, k), Vector{Int}(undef, length(X) - k)
            a = b = 0
            for y in X
                in_splitter[y] ? (Xₐ[a += 1] = y) : (Xₙ[b += 1] = y)
            end
            P[i] = Xₙ
            insert!(P, i + between, Xₐ)
            add_pivot(X, Xₐ, Xₙ)
            i += 1
        end
        for y in S
            in_splitter[y] = false
        end
    end

    function add_pivot(X, Xₐ, Xₙ)
        i = findfirst(isequal(X), pivots)
        if i !== nothing
            # every vertex of X was going to be used as a pivot, so every
            # vertex of both parts it just split into still has to be, and X
            # itself is no longer a part: L ← L ∪ {Xₐ, Xₙ} \ {X} in the
            # notation of Habib & Paul 2009, algorithm 2
            deleteat!(pivots, i)
            push!(pivots, Xₐ, Xₙ)
        else
            S, L = smaller_larger(Xₐ, Xₙ)
            push!(pivots, S)
            j = findfirst(isequal(X), modules)
            if j !== nothing
                modules[j] = L
            else
                push!(modules, L)
            end
        end
    end

    function partition_refinement!(P)
        while init_partition!(P)
            while !isempty(pivots)
                E = pop!(pivots)
                for y in E
                    in_pivot_part[y] = true
                end
                for x in E
                    # the neighbourhood of x minus the part it came from, built
                    # in one pass: O(deg x) given a representation that can list
                    # a vertex's neighbours, O(n) given one that cannot
                    S = [y for y in neighbors(G, x)
                         if in_V[y] && !in_pivot_part[y]]
                    refine!(P, S, x)
                end
                for y in E
                    in_pivot_part[y] = false
                end
            end
        end
    end

    function init_partition!(P)
        maximum(length, P) <= 1 && return false
        if isempty(modules)
            for (i, X) in enumerate(P)
                length(X) > 1 || continue
                x = get(first_pivot, X, first(X))
                # split X into the neighbours and non-neighbours of x, keeping
                # the order X had, in O(|X| + deg x) rather than by testing
                # every element of X against the graph
                for y in neighbors(G, x)
                    in_splitter[y] = true
                end
                A = [y for y in X if y != x && in_splitter[y]]
                N = [y for y in X if y != x && !in_splitter[y]]
                for y in neighbors(G, x)
                    in_splitter[y] = false
                end
                splice!(P, i, filter!(!isempty, [A, [x], N]))
                S, L = smaller_larger(A, N)
                center = x
                push!(pivots, S)
                push!(modules, L)
                break
            end
        else
            X = popfirst!(modules)
            x = first(X)
            push!(pivots, [x])
            first_pivot[X] = x
        end
        return true
    end

    partition_refinement!(P)
    return map(first, P)
end

function tournament_factorizing_permutation(
    G :: AbstractMatrix, # matrix representation of a tournament digraph
    V :: Vector{Int} = collect(1:checksquare(G)),
)
    n, P = length(V), [V]
    for x = V
        i = findfirst(C->x in C, P)
        C = P[i]
        B = filter(y->x != y && G[x,y] < G[y,x], C)
        A = filter(y->x != y && G[y,x] < G[x,y], C)
        splice!(P, i, filter!(!isempty, [B, [x], A]))
    end
    return map(first, P)
end

function intersect_permutation(
    V :: AbstractVector{E},   # vertices
    s :: StrongModuleTree{E}, # 1st strong module tree
    t :: StrongModuleTree{E}, # 2nd strong module tree
) where E
    Ms = strong_modules(s)
    Mt = strong_modules(t)
    # Ms ∪ Mt drops the modules the two trees share: a set never overlaps a copy
    # of itself, so keeping both copies would yield the same component twice
    U = filter!(overlap_components(Ms ∪ Mt)) do X
        (X in Ms || parent_node(s, X).kind != :prime) &&
        (X in Mt || parent_node(t, X).kind != :prime)
    end
    for x in V; push!(U, [x]); end
    R = Dict()
    for X in U
        S = parent_node(t, X)
        T = parent_node(s, X)
        union!(get!(()->Set{Int}(), R, (S, T)), X)
    end
    N = U ∪ map(sort!∘collect, values(R))
    N = filter!(X->length(X) > 1, N)
    T = Any[[] for x in V]
    for node in sort!(N, by=length, rev=true)
        an = Vector{Any}(node)
        for x in node
            push!(T[x], an)
        end
    end
    for x in V, i = 1:length(T[x])-1
        parent = T[x][i]
        child = T[x][i+1]
        child in parent && continue
        filter!(y->!(y in child), parent)
        push!(parent, child)
    end
    T = T[1][1]
    p = E[]
    record_vals(v::Vector) = foreach(record_vals, v)
    record_vals(x::E) = push!(p, x)
    record_vals(T)
    return p
end

function digraph_factorizing_permutation(
    G :: AbstractMatrix,
)
    n = checksquare(G)
    Gs = G .| G'
    Gd = G .& G'
    H = Gs .+ Gd
    ps = symgraph_factorizing_permutation(Gs)
    pd = symgraph_factorizing_permutation(Gd)
    s = StrongModuleTree(Gs, ps)
    t = StrongModuleTree(Gd, pd)
    p = intersect_permutation(1:n, s, t)
    h = StrongModuleTree(H, p)
    function sort_leaves!(h)
        for x in h.nodes
            x isa StrongModuleTree && sort_leaves!(x)
        end
        h.kind == :complete || return
        if h.edge == (1,1) # tournament node
            X = map(first_leaf, h.nodes)
            q = tournament_factorizing_permutation(G, X)
            o = Dict(x => i for (i, x) in enumerate(q))
            sort!(h.nodes, by=x->o[first_leaf(x)])
        else # 0/2-complete node
            d = Dict()
            c = (1:n)\leaves(h)
            for x in h.nodes
                i = first_leaf(x)
                push!(get!(d, G[c,i], []), x)
            end
            h.nodes .= [x for v in values(d) for x in v]
        end
    end
    sort_leaves!(h)
    return leaves(h)
end

function graph_factorizing_permutation(G::AbstractMatrix)
    a, b = extrema(G)
    0 ≤ a ≤ b ≤ 1 ||
        error("factoring multi-color two-structures not supported")
    issymmetric(G) ?
        symgraph_factorizing_permutation(G) :
        digraph_factorizing_permutation(G)
end

end # module

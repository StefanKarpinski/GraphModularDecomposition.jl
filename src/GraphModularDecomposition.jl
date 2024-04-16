module GraphModularDecomposition

export StrongModuleTree, strong_modules,
    graph_factorizing_permutation,
    digraph_factorizing_permutation,
    symgraph_factorizing_permutation,
    tournament_factorizing_permutation

using LinearAlgebra
using LinearAlgebra: checksquare
using SparseArrays

include("StrongModuleTrees.jl")
using .StrongModuleTrees

const \ = setdiff

function symgraph_factorizing_permutation(
    G :: AbstractMatrix,
    V :: Vector{Int} = collect(1:checksquare(G)),
)
    P = [V]
    center::Int = 0
    pivots::Vector{Vector{Int}} = []
    modules::Vector{Vector{Int}} = []
    first_pivot = Dict{Vector{Int},Int}()

    N_adj(x, X=V) = [y for y in X if y != x && G[x,y] != 0]
    N_non(x, X=V) = [y for y in X if y != x && G[x,y] == 0]

    smaller_larger(A, B) = length(A) <= length(B) ? (A, B) : (B, A)

    function refine!(P, S, x)
        i, between = 0, false
        while (i += 1) <= length(P)
            X = P[i]
            if center in X || x in X
                between = !between
                continue
            end
            Xₐ = X ∩ S
            isempty(Xₐ) && continue
            Xₙ = X \ Xₐ
            isempty(Xₙ) && continue
            P[i] = Xₙ
            insert!(P, i + between, Xₐ)
            add_pivot(X, Xₐ, Xₙ)
            i += 1
        end
    end

    function add_pivot(X, Xₐ, Xₙ)
        if X in pivots
            push!(pivots, Xₐ)
        else
            S, L = smaller_larger(Xₐ, Xₙ)
            push!(pivots, S)
            i = findfirst(isequal(X), modules)
            if i !== nothing
                modules[i] = L
            else
                push!(modules, L)
            end
        end
    end

    function partition_refinement!(P)
        while init_partition!(P)
            while !isempty(pivots)
                E = pop!(pivots)
                for x in E
                    S = N_adj(x) \ E
                    refine!(P, S, x)
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
                A, N = N_adj(x, X), N_non(x, X)
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

# compute whether A ∩ B, A \ B and B \ A are all non-empty
# this assumes that A and B are both sorted
function overlap(
    A :: AbstractVector{T},
    B :: AbstractVector{T},
) where {T}
    A === B && return false
    A_and_B = A_not_B = B_not_A = false
    m, n = length(A), length(B)
    (m ≤ 0 || n ≤ 0) && return false
    i = j = 1
    x, y = A[i], B[j]
    while i ≤ m && j ≤ n
        if x == y
            A_and_B = true # x ∈ A ∩ B
            x = get(A, i += 1, x)
            y = get(B, j += 1, y)
        elseif x < y
            A_not_B = true # x ∈ A \ B
            x = get(A, i += 1, x)
        else # y < x
            B_not_A = true # y ∈ B \ A
            y = get(B, j += 1, y)
        end
        A_not_B |= i ≤ m && j > n
        B_not_A |= i > m && j ≤ n
        A_and_B & A_not_B & B_not_A && return true
    end
    return false
end

function overlap_components(
    s :: StrongModuleTree,
    t :: StrongModuleTree,
    M = strong_modules(s),
    N = strong_modules(t),
)
    # TODO: implement linear time algorithm from this thesis:
    ## Dalhaus 1998: "Parallel algorithms for hierarchical clustering and
    ## applications to split decomposition and parity graph recognition"
    ## There's also a 2000 paper by the same name, but the PDF is jumbled

    O = M ∪ N
    n = length(O)
    R = SparseMatrixCSC{Int}(I, n, n)
    for i = 1:n-1, j = i+1:n
        overlap(O[i], O[j]) || continue
        R[i,j] = R[j,i] = 1
    end
    while true
        R′ = min.(1, R^2)
        R′ == R && break
        R .= R′
    end
    for (i, j) in map(Tuple, findall(!iszero, R))
        i < j || continue
        O[i] = O[j] = sort!(O[i] ∪ O[j])
    end
    return unique(O)
end

function intersect_permutation(
    V :: AbstractVector{E},   # vertices
    s :: StrongModuleTree{E}, # 1st strong module tree
    t :: StrongModuleTree{E}, # 2nd strong module tree
) where E
    Ms = strong_modules(s)
    Mt = strong_modules(t)
    U = filter!(overlap_components(s, t, Ms, Mt)) do X
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
    ps = graph_factorizing_permutation(Gs)
    pd = graph_factorizing_permutation(Gd)
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

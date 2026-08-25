module StrongModuleTrees

export StrongModuleTree, strong_modules, module_of,
    leaves, leaf_count, first_leaf, last_leaf, rand_leaf,
    nodes, node_count, parent_node, is_cograph

using ..GraphModularDecomposition
using ..GraphModularDecomposition: neighbors, is_symmetric
using SparseArrays: SparseMatrixCSC, nonzeros, nzrange, rowvals
using LinearAlgebra: checksquare

## core type definition ##

struct StrongModuleTree{T} <: AbstractVector{T}
    kind::Symbol
    edge::Tuple
    nodes::Vector{Union{T,StrongModuleTree{T}}}
end

## finding cutters ##

# A vertex w cuts the pair (u, v) when it is related differently to the two of
# them — when (G[w,u], G[u,w]) and (G[w,v], G[v,w]) differ.
#
# There are two ways to find the first and last such w, and which one is right
# depends on the representation, so there is a method per representation below.
#
# Scanning the permutation looks at every place before the cutter whether or not
# anything is there, so it is Θ(n) for a pair with no cutter and Θ(n²) over the
# whole permutation — but on a graph with edges everywhere it usually stops at
# the first place it looks, and it reads the matrix directly.
#
# Merging looks only at the places holding a relation, which is Θ(deg u + deg v)
# per pair and Θ(n + m) over the permutation — but it has to list every vertex's
# relations first, and that build costs the size of the representation. On a
# dense matrix the build alone is Θ(n²), against a scan that is usually Θ(1) per
# pair, so paying it there is a bad trade by three orders of magnitude.
#
# Both are linear in the representation they are chosen for.

transposed(G::SparseMatrixCSC) = copy(G')

# For each place in the permutation, the places related to it, in increasing
# order, with the labels in each direction. Stored as one flat array per field
# with a start index per place.
struct Relations
    start :: Vector{Int}
    place :: Vector{Int}
    in    :: Vector{Int}
    out   :: Vector{Int}
end

Base.getindex(r::Relations, j::Integer) = r.start[j]:r.start[j+1]-1

function Relations(G::SparseMatrixCSC, p::Vector{Int}, pos::Vector{Int})
    n = length(p)
    Gt = transposed(G)
    # collect each place's related places, deduplicating the two directions
    seen = falses(n)
    lists = Vector{Vector{Int}}(undef, n)
    for j = 1:n
        u, places = p[j], Int[]
        for w in Iterators.flatten((neighbors(G, u), neighbors(Gt, u)))
            k = pos[w]
            seen[k] && continue
            seen[k] = true
            push!(places, k)
        end
        for k in places
            seen[k] = false
        end
        lists[j] = places
    end
    total = sum(length, lists)
    start = zeros(Int, n + 2)
    start[1] = 1
    for j = 1:n
        start[j+1] = start[j] + length(lists[j])
    end
    start[n+2] = start[n+1]
    place, inlab, outlab = Vector{Int}(undef, total), Vector{Int}(undef, total),
                           Vector{Int}(undef, total)
    for j = 1:n
        u, at = p[j], start[j]
        sort!(lists[j])
        for k in lists[j]
            w = p[k]
            place[at], inlab[at], outlab[at] = k, G[w,u], G[u,w]
            at += 1
        end
    end
    return Relations(start, place, inlab, outlab)
end

# The first place before `j` that cuts the pair at places j and j+1, or 0.
function lower_cutter(r::Relations, j::Int)
    a, b = r[j], r[j+1]
    ia, ib = first(a), first(b)
    while true
        pa = ia <= last(a) ? r.place[ia] : typemax(Int)
        pb = ib <= last(b) ? r.place[ib] : typemax(Int)
        k = min(pa, pb)
        k < j || return 0
        if pa == k && pb == k
            (r.in[ia] != r.in[ib] || r.out[ia] != r.out[ib]) && return k
            ia += 1; ib += 1
        elseif pa == k
            return k # related to one and not the other
        else
            return k
        end
    end
end

# The last place after `j+1` that cuts the pair at places j and j+1, or 0.
function upper_cutter(r::Relations, j::Int)
    a, b = r[j], r[j+1]
    ia, ib = last(a), last(b)
    while true
        pa = ia >= first(a) ? r.place[ia] : 0
        pb = ib >= first(b) ? r.place[ib] : 0
        k = max(pa, pb)
        k > j + 1 || return 0
        if pa == k && pb == k
            (r.in[ia] != r.in[ib] || r.out[ia] != r.out[ib]) && return k
            ia -= 1; ib -= 1
        else
            return k
        end
    end
end

# scanning: for a matrix that stores every pair, so that reading one is O(1)
function find_cutters!(op, cl, lc, uc, G::AbstractMatrix, p::Vector{Int})
    n = length(p)
    for j = 1:n-1
        for i = 1:j-1
            G[p[i],p[j]] == G[p[i],p[j+1]] &&
            G[p[j],p[i]] == G[p[j+1],p[i]] && continue
            op[i] += 1
            cl[j] += 1
            lc[j] = i
            break
        end
        j += 1
        for i = n:-1:j+1
            G[p[i],p[j-1]] == G[p[i],p[j]] &&
            G[p[j-1],p[i]] == G[p[j],p[i]] && continue
            op[j] += 1
            cl[i] += 1
            uc[j-1] = i
            break
        end
    end
end

# merging: for a matrix that stores only the relations there are
function find_cutters!(op, cl, lc, uc, G::SparseMatrixCSC, p::Vector{Int})
    n = length(p)
    pos = zeros(Int, n)
    for (k, v) in enumerate(p)
        pos[v] = k
    end
    relations = Relations(G, p, pos)
    for j = 1:n-1
        i = lower_cutter(relations, j)
        if i != 0
            op[i] += 1
            cl[j] += 1
            lc[j] = i
        end
        i = upper_cutter(relations, j)
        if i != 0
            op[j+1] += 1
            cl[i] += 1
            uc[j] = i
        end
    end
end

## deciding a node's kind by counting crossings ##

"""
    crossing_kinds(G, v, root)

For each node of the nested structure `root`, whether it is `:complete` and with
which label, or `:prime` — or `nothing` when this cannot be answered this way.

`classify_nodes` decides a node by comparing each of its children against each
of the others. That is cheap for a prime node, which stops at the first pair
that disagrees, and quadratic for a degenerate one, which never disagrees: a
star's root has n-1 children with no edges among them, so a graph with n-1 edges
costs Θ(n²).

When the matrix is symmetric the question is only whether every pair of children
is related the same way, and that can be counted instead of compared. A pair of
vertices crosses between two children at exactly one node — the lowest one
holding both of them — so a single pass over the stored entries gives every node
the number of crossings it has and whether they all carry one label. Against the
number of pairs that could cross, which the children's sizes give, that decides
the node: none crossing means complete on the label of a non-edge, all crossing
under one label means complete on that label, and anything else is prime.

Only worth it for a matrix that stores its edges: reading a dense one costs the
Θ(n²) this is trying to avoid, and there Θ(Σk²) is no worse than the input.
"""
function crossing_kinds(G::SparseMatrixCSC, v::AbstractVector{<:Integer}, root::Vector)
    n = length(v)
    (checksquare(G) == n && sort(v) == 1:n && is_symmetric(G)) || return nothing

    place = Vector{Int}(undef, n)
    for (i, x) in enumerate(v)
        place[x] = i
    end

    # index the nodes: where each one starts and ends among the leaves, its
    # parent, the sum of the squares of its children's sizes, and for each leaf
    # the node directly above it
    ids = IdDict{Any,Int}()
    lo, hi, parent, sumsq = Int[], Int[], Int[], Int[]
    holder = zeros(Int, n)
    function index!(t::Vector, above::Int, at::Int)
        push!(lo, at); push!(hi, 0); push!(parent, above); push!(sumsq, 0)
        i = length(lo)
        ids[t] = i
        for x in t
            if x isa Vector
                at = index!(x, i, at)
                j = ids[x]
                sumsq[i] += (hi[j] - lo[j] + 1)^2
            else
                holder[at] = i
                sumsq[i] += 1
                at += 1
            end
        end
        hi[i] = at - 1
        return at
    end
    index!(root, 0, 1)

    # every stored entry crosses at the lowest node holding both its ends
    crossings = zeros(Int, length(lo))
    label = fill(-2, length(lo)) # -2 nothing yet, -1 more than one label
    rows, vals = rowvals(G), nonzeros(G)
    for column = 1:n, k in nzrange(G, column)
        row = rows[k]
        row < column && !iszero(vals[k]) || continue # each pair once
        a, b = minmax(place[row], place[column])
        i = holder[a]
        while hi[i] < b
            i = parent[i]
        end
        crossings[i] += 1
        label[i] = label[i] == -2 || label[i] == vals[k] ? vals[k] : -1
    end

    kinds = IdDict{Any,Tuple{Symbol,Tuple}}()
    for (t, i) in ids
        possible = ((hi[i] - lo[i] + 1)^2 - sumsq[i]) ÷ 2
        kinds[t] =
            crossings[i] == 0 ? (:complete, (0, 0)) :
            label[i] != -1 && crossings[i] == possible ?
                (:complete, (label[i], label[i])) : (:prime, ())
    end
    return kinds
end

crossing_kinds(G::AbstractMatrix, v::AbstractVector, root::Vector) = nothing

"""
    module_of(G, p, h)

For each node of `h`, whether its leaves are a module of `G` — or `nothing` when
that cannot be answered this way.

`h` is the tree of a coarser structure than `G`: `digraph_factorizing_permutation`
builds it over `Gs + Gd`, which records whether a pair is joined but not which
way round, so a node of `h` can be a module there and not in `G`. Deciding that
by comparing every leaf of a node against every vertex outside it costs
O(|node| · |outside|), which on a deep tree is worse than everything else in the
pipeline put together.

McConnell & de Montgolfier 2005, lemma 12: a node is a module of `G` exactly
when no cutter, taken in `G`, of a consecutive pair inside it lies outside it.
The cutters of consecutive pairs come from the same merge used to build the tree,
and fold up the tree in one pass, so this is O(1) a node afterwards.
"""
function module_of(G::SparseMatrixCSC, p::Vector{Int}, h::StrongModuleTree)
    n = length(p)
    checksquare(G) == n && sort(p) == 1:n || return nothing
    pos = zeros(Int, n)
    for (k, v) in enumerate(p)
        pos[v] = k
    end
    relations = Relations(G, p, pos)
    locut, hicut = fill(typemax(Int), n), zeros(Int, n)
    for j = 1:n-1
        i = lower_cutter(relations, j); i != 0 && (locut[j] = i)
        i = upper_cutter(relations, j); i != 0 && (hicut[j] = i)
    end
    answer = IdDict{Any,Bool}()
    # returns where the node ends, and the extreme cutters of the pairs inside it
    function fold(x, at::Int)
        x isa StrongModuleTree || return at + 1, typemax(Int), 0
        first = at
        low, high = typemax(Int), 0
        for (i, c) in enumerate(x.nodes)
            i > 1 && (low = min(low, locut[at-1]); high = max(high, hicut[at-1]))
            at, l, h = fold(c, at)
            low, high = min(low, l), max(high, h)
        end
        answer[x] = low >= first && high <= at - 1
        return at, low, high
    end
    fold(h, 1)
    return answer
end

module_of(G::AbstractMatrix, p::Vector{Int}, h::StrongModuleTree) = nothing

## constructing StrongModuleTrees ##

function StrongModuleTree(
    G::AbstractMatrix, # the graph/digraph/2-structure as a matrix
    p::Vector{Int} = graph_factorizing_permutation(G),
)
    # initialize data structures
    n = length(p)
    op = zeros(Int,n); op[1] = 1
    cl = zeros(Int,n); cl[n] = 1
    lc = collect(1:n-1)
    uc = collect(2:n)

    # count open and close parens in fracture tree
    # find lower and upper cutters for node pairs
    find_cutters!(op, cl, lc, uc, G, p)

    # remove non-module "dummy" nodes
    #
    # a node spanning i:j is a module iff no cutter of a pair inside it lies
    # outside it. rescanning lc and uc across the whole span to find that out
    # costs O(span) per node; instead every open paren carries the extremes of
    # everything seen since it opened, and passes them to its parent on the way
    # out, which is O(1) per node
    let s = Int[], lo = Int[], hi = Int[]
        for j = 1:n
            # the cutters of the pair (j-1, j) belong to every node still open
            if j > 1 && !isempty(s)
                lo[end] = min(lo[end], lc[j-1])
                hi[end] = max(hi[end], uc[j-1])
            end
            for _ = 1:op[j]
                push!(s, j)
                push!(lo, typemax(Int))
                push!(hi, typemin(Int))
            end
            for _ = 1:cl[j]
                i, l, u = pop!(s), pop!(lo), pop!(hi)
                if !isempty(s)
                    lo[end] = min(lo[end], l)
                    hi[end] = max(hi[end], u)
                end
                if i < j
                    i <= l && u <= j && continue
                end
                op[i] -= 1
                cl[j] -= 1
            end
        end
    end

    # create nodes for consecutive twins
    let s = Int[], t = Int[]
        l = 1
        for k = 1:n
            for _ = 1:op[k]+1
                push!(s, k) # matching node stack
                push!(t, l) # matching twin stack
                l = k
            end
            for c = cl[k]:-1:0
                i = pop!(t)
                j = pop!(s)
                l = i # continue twin chain by default
                i < j || continue
                if i <= lc[j-1] < uc[j-1] <= k
                    # this node and prev are twins
                    if c > 0
                        # not last parens ∴ last twin
                        op[i] += 1
                        cl[k] += 1
                        l = k + 1
                    end
                else # this node and prev aren't twins
                    if i < j-1
                        op[i] += 1
                        cl[j-1] += 1
                    end
                    l = j # this node starts new chain
                end
            end
        end
    end

    # remove singleton "dummy" nodes
    let s = Int[]
        for j = 1:n
            for _ = 1:op[j]; push!(s, j); end
            i′ = 0
            for _ = 1:cl[j]
                i = pop!(s)
                if i == i′
                    op[i] -= 1
                    cl[j] -= 1
                end
                i′ = i
            end
        end
    end
    op[1] -= 1
    cl[n] -= 1

    # construct and normalize the tree
    return StrongModuleTree(G, p, op, cl)
end

function StrongModuleTree(
    G::AbstractMatrix,    # graph/digraph/2-structure
    v::AbstractVector{T}, # modular permutation of vertices
    op::Vector{Int},      # open parens per vertex
    cl::Vector{Int},      # close parens per vertex
) where {T <: Any}
    # continues from end of StrongModuleTree(G, p)

    # returns the classified node together with its first leaf. the leaf comes
    # back up rather than being searched for again because first_leaf walks to
    # the bottom of the subtree: calling it inside the comparison loop below
    # costs O(depth) every time, which is quadratic on a deep tree
    function classify_nodes(t::Vector)
        n = length(t)
        nodes = Vector{Union{T,StrongModuleTree{T}}}(undef, n)
        leaf = Vector{T}(undef, n)
        for i = 1:n
            x = t[i]
            nodes[i], leaf[i] = x isa Vector ? classify_nodes(x) : (x, x)
        end
        # counted rather than compared, where counting applies
        if kinds !== nothing
            kind, edge = kinds[t]
            return StrongModuleTree{T}(kind, edge, nodes), leaf[1]
        end
        counts = zeros(Int, n)
        x, y = leaf[1], leaf[2]
        edge = (G[y,x], G[x,y])
        local a, b
        for i = 1:n, j = 1:n
            i == j && continue
            x, y = leaf[i], leaf[j]
            a, b = G[y,x], G[x,y]
            if edge == (a, b)
                counts[i] += 1
            elseif edge == (b, a)
                counts[j] += 1
            else
                break
            end
        end
        sort!(counts)
        kind = a == b && all(c -> c == n-1, counts) ? :complete :
            all(d -> d == 2, diff(counts)) ? :linear : :prime
        edge[1] <= edge[2] || (edge = reverse(edge))
        kind == :prime && (edge = ())
        return StrongModuleTree{T}(kind, edge, nodes), leaf[1]
    end

    function delete_weak_modules!(t::StrongModuleTree)
        i = 0
        while (i += 1) <= length(t)
            x = t[i]
            x isa StrongModuleTree || continue
            delete_weak_modules!(x)
            t.kind == x.kind != :prime && t.edge == x.edge || continue
            splice!(t.nodes, i, x.nodes)
            i += length(x)
        end
    end

    s = Any[[]]
    for (j, x) = enumerate(v)
        for _ = 1:op[j]
            t = []
            push!(s[end], t)
            push!(s, t)
        end
        push!(s[end], x)
        for _ = 1:cl[j]
            pop!(s)
        end
    end
    kinds = crossing_kinds(G, v, s[end])
    t, _ = classify_nodes(s[end])
    delete_weak_modules!(t)
    return t
end

## displaying StrongModuleTrees ##

function edge_string(t::StrongModuleTree, post::String="")
    t.kind == :prime    ? "" :
    t.kind == :complete ? "$(t.edge[1])$post" :
                          "$(join(t.edge,"/"))$post"
end

kind_string(t::StrongModuleTree) = "$(edge_string(t,"-"))$(t.kind)"
Base.summary(t::StrongModuleTree) = "$(length(t))-node $(kind_string(t)) $(typeof(t))"

function Base.show(io::IO, t::StrongModuleTree)
    if get(io, :compact, false)
        print(io,
            edge_string(t,"-"), t.kind, " ",
            node_count(t), "-node (",
            leaf_count(t), "-leaf) module: ",
            first_leaf(t)
        )
    else
        parens = t.kind == :prime ? "{}" : t.kind == :linear ? "[]" : "()"
        print(io, parens[1])
        for (i, x) in enumerate(t)
            print(io, x)
            i < length(t) && print(io, " ")
        end
        print(io, parens[2])
    end
end

## querying properties of StrongModuleTrees ##

node_count(t::StrongModuleTree) = length(t)
node_count(x::Any) = 1

leaf_count(t::StrongModuleTree) = sum(leaf_count, t.nodes)
leaf_count(v::Vector) = sum(leaf_count, v)
leaf_count(x::Any) = 1

first_leaf(t::StrongModuleTree) = first_leaf(first(t.nodes))
first_leaf(v::Vector) = first_leaf(first(v))
first_leaf(x::Any) = x

last_leaf(t::StrongModuleTree) = last_leaf(last(t.nodes))
last_leaf(v::Vector) = last_leaf(last(v))
last_leaf(x::Any) = x

rand_leaf(t::StrongModuleTree) = rand_leaf(rand(t.nodes))
rand_leaf(v::Vector) = rand_leaf(last(v))
rand_leaf(x::Any) = x

# Collected into one buffer rather than by concatenating each child's leaves
# into the parent's. Concatenating rebuilds every subtree once per level above
# it, so a single call costs O(size · depth) and asking every node of a deep
# tree for its leaves -- which is what strong_modules does -- costs far more
# than the leaves are worth. Theorem 8 of McConnell & de Montgolfier 2005 bounds
# the leaves of all the strong modules together by 2m + 3n, so gathering them
# should cost the size of the graph and no more.
function leaves!(L::Vector, t::StrongModuleTree)
    for x in t.nodes
        x isa StrongModuleTree ? leaves!(L, x) : push!(L, x)
    end
    return L
end

leaves(t::StrongModuleTree{T}) where {T} = leaves!(T[], t)
leaves(x::Any) = [x]

function nodes!(v::Vector{StrongModuleTree{T}}, t::StrongModuleTree{T}) where T
    for x in t.nodes
        x isa StrongModuleTree || continue
        nodes!(v, x)
    end
    push!(v, t)
    return v
end
nodes(t::StrongModuleTree) = nodes!(typeof(t)[], t)

strong_modules(t::StrongModuleTree) = map(sort!∘leaves, nodes(t))

"""
    TreeIndex(t)

An index over the nodes of `t`: where each leaf falls in the tree's leaf order,
which range of that order each node spans, and each node's parent.

`parent_node` answers "which node sits just above this set of leaves" by
materializing the leaf set of every child it passes, which costs a subtree walk
per step. The same question is a lowest-common-ancestor query: the node wanted
is the lowest one whose leaf range covers the set, and a set's range is fixed by
its lowest and highest leaf alone. Building the index costs one pass over the
tree, after which a query is O(|S| + depth) and touches no leaf sets at all.
"""
struct TreeIndex{T}
    leaf   :: Dict{T,Int}                   # leaf -> its place in the leaf order
    node   :: Vector{StrongModuleTree{T}}
    parent :: Vector{Int}                   # node -> its parent, 0 for the root
    first  :: Vector{Int}                   # node -> its first leaf's place
    last   :: Vector{Int}                   # node -> its last leaf's place
    holder :: Vector{Int}                   # leaf place -> the node just above it
end

function TreeIndex(t::StrongModuleTree{T}) where {T}
    ix = TreeIndex{T}(Dict{T,Int}(), StrongModuleTree{T}[], Int[], Int[], Int[], Int[])
    function visit(x::StrongModuleTree{T}, parent::Int)
        push!(ix.node, x); push!(ix.parent, parent)
        push!(ix.first, length(ix.leaf) + 1); push!(ix.last, 0)
        i = length(ix.node)
        for y in x.nodes
            if y isa StrongModuleTree
                visit(y, i)
            else
                ix.leaf[y] = length(ix.leaf) + 1
                push!(ix.holder, i)
            end
        end
        ix.last[i] = length(ix.leaf)
        return i
    end
    visit(t, 0)
    return ix
end

"""
    containing_node(ix, S) -> (node, exact)

The lowest node of the indexed tree whose leaves contain every element of `S`,
and whether its leaves are exactly `S` — which is what asking whether `S` is one
of the tree's strong modules amounts to.
"""
function containing_node(ix::TreeIndex, S)
    lo = hi = 0
    for x in S
        p = ix.leaf[x]
        lo == 0 && (lo = hi = p)
        p < lo && (lo = p)
        p > hi && (hi = p)
    end
    lo == 0 && return (1, false) # nothing to contain: the root, and not exactly
    i = ix.holder[lo]
    while (ix.first[i] > lo || ix.last[i] < hi) && ix.parent[i] != 0
        i = ix.parent[i]
    end
    return (i, ix.last[i] - ix.first[i] + 1 == length(S))
end

"""
    parent_index(ix, S)

The node [`parent_node`](@ref) would return: the lowest node whose leaves
*strictly* contain `S`, so a node whose leaves are exactly `S` yields its parent.
"""
function parent_index(ix::TreeIndex, S)
    i, exact = containing_node(ix, S)
    exact && ix.parent[i] != 0 ? ix.parent[i] : i
end

function parent_node(t::StrongModuleTree, S::Vector)
    for x in t.nodes
        x isa StrongModuleTree || continue
        L = leaves(x)
        isempty(S ∩ L) && continue
        S ⊆ L && length(S) < length(L) && return parent_node(x, S)
        break
    end
    return t
end

## StrongModuleTree as a collection ##

Base.size(t::StrongModuleTree) = (length(t),)
Base.length(t::StrongModuleTree) = length(t.nodes)
Base.eltype(t::StrongModuleTree) = eltype(t.nodes)
Base.getindex(t::StrongModuleTree, i::Int) = t.nodes[i]

Base.getindex(v::Vector{T}, t::StrongModuleTree) where {T} =
    StrongModuleTree{T}(t.kind, t.edge, map(x->v[x], t.nodes))

function Base.sort!(t::StrongModuleTree; lt=isless, by=first_leaf, rev::Bool=false)
    for x in t.nodes
        x isa StrongModuleTree || continue
        sort!(x, lt=lt, by=by, rev=rev)
    end
    sort!(t.nodes, lt=lt, by=by, rev=rev)
    return t
end

## predicates on StrongModuleTrees

is_cograph(t::StrongModuleTree) = t.kind != :prime && all(is_cograph, t.nodes)
is_cograph(x::Any) = true

end # module

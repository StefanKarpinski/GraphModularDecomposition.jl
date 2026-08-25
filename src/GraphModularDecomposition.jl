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

"""
    is_symmetric(G)

Whether `G` equals its own transpose.

`LinearAlgebra.issymmetric` is not to be trusted here: on a sparse matrix it can
report symmetry that is not there, because while looking for the mirror of a
stored entry it can read past the end of the column it is searching and find an
entry belonging to the next one. JuliaSparse/SparseArrays.jl#748, and #636
before it, which was fixed only for the case of an entirely empty column.
Getting this wrong runs the undirected algorithm on a digraph and returns a
permutation that is not a factorizing permutation at all.

Comparing against a materialized transpose is correct and costs the size of the
representation: O(n + nnz) for a sparse matrix, which is in fact quicker than
`issymmetric`, and O(n²) for a dense one, where `issymmetric` is both correct
and allocation-free and so is left in place.
"""
is_symmetric(G::AbstractMatrix) = issymmetric(G)
is_symmetric(G::SparseMatrixCSC) = G == permutedims(G)

include("StrongModuleTrees.jl")
using .StrongModuleTrees
using .StrongModuleTrees: TreeIndex, containing_node, parent_index


function symgraph_factorizing_permutation(
    G :: AbstractMatrix,
    V :: Vector{Int} = collect(1:checksquare(G)),
)
    n, m = checksquare(G), length(V)
    m <= 1 && return copy(V)

    # An ordered partition of V whose parts are consecutive ranges of `elems`.
    # Splitting a part rearranges only its own range, so no other part moves and
    # no part's position shifts; two parts therefore compare by position in O(1),
    # and the part holding a vertex is one lookup away. That is what lets refine!
    # reach the parts a pivot set meets through the pivot set's own elements,
    # instead of walking the whole partition looking for them.
    elems = copy(V)
    pos, part = zeros(Int, n), zeros(Int, n)
    for (i, v) in enumerate(elems)
        pos[v], part[v] = i, 1
    end
    lo, hi = [1], [m]
    hits = [Int[]]     # scratch per part: which of the pivot set's vertices it holds
    first_pivot = [0]  # part -> the vertex it should be refined by first, 0 if none
    piv_pos = [0]      # part -> its live index in `pivots`, 0 if absent
    mod_pos = [0]      # part -> its live index in `modules`, 0 if absent
    nparts = 1

    function new_part!(l, h)
        push!(lo, l); push!(hi, h)
        push!(hits, Int[]); push!(first_pivot, 0)
        push!(piv_pos, 0); push!(mod_pos, 0)
        nparts += 1
        return length(lo)
    end

    part_size(p) = hi[p] - lo[p] + 1

    # `pivots` is Habib & Paul's L, every vertex of whose parts is used as a
    # pivot, and `modules` is their K, which contributes one vertex per part.
    # Both hold part numbers, except that a negative entry in `pivots` denotes
    # the one-vertex set {-e}, which K's rule produces and which is not a part.
    # An entry is removed by clearing its position index and leaving the entry
    # itself to be skipped on the way out, so removing from the middle is O(1).
    pivots, modules = Int[], Int[]
    npivots, nmodules, mread = 0, 0, 1

    function push_pivot!(p)
        push!(pivots, p)
        p > 0 && (piv_pos[p] = length(pivots))
        npivots += 1
    end

    function pop_pivot!()
        while !isempty(pivots)
            k = length(pivots)
            e = pop!(pivots)
            if e < 0
                npivots -= 1
                return e
            elseif piv_pos[e] == k # a live entry rather than a superseded one
                piv_pos[e] = 0
                npivots -= 1
                return e
            end
        end
        return 0
    end

    function drop_pivot!(p)
        piv_pos[p] == 0 && return
        piv_pos[p] = 0
        npivots -= 1
    end

    function push_module!(p)
        push!(modules, p)
        mod_pos[p] = length(modules)
        nmodules += 1
    end

    function replace_module!(old, new)
        j = mod_pos[old]
        mod_pos[old] = 0
        modules[j] = new
        mod_pos[new] = j
    end

    function popfirst_module!()
        while mread <= length(modules)
            p = modules[mread]
            mread += 1
            if mod_pos[p] == mread - 1
                mod_pos[p] = 0
                nmodules -= 1
                return p
            end
        end
        return 0
    end

    center::Int = 0
    in_splitter = falses(n)   # the neighbourhood being split by, in init_partition!
    in_pivot_part = falses(n) # the part the current pivot vertices came from
    in_V = falses(n)          # V may be a subset of the vertices of G
    for y in V
        in_V[y] = true
    end

    # Move `ys`, all of which lie in part p, to the front or the back of p's
    # range and split them off as a part of their own. Only the elements of `ys`
    # move, so this costs O(|ys|) and not O(|p|); the price is that they end up
    # in the order `ys` had rather than the order p had, which changes which
    # factorizing permutation comes out but not which trees it factors.
    # Returns the two parts, whose vertices did and did not belong to `ys`.
    function split_part!(p, ys, last)
        k, sz = length(ys), part_size(p)
        (k == 0 || k == sz) && return (0, 0)
        cursor = last ? hi[p] : lo[p]
        for y in ys
            j, u = pos[y], elems[cursor]
            elems[j], pos[u] = u, j
            elems[cursor], pos[y] = y, cursor
            cursor += last ? -1 : 1
        end
        ylo, yhi = last ? (hi[p] - k + 1, hi[p]) : (lo[p], lo[p] + k - 1)
        rlo, rhi = last ? (lo[p], hi[p] - k) : (lo[p] + k, hi[p])
        # the new part number goes to the smaller side, so that relabelling its
        # vertices costs O(min(k, sz - k)) rather than O(sz)
        if k <= sz - k
            q = new_part!(ylo, yhi)
            for j = ylo:yhi
                part[elems[j]] = q
            end
            lo[p], hi[p] = rlo, rhi
            return (q, p)
        else
            q = new_part!(rlo, rhi)
            for j = rlo:rhi
                part[elems[j]] = q
            end
            lo[p], hi[p] = ylo, yhi
            return (p, q)
        end
    end

    # Habib & Paul 2009, algorithm 2: a part of L splits into two parts of L,
    # and a part outside L contributes its smaller half to L and leaves its
    # larger half to K.
    function record_split!(p, a, b)
        if piv_pos[p] != 0
            drop_pivot!(p)
            push_pivot!(a)
            push_pivot!(b)
        else
            small, large = part_size(a) <= part_size(b) ? (a, b) : (b, a)
            push_pivot!(small)
            if mod_pos[p] != 0
                replace_module!(p, large)
            else
                push_module!(large)
            end
        end
    end

    touched = Int[]

    function refine!(S, x)
        # the parts holding the center and the pivot are never split, and their
        # positions are what everything else is placed relative to
        pc, px = part[center], part[x]
        empty!(touched)
        for y in S
            p = part[y]
            (p == pc || p == px) && continue
            isempty(hits[p]) && push!(touched, p)
            push!(hits[p], y)
        end
        for p in touched
            # a part lying between the center and the pivot keeps its neighbours
            # of x on the far side from the center; every other part keeps them
            # on the near side
            between = lo[pc] < lo[p] < lo[px] || lo[px] < lo[p] < lo[pc]
            a, b = split_part!(p, hits[p], between)
            empty!(hits[p])
            a == 0 && continue
            first_pivot[p] = 0 # p is not the part it was
            record_split!(p, a, b)
        end
    end

    # positions before this one hold singleton parts, and parts only shrink, so
    # the search for a part still worth refining never has to start over
    unrefined = 1

    function next_part()
        while unrefined <= m
            p = part[elems[unrefined]]
            part_size(p) > 1 && return p
            unrefined = hi[p] + 1
        end
        return 0
    end

    function init_partition!()
        nparts == m && return false
        if nmodules == 0
            p = next_part()
            p == 0 && return false
            x = first_pivot[p] != 0 ? first_pivot[p] : elems[lo[p]]
            first_pivot[p] = 0 # p is about to stop being the part it was
            # Split p into the neighbours of x, then x, then the rest, by moving
            # only x and its neighbours to the front of p's range. Whatever is
            # left is already in place and becomes the third part, so this costs
            # O(deg x) rather than O(|p|).
            #
            # That matters because a part whose vertices are all twins can never
            # be refined, so this loop is what takes it apart, one vertex per
            # call. A star's leaves are such a part, and scanning it each time
            # made a graph with n-1 edges cost Θ(n²).
            l, h = lo[p], hi[p]
            cursor = l
            for y in neighbors(G, x)
                (y != x && part[y] == p) || continue
                j = pos[y]
                u = elems[cursor]
                elems[j], pos[u] = u, j
                elems[cursor], pos[y] = y, cursor
                cursor += 1
            end
            j = pos[x] # x sits at or after the neighbours, never among them
            u = elems[cursor]
            elems[j], pos[u] = u, j
            elems[cursor], pos[x] = x, cursor
            xj = cursor
            na, nn = xj - l, h - xj
            # the old number stays with the larger side, so only the smaller
            # side's vertices are relabelled
            px = new_part!(xj, xj)
            part[x] = px
            pa = pn = 0
            if na > 0 && nn > 0
                if na >= nn
                    pa = p; lo[p], hi[p] = l, xj - 1
                    pn = new_part!(xj + 1, h)
                    for j = xj+1:h; part[elems[j]] = pn; end
                else
                    pn = p; lo[p], hi[p] = xj + 1, h
                    pa = new_part!(l, xj - 1)
                    for j = l:xj-1; part[elems[j]] = pa; end
                end
            elseif na > 0
                pa = p; lo[p], hi[p] = l, xj - 1
            else
                pn = p; lo[p], hi[p] = xj + 1, h
            end
            drop_pivot!(p)
            mod_pos[p] != 0 && (mod_pos[p] = 0; nmodules -= 1)
            center = x
            small, large = na <= nn ? (pa, pn) : (pn, pa)
            small != 0 && push_pivot!(small)
            large != 0 && push_module!(large)
        else
            p = popfirst_module!()
            p == 0 && return false
            x = elems[lo[p]]
            push_pivot!(-x) # the one-vertex set {x}
            first_pivot[p] = x
        end
        return true
    end

    while init_partition!()
        while npivots > 0
            e = pop_pivot!()
            e == 0 && break
            E = e < 0 ? [-e] : elems[lo[e]:hi[e]]
            for y in E
                in_pivot_part[y] = true
            end
            for x in E
                # the neighbourhood of x minus the part it came from, built in
                # one pass: O(deg x) given a representation that can list a
                # vertex's neighbours, O(n) given one that cannot
                S = [y for y in neighbors(G, x) if in_V[y] && !in_pivot_part[y]]
                refine!(S, x)
            end
            for y in E
                in_pivot_part[y] = false
            end
        end
    end
    return elems # every part is a singleton, so this is the permutation
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
    # Index both trees once. Everything below asks, of a great many sets, which
    # node sits above the set and whether the set is exactly one of the tree's
    # modules; unindexed, the first walks the tree materializing the leaf set of
    # every child it passes and the second scans the whole list of modules.
    si, ti = TreeIndex(s), TreeIndex(t)
    # Ms ∪ Mt drops the modules the two trees share: a set never overlaps a copy
    # of itself, so keeping both copies would yield the same component twice
    U = filter!(overlap_components(Ms ∪ Mt)) do X
        # `X in Ms` is `containing_node` reporting an exact match, and when it
        # reports none the node it found is what parent_node would have returned
        i, exact = containing_node(si, X)
        (exact || si.node[i].kind != :prime) || return false
        j, exact = containing_node(ti, X)
        return exact || ti.node[j].kind != :prime
    end
    for x in V; push!(U, [x]); end
    # keyed by node numbers rather than by the nodes themselves, which hash
    # their whole subtree
    R = Dict{Tuple{Int,Int},Set{Int}}()
    for X in U
        key = (parent_index(ti, X), parent_index(si, X))
        union!(get!(()->Set{Int}(), R, key), X)
    end
    N = U ∪ map(sort!∘collect, values(R))
    N = filter!(X->length(X) > 1, N)
    # N is laminar, so it is the node set of a tree, and sorting by decreasing
    # size makes the sets holding a given leaf that leaf's ancestors, outermost
    # first. That gives every set its depth in one pass over the sets. Nesting
    # them by searching and filtering instead costs O(|parent|·|child|) for each
    # parent and child, and pays it once per leaf of the child.
    sort!(N, by=length, rev=true)
    chain = [Int[] for _ in V] # leaf -> its ancestors in N, outermost first
    for (i, X) in enumerate(N), x in X
        push!(chain[x], i)
    end
    depth = zeros(Int, length(N))
    for c in chain, (k, i) in enumerate(c)
        depth[i] = k
    end
    # Walk the tree, taking each node's own leaves before any of its children.
    # Grouping them is not cosmetic: a module can fail to be a node of N and
    # still have to come out contiguous, and the only thing keeping such a
    # module together is that the leaves belonging to a node directly are not
    # separated by that node's children.
    p, entered = E[], falses(length(N))
    function emit(i)
        for x in N[i]
            length(chain[x]) == depth[i] && push!(p, x) # i holds x directly
        end
        for x in N[i]
            length(chain[x]) == depth[i] && continue
            c = chain[x][depth[i]+1]
            entered[c] && continue
            entered[c] = true
            emit(c)
        end
    end
    emit(first(chain[first(V)]))
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
            # No edges run between this node's children, or all of them are
            # reciprocal, so which unions of children are modules of G is
            # decided entirely outside the node: a union is a module exactly
            # when its children relate identically to every outside vertex, in
            # both directions. Group them by that, so each such union comes out
            # contiguous.
            #
            # This is the relation R_X of McConnell & de Montgolfier 2005,
            # "Linear-time modular decomposition of directed graphs", lemma 10
            # and step 4 of their algorithm 3: S ~ S' when N⁺(S) and N⁺(S')
            # agree outside the node and N⁻(S) and N⁻(S') do too, taken over
            # the children that are modules of G. Their labels for H are not
            # these -- theirs number a reciprocal pair 1 and a one-way pair 2,
            # H = Gs + Gd numbers them the other way round -- so their
            # "0-complete and 1-complete" is this branch and their "2-complete"
            # is the tournament branch above, which takes arbitrary
            # representatives on the strength of their lemma 11. Representatives
            # are sound there and not here, which is what went wrong.
            #
            # Both halves of this matter. Looking only at G[v,x], as this used
            # to, misses unions that are separated by the other direction. And
            # taking the relation from one leaf of a child misrepresents the
            # child: being a module of H does not make it a module of G -- H
            # records only whether a pair is joined, not which way -- so a
            # child's own leaves can disagree, and then no leaf speaks for it.
            # Such a child joins no union and is set aside.
            outside = (1:n)\leaves(h)
            # this path only ever sees 0/1 entries -- G .| G' and G .& G' above
            # are bitwise and mean nothing otherwise -- so one number per
            # outside vertex distinguishes all four ways a pair can be related
            relation(x, v) = G[v,x] + 2G[x,v]
            index = Dict{Vector{Int},Int}()
            groups, aside = Vector{Any}[], Any[]
            for x in h.nodes
                L = leaves(x)
                i = first(L)
                signature = [relation(i, v) for v in outside]
                if all(relation(y, v) == relation(i, v) for y in L, v in outside)
                    k = get!(index, signature, length(groups) + 1)
                    k > length(groups) ? push!(groups, Any[x]) : push!(groups[k], x)
                else
                    push!(aside, x)
                end
            end
            h.nodes .= vcat([x for g in groups for x in g], aside)
        end
    end
    sort_leaves!(h)
    return leaves(h)
end

function graph_factorizing_permutation(G::AbstractMatrix)
    a, b = extrema(G)
    0 ≤ a ≤ b ≤ 1 ||
        error("factoring multi-color two-structures not supported")
    is_symmetric(G) ?
        symgraph_factorizing_permutation(G) :
        digraph_factorizing_permutation(G)
end

end # module

"""
Overlap classes of a family of sets, in linear time.

Two sets `A` and `B` *overlap* if `A ∩ B`, `A \\ B` and `B \\ A` are all
non-empty. The *overlap graph* of a family `F` has the sets of `F` as its
vertices and joins two sets that overlap; its connected components are the
*overlap classes* of `F`.

The overlap graph can have Θ(|F|²) edges — consider `{x,y₁}, {x,y₂}, …` — so it
cannot be built explicitly in linear time. Dahlhaus' algorithm sidesteps this by
building a *different* graph `D` on the same vertex set which has at most
Σ|X| edges and, although it is in general not a subgraph of the overlap graph,
provably has the same connected components.

Implemented from:

  * Dahlhaus 1998/2000: "Parallel algorithms for hierarchical clustering and
    applications to split decomposition and parity graph recognition",
    J. Algorithms 36(2):205-240, section 4.
  * Charbit, Habib, Limouzy, de Montgolfier, Raffinot & Rao 2007:
    "A Note on Computing Set Overlap Classes", arXiv:0711.4573.

The second paper restates Dahlhaus' section 4 much more legibly and, crucially,
replaces its offline-lowest-common-ancestor subroutine — linear only in theory —
with a second pass of ordered partition refinement, which is what makes the
Θ(n + Σ|X|) bound achievable in practice. This module follows that variant.
"""
module OverlapComponents

export overlap, overlap_classes, overlap_components

## the overlap predicate ##

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

## linear-time helpers ##

# stable counting sort of `items` by `key[item] ∈ 0:kmax`, descending
function sortperm_desc(items, key::Vector{Int}, kmax::Int)
    count = zeros(Int, kmax + 1)
    for i in items
        count[key[i]+1] += 1
    end
    next, acc = Vector{Int}(undef, kmax + 1), 1
    for k = kmax:-1:0
        next[k+1] = acc
        acc += count[k+1]
    end
    order = Vector{Int}(undef, length(items))
    for i in items
        order[next[key[i]+1]] = i
        next[key[i]+1] += 1
    end
    return order
end

# union-find over 1:m with path halving and union by size
function uf_find(parent::Vector{Int}, x::Int)
    while parent[x] != x
        parent[x] = parent[parent[x]]
        x = parent[x]
    end
    return x
end

function uf_union!(parent::Vector{Int}, size::Vector{Int}, x::Int, y::Int)
    x, y = uf_find(parent, x), uf_find(parent, y)
    x == y && return false
    size[x] < size[y] && ((x, y) = (y, x))
    parent[y] = x
    size[x] += size[y]
    return true
end

## ordered partition refinement ##

# An ordered partition of the ground set, represented as in Charbit et al.
# figure 4: a permutation `elems` of the ground set, its inverse `pos`, and a
# collection of parts, each a contiguous range `lo[p]:hi[p]` of `elems`.
struct Partition
    elems :: Vector{Int} # ground set in partition order
    pos   :: Vector{Int} # pos[v] = index of v in elems
    part  :: Vector{Int} # part[v] = id of the part containing v
    lo    :: Vector{Int} # lo[p]:hi[p] = extent of part p within elems
    hi    :: Vector{Int}
    count :: Vector{Int} # scratch: pivot elements seen so far in part p
end

function Partition(elems::Vector{Int})
    n = length(elems)
    pos = Vector{Int}(undef, n)
    for (i, v) in enumerate(elems)
        pos[v] = i
    end
    Partition(elems, pos, ones(Int, n), [1], [n], [0])
end

Partition(n::Integer) = Partition(collect(1:n))

# split the last `k` elements of part `p` off as a new part, returning its id
function split_part!(P::Partition, p::Int, k::Int)
    push!(P.lo, P.hi[p] - k + 1)
    push!(P.hi, P.hi[p])
    push!(P.count, 0)
    q = length(P.lo)
    for j = P.lo[q]:P.hi[q]
        P.part[P.elems[j]] = q
    end
    P.hi[p] = P.lo[q] - 1
    return q
end

## step 1: lexicographically order the ground set ##

# Refine the trivial partition of the ground set by each set of `F` in `order`,
# splitting every part `C` it meets into `C \ X` followed by `C ∩ X`.
#
# Read `F` as a boolean matrix whose rows are the sets in `order` and whose
# columns are ground set elements. This is then an MSD radix sort of the
# columns: the resulting partition order is the lexicographic order of the
# columns with the *first* row of `order` most significant — exactly the column
# order Dahlhaus' algorithm needs (Charbit et al. lemma 7).
function lex_order(F, order, n::Int)
    P = Partition(n)
    met = Int[]
    for i in order
        empty!(met)
        for v in F[i]
            p = P.part[v]
            P.count[p] == 0 && push!(met, p)
            P.count[p] += 1
            # swap v down into the next free slot at the tail of its part; v is
            # never already past that slot, since the slots beyond it hold the
            # elements of F[i] moved so far and each is moved exactly once
            j = P.hi[p] - P.count[p] + 1
            k, u = P.pos[v], P.elems[j]
            P.elems[k], P.pos[u] = u, k
            P.elems[j], P.pos[v] = v, j
        end
        for p in met
            k, P.count[p] = P.count[p], 0
            k < P.hi[p] - P.lo[p] + 1 && split_part!(P, p, k)
        end
    end
    return P
end

## step 3: compute Max ##

"""
    dahlhaus_max(F, order, sizes, elems, pos, n) -> (maxset, maxsize)

For each `X ∈ F`, `Max(X)` is the largest set of `F` — in `order`, which lists
`F` by decreasing size — that overlaps `X` and is at least as large as `X`.
It is undefined when no set overlapping `X` is that large; `maxset[i]` is then
`0` and `maxsize[i]` is `0`.

`elems`/`pos` must be the lexicographic ground set order from [`lex_order`](@ref).
Writing `left(X)` and `right(X)` for the first and last positions of `X` in that
order, `Max(X)` is the set highest in `order` that contains `right(X)` but not
`left(X)` (Charbit et al. lemma 6). So we replay the refinement of `lex_order`,
which by construction never moves an element, and watch for the first split that
separates `left(X)` from `right(X)`: the set that caused it is `Max(X)`.

Sets are dropped once every set at least as large as them has been processed, so
that any split found is by a set of size ≥ |X|, as `Max` requires.
"""
function dahlhaus_max(F, order, sizes::Vector{Int}, elems::Vector{Int},
                      pos::Vector{Int}, n::Int)
    m = length(F)
    maxset, maxsize = zeros(Int, m), zeros(Int, m)

    # left and right extent of each set in the lexicographic order
    left, right = fill(n + 1, m), zeros(Int, m)
    for i = 1:m, v in F[i]
        left[i] = min(left[i], pos[v])
        right[i] = max(right[i], pos[v])
    end

    # AM[j] lists the sets whose right extent is position j, in increasing order
    # of left extent, as a doubly linked list so that any set can also be
    # dropped in O(1) once it is too large to be anyone's Max
    live = falses(m)
    head, nxt, prv = zeros(Int, n), zeros(Int, m), zeros(Int, m)
    pending = [i for i = 1:m if sizes[i] > 1] # smaller sets overlap nothing
    for i in sortperm_desc(pending, left, n + 1) # decreasing left, pushed to
        j = right[i]                             # the front: so AM[j] increases
        nxt[i], prv[i], live[i] = head[j], 0, true
        head[j] != 0 && (prv[head[j]] = i)
        head[j] = i
    end

    function drop!(i::Int)
        live[i] || return
        live[i] = false
        j = right[i]
        prv[i] == 0 ? (head[j] = nxt[i]) : (nxt[prv[i]] = nxt[i])
        nxt[i] != 0 && (prv[nxt[i]] = prv[i])
    end

    # `order` lists the sets by decreasing size, so each size class is one of
    # its runs; `first_of_size` tracks where the run being processed began
    P, met, first_of_size = Partition(elems), Int[], 1
    for (t, i) in enumerate(order)
        empty!(met)
        for v in F[i]
            p = P.part[v]
            P.count[p] == 0 && push!(met, p)
            P.count[p] += 1
        end
        for p in met
            k, P.count[p] = P.count[p], 0
            k < P.hi[p] - P.lo[p] + 1 || continue
            # the elements of F[i] already occupy the tail of the part: parts
            # are ranges of the lexicographic order, within which the members
            # of the most significant set not yet processed come last
            q = split_part!(P, p, k)
            boundary = P.hi[p]
            for j = P.lo[q]:P.hi[q]
                while head[j] != 0 && left[head[j]] ≤ boundary
                    x = head[j]
                    maxset[x], maxsize[x] = i, sizes[i]
                    drop!(x)
                end
            end
        end
        # once the last set of this size is done, nothing left can be its Max
        if t == m || sizes[order[t+1]] != sizes[i]
            for u = first_of_size:t
                drop!(order[u])
            end
            first_of_size = t + 1
        end
    end
    return maxset, maxsize
end

## the whole algorithm ##

ground_set_size(F) = mapreduce(X -> isempty(X) ? 0 : Int(maximum(X)), max, F, init=0)

"""
    overlap_classes(F, n = maximum element of F) -> (class, nclasses)

Assign each set of the family `F` the index of its overlap class, in
Θ(n + Σ|X|) time. The sets of `F` must be duplicate-free subsets of `1:n`.

Sets that overlap nothing form classes of their own, and equal sets are placed
in the same class only if something else connects them: a set never overlaps
itself or a copy of itself.
"""
function overlap_classes(
    F :: AbstractVector{<:AbstractVector{<:Integer}},
    n :: Integer = ground_set_size(F),
)
    m = length(F)
    m == 0 && return Int[], 0
    n = Int(n)
    sizes = [length(X) for X in F]

    # the algorithm silently misbehaves on out of range or repeated elements
    # rather than failing, so rule both out up front -- still one linear pass
    seen = falses(n)
    for X in F
        for v in X
            1 ≤ v ≤ n ||
                throw(ArgumentError("element $v is outside the ground set 1:$n"))
            seen[v] && throw(ArgumentError("element $v is repeated within a set"))
            seen[v] = true
        end
        for v in X
            seen[v] = false
        end
    end

    # LF order: the sets by decreasing size, ties broken arbitrarily
    order = sortperm_desc(1:m, sizes, n)
    P = lex_order(F, order, n)
    _, maxsize = dahlhaus_max(F, order, sizes, P.elems, P.pos, n)

    # SL[v] lists the sets containing v by increasing size, as a CSR array
    start = zeros(Int, n + 2)
    for i = 1:m, v in F[i]
        start[v+1] += 1
    end
    start[1] = 1
    for v = 1:n+1
        start[v+1] += start[v]
    end
    cursor, SL = copy(start), Vector{Int}(undef, start[n+2] - 1)
    for i in Iterators.reverse(order), v in F[i]
        SL[cursor[v]] = i
        cursor[v] += 1
    end

    # Dahlhaus' graph: walking SL[v] from the smallest set up, join consecutive
    # sets whenever some set already passed has a Max at least as large as the
    # later of the two. Both then overlap that set or its Max, which overlap
    # each other, so all four share an overlap class (Charbit et al. lemma 2);
    # and conversely every overlapping pair is joined by such a walk.
    parent, weight = collect(1:m), ones(Int, m)
    for v = 1:n
        best = prev = 0
        for t = start[v]:start[v+1]-1
            x = SL[t]
            prev != 0 && sizes[x] ≤ best && uf_union!(parent, weight, prev, x)
            best = max(best, maxsize[x])
            prev = x
        end
    end

    class, nclasses = zeros(Int, m), 0
    for i = 1:m
        r = uf_find(parent, i)
        class[r] == 0 && (class[r] = nclasses += 1)
        class[i] = class[r]
    end
    return class, nclasses
end

"""
    overlap_components(F, n = maximum element of F) -> Vector{Vector{Int}}

The union of each overlap class of the family `F`, as a sorted vector, in
Θ(n + Σ|X|) time. The sets of `F` must be duplicate-free subsets of `1:n`.

Any two of the unions returned are disjoint, nested, or equal: a class whose
union coincides with a larger set that overlaps nothing is reported separately.
"""
function overlap_components(
    F :: AbstractVector{<:AbstractVector{<:Integer}},
    n :: Integer, # the ground set is 1:n; see the one-argument methods below
)
    n = Int(n)
    class, nclasses = overlap_classes(F, n)
    nclasses == 0 && return Vector{Int}[]

    # bucket the (class, element) incidences by element and then by class, so
    # that each class ends up with its elements in increasing order
    N = sum(length, F)
    cs, es = Vector{Int}(undef, N), Vector{Int}(undef, N)
    t = 0
    for i in eachindex(F), v in F[i]
        t += 1
        cs[t], es[t] = class[i], v
    end
    for (key, kmax) in ((es, n), (cs, nclasses))
        count = zeros(Int, kmax + 1)
        for k in key
            count[k+1] += 1
        end
        next, acc = Vector{Int}(undef, kmax + 1), 1
        for k = 0:kmax
            next[k+1] = acc
            acc += count[k+1]
        end
        cs′, es′ = similar(cs), similar(es)
        for t = 1:N
            cs′[next[key[t]+1]], es′[next[key[t]+1]] = cs[t], es[t]
            next[key[t]+1] += 1
        end
        copyto!(cs, cs′)
        copyto!(es, es′)
    end

    components = [Int[] for _ = 1:nclasses]
    for t = 1:N
        U = components[cs[t]]
        (isempty(U) || U[end] != es[t]) && push!(U, es[t])
    end
    return components
end

"""
    overlap_components(F) -> Vector{Vector}

The union of each overlap class of the family `F`, as a sorted vector. The sets
may be over any ground set; integer elements are used as ground set indices
directly, anything else is relabeled first, at a cost of O(Σ|X| log Σ|X|).
"""
function overlap_components(F::AbstractVector{<:AbstractVector{<:Integer}})
    lo = mapreduce(X -> isempty(X) ? typemax(Int) : Int(minimum(X)), min, F, init=typemax(Int))
    lo ≥ 1 && return overlap_components(F, ground_set_size(F))
    shift(X, d) = [x + d for x in X]
    U = overlap_components([shift(X, 1 - lo) for X in F], ground_set_size(F) + 1 - lo)
    return [shift(X, lo - 1) for X in U]
end

function overlap_components(F::AbstractVector{<:AbstractVector{T}}) where {T}
    elements = sort!(unique!(collect(T, Iterators.flatten(F))))
    index = Dict(x => i for (i, x) in enumerate(elements))
    U = overlap_components([[index[x] for x in X] for X in F], length(elements))
    return [elements[X] for X in U]
end

end # module

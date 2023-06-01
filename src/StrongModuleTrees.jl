module StrongModuleTrees

export StrongModuleTree, strong_modules,
    leaves, leaf_count, first_leaf, last_leaf, rand_leaf,
    nodes, node_count, parent_node, is_cograph

using ..GraphModularDecomposition

## core type definition ##

struct StrongModuleTree{T} <: AbstractVector{T}
    kind::Symbol
    edge::Tuple
    nodes::Vector{Union{T,StrongModuleTree{T}}}
end

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

    # remove non-module "dummy" nodes
    let s = Int[]
        for j = 1:n
            for _ = 1:op[j]; push!(s, j); end
            for _ = 1:cl[j]
                i = pop!(s)
                if i < j
                    l = minimum(lc[k] for k = i:j-1)
                    u = maximum(uc[k] for k = i:j-1)
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

    function classify_nodes(t::Vector)
        n = length(t)
        counts = zeros(Int, n)
        x, y = first_leaf(t[1]), first_leaf(t[2])
        edge = (G[y,x], G[x,y])
        local a, b
        for i = 1:n, j = 1:n
            i == j && continue
            x, y = first_leaf(t[i]), first_leaf(t[j])
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
        StrongModuleTree{T}(kind, edge, map(x->x isa Vector ? classify_nodes(x) : x, t))
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
    t = classify_nodes(s[end])
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

function leaves(t::StrongModuleTree{T}) where T
    L = T[]
    for x in t.nodes
        append!(L, leaves(x))
    end
    return L
end
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

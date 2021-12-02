using Combinatorics
using GraphModularDecomposition
using GraphModularDecomposition: overlap
using IterTools
using LinearAlgebra
using LinearAlgebra: checksquare
using Random
using SparseArrays
using Test

const \ = setdiff

function is_module(G::AbstractMatrix, S::Vector{Int})
    n = checksquare(G)
    !isempty(S) &&
    all(G[i,k] == G[j,k] && G[k,i] == G[k,j]
        for i in S for j in S for k in 1:n if k ∉ S)
end

function all_modules(G::AbstractMatrix)
    n = checksquare(G)
    filter!(collect(subsets(1:n))) do S
        1 < length(S) && is_module(G, S)
    end
end

function all_strong_modules(G::AbstractMatrix)
    modules = all_modules(G)
    filter(modules) do A
        all(modules) do B
            !overlap(A, B)
        end
    end
end

function is_modular_permutation(G::AbstractMatrix, p::Vector{Int};
    modules = all_strong_modules(G))
    isempty(modules) && return true
    diffs = map(modules) do M
        diff(findall(in(M), p))
    end
    maximum(maximum, diffs) == 1
end

findin_partition(P, S) = sort!(map(x->findfirst(X->x in X, P), S))

function is_modular_partition(G::AbstractMatrix, P::Vector{Vector{Int}};
    modules = all_strong_modules(G))
    sort!(vcat(P...)) == collect(1:size(G,2)) || error("not a partition")
    diffs = map(modules) do M
        diff(findin_partition(P, M))
    end
    maximum(maximum, diffs) == 1
end

is_tournament(G::AbstractMatrix) = all(G + G' + I .== 1)

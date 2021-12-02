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

is_module(G::AbstractMatrix, S::Vector{Int}) = !isempty(S) &&
    all(G[i,k] == G[j,k] && G[k,i] == G[k,j]
        for i in S for j in S for k in axes(G,2)\S)

function all_modules(G::AbstractMatrix)
    n = checksquare(G)
    filter!(S->1 < length(S) && is_module(G, S), collect(subsets(1:n)))
end

function all_strong_modules(G::AbstractMatrix)
    modules = all_modules(G)
    filter(A -> all(B -> !overlap(A, B), modules), modules)
end

function is_modular_permutation(G::AbstractMatrix, p::Vector{Int};
    modules = all_strong_modules(G))
    isempty(modules) && return true
    diffs = map(M->diff(findall(in(M), p)), modules)
    maximum(maximum, diffs) == 1
end

findin_partition(P, S) = sort!(map(x->findfirst(X->x in X, P), S))

function is_modular_partition(G::AbstractMatrix, P::Vector{Vector{Int}};
    modules = all_strong_modules(G))
    sort!(vcat(P...)) == collect(1:size(G,2)) || error("not a partition")
    diffs = map(M->diff(findin_partition(P, M)), modules)
    maximum(maximum, diffs) == 1
end

is_tournament(G::AbstractMatrix) = G + G' + I == fill(true, size(G))

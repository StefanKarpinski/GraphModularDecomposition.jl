using GraphModularDecomposition
using LinearAlgebra
using Random
using SparseArrays
using Test

include("testutils.jl")
using .TestUtils

@testset "factorizing permutations for symmetric graphs" begin
    for _ = 1:100
        n = rand(3:10)
        G = rand(0:1, n, n)
        G .= G .⊻ G'
        @test G == G'
        p = symgraph_factorizing_permutation(G)
        @test is_modular_permutation(G, p)
        T = sort!(StrongModuleTree(G, p))
        @test sort!(strong_modules(T)) == sort!(all_strong_modules(G))
        for _ = 1:10
            p′ = symgraph_factorizing_permutation(G, shuffle(1:n))
            @test is_modular_permutation(G, p′)
            T′ = sort!(StrongModuleTree(G, p′))
            @test T == T′
        end
    end
end

@testset "factorizing permutations for tournament graphs" begin
    for _ = 1:1000
        n = rand(3:10)
        T = Int[i != j && rand(Bool) for i=1:n, j=1:n]
        T .= T .⊻ T' .⊻ tril(ones(Int,n,n),-1)
        @test is_tournament(T)
        p = tournament_factorizing_permutation(T)
        modules = all_modules(T)
        @test is_modular_permutation(T, p, modules=modules)
    end
end

@testset "factorizing permutations for directed graphs" begin
    for _ = 1:1000
        n = rand(3:10)
        G = Int[i != j && rand(Bool) for i=1:n, j=1:n]
        p = digraph_factorizing_permutation(G)
        @test is_modular_permutation(G, p)
        p = graph_factorizing_permutation(G)
        @test is_modular_permutation(G, p)
    end
end

function test_permutations(G::AbstractMatrix, T::StrongModuleTree, N::Integer=100)
    n = LinearAlgebra.checksquare(G)
    for _ = 1:N
        p = shuffle(1:n)
        G′ = G[p,p]
        @test sort!(p[StrongModuleTree(G′)]) == T
    end
end

@testset "symmetric graph example [Wikipedia]" begin
    G = sparse(
        [1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4, 4, 5,
         5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 6, 6, 7, 7, 7, 7, 7, 7, 7, 7,
         8, 8, 8, 8, 8, 9, 9, 9, 9, 9, 10, 10, 10, 10, 11, 11, 11, 11],
        [2, 3, 4, 1, 4, 5, 6, 7, 1, 4, 5, 6, 7, 1, 2, 3, 5, 6, 7, 2,
         3, 4, 6, 7, 2, 3, 4, 5, 8, 9, 10, 11, 2, 3, 4, 5, 8, 9, 10,
         11, 6, 7, 9, 10, 11, 6, 7, 8, 10, 11, 6, 7, 8, 9, 6, 7, 8, 9],
        1
    )
    T = sort!(StrongModuleTree(G))
    @test repr(T) == "{1 ((2 3) 4) 5 (6 7) (8 9 (10 11))}"
    test_permutations(G, T)
end

@testset "directed graph example [Capelle, Habib, Montgolfier 2002]" begin
    G = sparse(
        [1, 1, 1, 2, 3, 3, 3, 4, 5, 5, 6, 6, 6, 6, 7, 7, 7, 8, 8, 10, 10,
         10, 10, 11, 11, 11, 12, 12, 12, 12, 13, 13, 13, 13, 14, 14, 14],
        [3, 4, 5, 1, 2, 4, 5, 2, 4, 2, 7, 8, 9, 10, 8, 9, 10, 9, 10, 11,
         12, 13, 14, 9, 10, 14, 9, 10, 11, 13, 9, 10, 11, 12, 9, 10, 11],
        1
    )
    T = sort!(StrongModuleTree(G))
    @test repr(T) == "({1 2 [3 4 5]} {[6 7 8] 9 10 {11 (12 13) 14}})"
    test_permutations(G, T)
end

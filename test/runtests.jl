using GraphModularDecomposition
using LinearAlgebra
using Random
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

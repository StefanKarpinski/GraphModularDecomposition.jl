# Empirical check that OverlapComponents runs in time linear in the size of its
# input, n + Σ|X|. Too slow for the standard test suite, so it lives here:
#
#     julia bench/overlap_scaling.jl               # wall clock scaling
#     julia bench/overlap_scaling.jl --instructions # retired instruction counts
#
# The module has no dependencies, so this needs no project environment.
#
# The wall clock tables report time per unit of input, which is flat iff the
# implementation is linear, and the log-log slope of time against input size,
# which is 1.0 iff it is linear. At these sizes the working set outgrows every
# level of cache, so measured time keeps a memory hierarchy ramp in it that no
# amount of algorithmic linearity removes; expect slopes a little above 1.
#
# The --instructions mode settles that ambiguity where perf(1) is available.
# Retired instructions are unaffected by cache misses -- a miss costs cycles,
# not instructions -- so instructions per unit of input is flat exactly when the
# work is linear. Each measurement runs the algorithm once and again five times
# in a fresh process and takes the difference, which cancels out interpreter
# startup, compilation and building the input family.

include(joinpath(@__DIR__, "..", "src", "OverlapComponents.jl"))
using .OverlapComponents
using Printf, Random, Statistics

## family shapes, each parameterized by a single size knob ##

# a star: every set meets set 1 and nothing else, so the overlap graph has
# Θ(m²) edges while the input has only 2m. The case Dahlhaus' detour exists for.
star(s) = (s, [[1, i] for i = 2:s])

# a chain of sliding windows: a single overlap class spanning everything
windows(s, w=8) = (s, [collect(i:i+w-1) for i = 1:s-w])

# random fixed size subsets: many classes, many near misses
function random_sets(s, k=8)
    rng = MersenneTwister(s)
    (s, map(1:s) do _
        X = Int[]
        while length(X) < min(k, s)
            v = rand(rng, 1:s)
            v in X || push!(X, v)
        end
        sort!(X)
    end)
end

# a nested chain: laminar, nothing overlaps anything, but Σ|X| is quadratic in
# the ground set -- so this shape scales the input size without scaling n
nested(s) = (s, [collect(1:i) for i = 1:s])

# two balanced hierarchies over one ground set: the shape that modular
# decomposition of a digraph actually feeds to overlap components
function hierarchy(V, out=Vector{Int}[])
    push!(out, V)
    length(V) > 1 || return out
    k = length(V) ÷ 2
    hierarchy(V[1:k], out)
    hierarchy(V[k+1:end], out)
    return out
end

function two_trees(s)
    rng = MersenneTwister(s)
    (s, [hierarchy(randperm(rng, s)); hierarchy(randperm(rng, s))])
end

# ladders chosen so that the largest input in each shape is a few million
const SHAPES = [
    ("star (Θ(m²) overlap graph edges)", star,        [1, 2, 4, 8, 16, 32] .* 10^5),
    ("sliding windows",                  windows,     [1, 2, 4, 8, 16, 32] .* 10^5),
    ("random 8-element subsets",         random_sets, [25, 50, 100, 200, 400, 800] .* 10^3),
    ("nested chain (Σ|X| = Θ(n²))",      nested,      [750, 1500, 3000, 6000, 12000]),
    ("two balanced hierarchies",         two_trees,   [25, 50, 100, 200, 400] .* 10^3),
]

shape_named(name) = only(s for (_, s, _) in SHAPES if string(s) == name)

## the quadratic baseline this replaced ##

function pairwise_components(F::Vector{Vector{Int}})
    m = length(F)
    S = map(sort, F)
    class, nclasses = zeros(Int, m), 0
    for i = 1:m
        class[i] == 0 || continue
        nclasses += 1
        stack = [i]
        while !isempty(stack)
            x = pop!(stack)
            class[x] == 0 || continue
            class[x] = nclasses
            for j = 1:m
                class[j] == 0 && overlap(S[x], S[j]) && push!(stack, j)
            end
        end
    end
    U = [Int[] for _ = 1:nclasses]
    for (i, X) in enumerate(F)
        append!(U[class[i]], X)
    end
    return sort!(map(sort!∘unique!, U))
end

## wall clock scaling ##

function best_time(f, args...; trials=5)
    f(args...) # warm up and compile
    minimum(1:trials) do _
        GC.gc()
        @elapsed f(args...)
    end
end

# least squares slope of log(y) against log(x): 1.0 means linear
function fitted_exponent(x, y)
    lx, ly = log.(x), log.(y)
    sum((lx .- mean(lx)) .* (ly .- mean(ly))) / sum((lx .- mean(lx)).^2)
end

function report(name, shape, scales)
    @printf("\n%s\n", name)
    @printf("  %10s %12s %12s %10s %9s %7s\n",
            "n", "sets", "input size", "time", "ns/input", "slope")
    sizes, times, previous = Float64[], Float64[], nothing
    for s in scales
        n, F = shape(s)
        N = n + sum(length, F)
        t = best_time(overlap_components, F, n)
        push!(sizes, N); push!(times, t)
        slope = previous === nothing ? NaN :
            log(t / previous[2]) / log(N / previous[1])
        @printf("  %10d %12d %12d %10.4fs %9.1f %7s\n", n, length(F), N, t,
                1e9t / N, isnan(slope) ? "-" : @sprintf("%.2f", slope))
        flush(stdout)
        previous = (N, t)
    end
    @printf("  fitted exponent over the whole ladder: %.2f\n",
            fitted_exponent(sizes, times))
    flush(stdout)
end

function crossover(shape, scales)
    @printf("\nstar: linear vs quadratic\n")
    @printf("  %10s %12s %12s %12s %10s\n",
            "n", "input size", "Dahlhaus", "pairwise", "speedup")
    for s in scales
        n, F = shape(s)
        N = n + sum(length, F)
        fast = best_time(overlap_components, F, n)
        slow = best_time(pairwise_components, F, trials=1)
        @printf("  %10d %12d %11.4fs %11.4fs %9.0fx\n", n, N, fast, slow, slow/fast)
        flush(stdout)
    end
end

function verify()
    print("verifying against the pairwise baseline ... ")
    for (name, shape, _) in SHAPES, s in (50, 100, 200)
        n, F = shape(s)
        sort(overlap_components(F, n)) == pairwise_components(F) ||
            error("mismatch on $name at scale $s")
    end
    println("ok")
    flush(stdout)
end

## retired instruction counts, via perf(1) ##

# run the algorithm `reps` times on one family and report nothing: this process
# is measured from the outside
function run_only(name, size, reps)
    n, F = shape_named(name)(size)
    total = 0
    for _ = 1:reps
        total += length(overlap_components(F, n))
    end
    total > 0 || error("no components")
end

function perf_instructions(name, size, reps)
    self = joinpath(@__DIR__, basename(@__FILE__))
    out = read(pipeline(`perf stat -x, -e instructions --
                         $(Base.julia_cmd()) $self --run $name $size $reps`,
                        stderr = `cat`), String)
    for line in split(out, '\n')
        fields = split(line, ',')
        length(fields) ≥ 3 && fields[3] == "instructions" &&
            return parse(Float64, fields[1])
    end
    error("could not read an instruction count from perf:\n$out")
end

function instruction_report()
    Sys.which("perf") === nothing && error("--instructions needs perf(1)")
    println("retired instructions per unit of input, measured by perf(1)")
    println("(a flat column is linear work, independent of the memory hierarchy)")
    for (name, shape, scales) in SHAPES
        @printf("\n%s\n", name)
        @printf("  %10s %12s %18s %12s\n", "n", "input size", "instructions",
                "instr/input")
        sizes, counts = Float64[], Float64[]
        for s in scales
            n, F = shape(s)
            N = n + sum(length, F)
            one = perf_instructions(string(shape), s, 1)
            six = perf_instructions(string(shape), s, 6)
            instructions = (six - one) / 5
            push!(sizes, N); push!(counts, instructions)
            @printf("  %10d %12d %18.0f %12.1f\n", n, N, instructions,
                    instructions / N)
            flush(stdout)
        end
        @printf("  fitted exponent over the whole ladder: %.3f\n",
                fitted_exponent(sizes, counts))
        flush(stdout)
    end
end

## entry points ##

if length(ARGS) ≥ 4 && ARGS[1] == "--run"
    run_only(ARGS[2], parse(Int, ARGS[3]), parse(Int, ARGS[4]))
elseif length(ARGS) ≥ 1 && ARGS[1] == "--instructions"
    instruction_report()
else
    verify()
    for (name, shape, ladder) in SHAPES
        report(name, shape, ladder)
    end
    crossover(star, [500, 1_000, 2_000, 4_000])
end

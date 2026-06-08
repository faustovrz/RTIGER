# Reproducible benchmark / equivalence harness for RTIGER's Julia core.
# Usage: julia harness.jl <path-to-julia-src-dir> <nsamples> <T> <rigidity> <seed> <tag>
# Sources the two core .jl files from the given dir, generates synthetic
# observations, runs `fit`, prints timings, and writes the fitted params +
# viterbi path to bench/out_<tag>.txt so two implementations can be diffed.

using Random, LinearAlgebra, Printf

srcdir   = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "..", "inst", "julia")
nsamples = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 1
T        = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 1000
rigidity = length(ARGS) >= 4 ? parse(Int, ARGS[4]) : 50
seed     = length(ARGS) >= 5 ? parse(Int, ARGS[5]) : 42
tag      = length(ARGS) >= 6 ? ARGS[6] : "baseline"

include(joinpath(srcdir, "AuxilaryFunctions.jl"))
include(joinpath(srcdir, "rHMM_methods.jl"))

# --- synthetic data generation (3 states: mat/het/pat) ---------------------
# Build a realistic r-rigid-ish state path, then emit [k, n] counts.
function gen_obs(nsamples::Int, T::Int, rigidity::Int, seed::Int; coverage=1.0, nchr=1)
    rng = MersenneTwister(seed)
    pstate = [0.9, 0.5, 0.1]            # ref-allele fraction per state
    O = Dict{Any,Any}()
    for s in 1:nsamples
        chrs = Dict{Any,Any}()
        for c in 1:nchr
            # random segment path with min length ~ rigidity
            states = Int[]
            cur = rand(rng, 1:3)
            while length(states) < T
                seglen = rigidity + rand(rng, 0:3*rigidity)
                append!(states, fill(cur, seglen))
                cur = rand(rng, setdiff(1:3, cur))
            end
            states = states[1:T]
            obs = Matrix{Int64}(undef, T, 2)
            for t in 1:T
                n = rand(rng, 0:1) == 0 ? 0 : 1 + rand(rng, 0:max(0, round(Int, 2coverage)))
                p = pstate[states[t]]
                k = 0
                for _ in 1:n
                    k += rand(rng) < p ? 1 : 0
                end
                obs[t, 1] = k
                obs[t, 2] = n
            end
            chrs["Chr$c"] = obs
        end
        O["S$s"] = chrs
    end
    return O
end

obs = gen_obs(nsamples, T, rigidity, seed)
param0 = erstePara(3, rigidity; aBeta=[0.1,1.0,1.9], bBeta=[1.9,1.0,0.1], equalstart=true)

# warmup (compile) on a tiny copy
let
    warm = gen_obs(1, 200, min(rigidity, 20), 1)
    p = erstePara(3, min(rigidity, 20); aBeta=[0.1,1.0,1.9], bBeta=[1.9,1.0,0.1], equalstart=true)
    fit(warm, nothing, p, 3, 1e-3, false, true, false, 20, nothing, true)
end

@printf("=== RUN tag=%s nsamples=%d T=%d rigidity=%d ===\n", tag, nsamples, T, rigidity)
GC.gc()
t0 = time()
res = fit(obs, nothing, deepcopy(param0), 30, 1e-4, false, true, false, 20, nothing, true)
elapsed = time() - t0
@printf("fit total: %.4f s  (iterations=%d)\n", elapsed, res[:numberofiterations])
@printf("alloc:    %.2f MiB\n", (@allocated fit(obs, nothing, deepcopy(param0), 30, 1e-4, false, true, false, 20, nothing, true)) / 2^20)

# --- dump results for equivalence check ------------------------------------
open(joinpath(@__DIR__, "out_$(tag).txt"), "w") do f
    p = res[:parameterSet]
    println(f, "alpha=", round.(vec(p[:paraBetaAlpha]); digits=6))
    println(f, "beta=",  round.(vec(p[:paraBetaBeta]);  digits=6))
    println(f, "pi=",    round.(vec(p[:pi]);            digits=6))
    println(f, "transition=", round.(vec(p[:transition]); digits=6))
    println(f, "iters=", res[:numberofiterations])
    for s in sort(collect(keys(res[:viterbiPath])))
        for c in sort(collect(keys(res[:viterbiPath][s])))
            println(f, "vit[$s][$c]=", res[:viterbiPath][s][c])
        end
    end
end
println("wrote out_$(tag).txt")

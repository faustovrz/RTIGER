# Time each EM component separately at realistic per-sample size.
using Random, LinearAlgebra, Printf
srcdir = length(ARGS) >= 1 ? ARGS[1] : joinpath(@__DIR__, "..", "..", "inst", "julia")
T      = length(ARGS) >= 2 ? parse(Int, ARGS[2]) : 50000
r      = length(ARGS) >= 3 ? parse(Int, ARGS[3]) : 250
include(joinpath(srcdir, "AuxilaryFunctions.jl"))
include(joinpath(srcdir, "rHMM_methods.jl"))

function genstates(T, r, rng)
    states = Int[]; cur = 1
    while length(states) < T
        append!(states, fill(cur, r + rand(rng, 0:3r)))
        cur = rand(rng, setdiff(1:3, cur))
    end
    return states[1:T]
end
rng = MersenneTwister(7)
states = genstates(T, r, rng)
pstate = [0.9,0.5,0.1]
obs = Matrix{Int64}(undef, T, 2)
for t in 1:T
    n = rand(rng,0:1)==0 ? 0 : 1+rand(rng,0:2)
    k = sum(rand(rng) < pstate[states[t]] for _ in 1:n; init=0)
    obs[t,1]=k; obs[t,2]=n
end

a=[0.1,1.0,1.9]; b=[1.9,1.0,0.1]
A = 0.1*ones(3,3)+10*Diagonal(ones(3)); A=A./sum(A,dims=2)
logA=log.(A); logpi=log.(fill(1/3,3)); s=3

bench(name, f) = (f(); GC.gc(); t=time(); for _ in 1:3; f(); end; el=(time()-t)/3;
    al=(@allocated f())/2^20; @printf("%-14s %8.4f s   %9.2f MiB\n", name, el, al))

@printf("=== phase profile  T=%d r=%d ===\n", T, r)
local lp, lP, al, be, ze
bench("getlogpsi",  () -> (global lp = getlogpsi(obs,a,b)))
bench("productpsi", () -> (global lP = productpsi(lp,r)))
bench("forward",    () -> (global al = forward(T,r,s,logpi,lP,logA,lp)))
bench("backward",   () -> (global be = backward(T,r,s,lP,logA,lp)))
bench("zeta",       () -> (global ze = zeta(al,be,logA,lP,lp,r)))
bench("gamma",      () -> gamma(ze,al,be,r))
bench("viterbi",    () -> viterbi(logpi,lP,lp,logA,r))

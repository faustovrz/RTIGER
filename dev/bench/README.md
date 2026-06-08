# Julia core benchmark / equivalence harness

Tooling to (a) measure the cost of RTIGER's Julia EM core and (b) prove that
performance refactors produce identical output.

## Files
- `harness.jl` — generates synthetic 3-state (mat/het/pat) observations, runs
  the full `fit`, prints wall time + allocations, and dumps fitted parameters
  and the Viterbi path to `out_<tag>.txt` for diffing.
- `profile_phases.jl` — times each EM component (getlogpsi, productpsi,
  forward, backward, zeta, gamma, viterbi) separately at a realistic
  per-sample size.
- `out_baseline.txt` — reference output from upstream `main` (the equivalence
  oracle). Any refactor must reproduce this.

## Usage
```bash
# from this directory; arg 1 is the julia source dir (defaults to ../../inst/julia)
julia harness.jl ../../inst/julia 1 1000 50 42 baseline
julia profile_phases.jl ../../inst/julia 50000 250

# after a refactor, regenerate and diff:
julia harness.jl ../../inst/julia 1 1000 50 42 optimized
diff out_baseline.txt out_optimized.txt   # must be empty (Viterbi path identical)
```

## Baseline (upstream main, this machine)
- `fit` 1 sample × 1000 markers, r=50: ~2.0 s, ~427 MiB allocated, 7 iters.
- Per sample per EM iteration at 50k markers, r=250: ~94 ms, ~282 MiB.
  Dominated by `forward` (111 MiB) and `backward` (104 MiB) — per-position
  slicing and array-returning log-sum-exp helpers in the hot loop.

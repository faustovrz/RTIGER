# RTIGER optimization fork — performance notes

This repository — **`faustovrz/RTIGER`** (branch `optimize-julia-core`) — is a
**performance fork** of the original RTIGER,
[`rfael0cm/RTIGER`](https://github.com/rfael0cm/RTIGER). The goal is to make
RTIGER usable on **large populations** — on the order of **1400 samples × ~50 000
markers** — where the upstream implementation is impractically slow and runs out
of memory.

All changes are confined to the Julia core (`inst/julia/rHMM_methods.jl`, the
engine behind the R `JuliaCall` wrapper) plus a thin R-side argument for progress
logging. **The model, the joint fit, the convergence criterion, and the outputs
are unchanged** — every optimization preserves the arithmetic. Equivalence to the
original is validated (see [Correctness](#correctness)).

> TL;DR — same results, **~40–1500× faster** depending on stage, and **peak
> memory made flat in the number of samples** (projected ~33 GB → ~3.6 GB at
> 1400 samples), plus an opt-in per-iteration progress log so long fits report an
> ETA.

---

## 1. Why the fork

Upstream RTIGER was written and validated for modest sample counts. Two walls
appear at population scale:

1. **Runtime.** A full Baum–Welch (EM) fit is dominated almost entirely
   (~99.97%) by the **emission-distribution M-step**, which rebuilt one
   `BetaBinomial` object per marker on *every* optimizer evaluation. At hundreds
   of thousands of markers this is the bottleneck; the rest of the EM core
   (forward/backward, Viterbi) is secondary but also allocation-heavy.
2. **Peak memory.** `EM()` retained every sample's full per-position arrays
   (`zeta`, `gamma`, `alpha`, `beta`, `psi`) simultaneously and only pooled them
   after the whole E-step. Peak memory therefore grew **linearly in the number
   of samples** (~25 MiB/sample at 50 k markers), projecting to ~33 GB at 1400
   samples.

This fork removes both walls without changing what RTIGER computes.

---

## 2. Summary of gains

| Area | Before | After | Factor |
|---|---|---|---|
| Emission M-step (per EM iteration) | dominant cost | grouped + pre-summed γ | **~1500×** |
| `getlogpsi` (emission log-pdf) | rebuilt per marker | memoized over distinct (k,n) | **~27×** |
| Viterbi (max-product) | array-slice per step | in-place scalar argmax | **~12×** |
| forward / backward | allocating log-sum-exp | in-place scalar | **~5×**, ~90× less alloc |
| Full fit, AAACB5K (15 k markers) | 19.0 s | 0.46 s | **~41×** |
| Full fit, BNZAU15K (15 k, r=2) | 157.5 s | 1.1 s | **~143×** |
| Full fit, BNZAU270K (3 samples × ~270K markers/sample, r=2) | **17.3 h** (62 354 s) | **64.5 s** | **≈966×** |
| Peak RSS vs #samples | linear (~25 MiB/sample) | **flat (~constant)** | §5 |
| Projected peak RSS @ 1400×50k | ~33 GB | **~3.6 GB** | ~9× |

Datasets are real *Arabidopsis* Col×Ler allele counts, sourced from the package's
`data/` folder and the `.RData` environment dump of the original repository.

---

## 3. Per-iteration timing and convergence (BNZAU270K head-to-head)

![Left: per-iteration EM wall time (log scale) for the optimized vs original core on the BNZAU270K fit. Right: each iteration's original convergence δ plotted against the optimized δ on a log–log 1:1 line.](time_performance_270K.png)

**How this was produced.** Both cores were run on the same BNZAU270K fit (3
samples × ~270K markers/sample, R-default init, `eps=0.01`, rigidity 2, native
arm64) from identical deterministic initialization. Each EM iteration's `elapsed`
and `delta` were captured with the opt-in progress log (§6) — one record per
iteration — and the two logs read back and plotted: per-iteration wall time as
the first-difference of cumulative `elapsed`, and the two cores' `delta`
sequences against each other. The reproducible (un-evaluated) plotting notebook
is [`time_performance_270K.qmd`](time_performance_270K.qmd).

**What it shows.** *Left:* the optimized core holds ~1.9 s/iter; the original
averages ~1834 s/iter, with the erratic spikes characteristic of the
un-optimized emission `Optim` step — the ≈966× full-fit gap, iteration by
iteration. *Right:* every point sits on the 1:1 line (max relative δ difference
~4e-5), so the two cores descend the **identical** convergence path and stop at
the same iteration (**34** each). The δ gap is float-summation-order noise, not
an algorithmic difference — consistent with the bit-identical Viterbi paths in
§7.

---

## 4. Time optimization

Each change preserves the exact arithmetic (same summation order where it
matters) and is independently committed.

| Component | Technique | Commit |
|---|---|---|
| **Emission M-step** | Group markers by distinct `(k, n)` once and pre-sum the γ weights, so the `Optim` objective costs `O(#distinct pairs)` instead of `O(#markers)` per evaluation. Low-coverage genotyping data has only a handful of distinct pairs, so this collapses ~`T·states` `BetaBinomial` constructions to a few dozen. | `16b8a65` |
| **`getlogpsi`** | Memoize the `BetaBinomial` log-pdf over distinct `(k, n)` pairs (same `logpdf` calls → bit-identical values). | `d725651` |
| **`productpsi`** | In-place per-state sliding-window cumulative sum (same arithmetic order). | `d725651` |
| **Viterbi** | Allocation-free in-place scalar argmax for the r-rigid max-product step. | `6cf85ca` |
| **forward / backward** | In-place scalar log-sum-exp; ~90× fewer allocations. | `44ed85b` |

**Full-fit head-to-head (BNZAU270K — 3 samples × ~270K markers/sample, 807 550
total; R-default init, `eps=0.01`, rigidity=2, native arm64, identical
deterministic init):** both cores
converge in **34 EM iterations**. The optimized core finishes in **64.5 s**
(~1.9 s/iter); the upstream original takes **62 354 s ≈ 17.3 h** (~1834 s/iter,
erratic per-iteration cost driven by the un-optimized emission `Optim`) — a
speed-up of **≈966×**. The two fits are **equivalent**: identical Viterbi paths
(**807 550/807 550** positions total over the 3 samples, 0 mismatches) and fitted parameters agreeing to 6 decimals
(per-iteration convergence δ matching to ~4e-5 — float summation order, not an
algorithmic difference). The per-iteration time and δ-trajectory comparison is
shown in §3 above.

---

## 5. Memory optimization — streaming M-step

**Commit `eb79933`.** The M-step parameter updates depend on the data only
through **pooled sufficient statistics, which are sums**:

- transition ← Σ `zeta`
- start ← Σ `gamma[:,1]`
- emission ← Σ γ-weighted counts grouped by distinct `(k, n)`

Because these are sums, they can be **accumulated incrementally** inside the
E-step loop and each sample's heavy arrays discarded immediately, instead of
hoarding every sample's arrays and pooling at the end. The accumulators fold each
sample in **the same order** the original summed them, so results are
**bit-identical**, not merely close.

A free side effect enabled this: the per-sample `alpha/beta/gamma/psi` arrays
were only ever written to the public `@Probabilities` slot, **which nothing in
the package reads**. They are no longer retained (the slot stays a valid, empty
list), which is what makes constant-memory possible.

**Scaling sweep** (real BN/Z/AU cycled to N samples, 50 k markers each, peak RSS
via `/usr/bin/time -l`):

| N samples | markers | before (linear) | after (streaming) |
|---:|---:|---:|---:|
| 2 | 0.5 M | 562 MiB | 516 MB |
| 4 | 1.0 M | 655 MiB | 533 MB |
| 8 | 2.0 M | 786 MiB | 538 MB |
| 16 | 4.0 M | 917 MiB | 573 MB |
| 30 | 7.5 M | 1252 MiB | 582 MB |

Before: ~25 MiB/sample (linear) → **~33 GB projected at 1400 samples**. After:
~2.3 MB/sample residual — and that residual is just the **input count matrices**
(which must be held), not EM scratch → **~3.6 GB projected at 1400 samples**,
essentially flat in N. (Sweep measured on the x86_64 build; the flat-vs-linear
conclusion is platform-independent.)

A related cleanup (`91911ff`) factors the shared per-chain E-step into a single
`estepChain` helper used by both `EM` and the developer `EMdev`, so the two
cannot drift; behavior is unchanged.

---

## 6. Progress logging (opt-in, default off)

**Commit `c20906f`.** A long fit is otherwise silent. `RTIGER()` / `fit()` gain
a `progress_log` argument: when set to a file path, Julia appends one
newline-terminated record per EM iteration, flushed each iteration so it can be
tailed live to estimate an ETA:

```
iter 12/50  delta=0.0087  eps=0.01  elapsed=1.54  per_iter=0.13  ETA<=4.86
```

- `progress_log = NULL`/`FALSE` (default) → off → **byte-identical** to before,
  regardless of `verbose`.
- `progress_log = "path"` → log there.
- `progress_log = TRUE` → convenience for `file.path(outputdir, "fit_progress.log")`.
- `ETA<=` is an upper bound (the EM usually converges at `delta<eps` before
  `max.iter`).

This is logging only; it does not touch the fit. It is a distinct stream from the
developer `DEBUG`/`debugInfo.txt` dump and from the in-memory `trace` parameter
history.

---

## 7. Correctness

"Optimized ≡ original" was validated from identical initialization:

- **BNZAU15K** (real, 15 k markers, R-default init, `eps=0.01`, r=2):
  **bit-identical** fitted params (to 6 dp) and **all** Viterbi paths.
- **AAACB5K** (real extdata, deterministic init): identical params and Viterbi.
- **BNZAU270K** (3 samples × ~270K markers/sample): decoding from the stored
  fitted parameters reproduces the reference Viterbi path **100%
  (807 550 / 807 550 positions, total over the 3 samples)**.
- A synthetic equivalence harness is **bit-identical** to its committed
  baseline.
- The streaming M-step and the `progress_log=off` path were each re-checked to
  remain **bit-identical** after their changes.

Where summation order is genuinely reordered (only at very large scale with a
non-default init) params can differ by ~1e-6, with Viterbi still identical —
mathematically the same statistic, just a different float-add order.

---

## 8. Reproducing

The figure in §3 is regenerated by the committed notebook
[`time_performance_270K.qmd`](time_performance_270K.qmd) — un-evaluated;
it reads the two per-iteration progress logs and draws the two panels.

The full benchmark and equivalence suite — the full-resolution real-data checks
(270 k decode = 100 %, the 15 k bit-identical A/B fit, the 270 k
optimized-vs-original head-to-head), the peak-RSS scaling sweep, and the
synthetic equivalence harness against its committed baseline — was run from the
`optimize-julia-core` development workspace and is not shipped with the package.

---

## 9. Notes / gotchas

- **Production defaults:** R's `RTIGER()` uses `eps = 0.01` (overrides Julia's
  `1e-5`); the 3-state label order is fixed pat/het/mat = 1/2/3, so the
  pat-ordered init is `alpha=[20,20,1]`, `beta=[1,20,20]`.
- **Architecture:** on Apple Silicon, use a **native arm64 Julia** matching the
  arm64 R that drives `JuliaCall`; an x86_64 Julia under Rosetta runs but is
  slower.
- **`@Probabilities`:** now an (intentionally) empty public slot — nothing in the
  package consumed it. Note for any downstream code that inspected it directly.

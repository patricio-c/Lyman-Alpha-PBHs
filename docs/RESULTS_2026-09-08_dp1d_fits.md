# The additive difference at z=3, with error bars — first falsifiable fit

**Run:** 2026-09-08, `stages/10_dp1d.py` on `cache_cdm.npz` / `cache_fct.npz`
(40 Mpc/h pair, 1536 sightlines, 2048 pixels, `relos.py --uniform 512
--seed 12345`), 2000 paired bootstraps, common `tau_eff = 0.42461`
(A = 1.000011 for CDM, 1.090577 for FCT).
Output: `figures/dp1d_40.{png,txt}` and `figures/dp1d_40_boot.npz`.

---

## 1. Headline

**The additive difference `ΔP1D = P1D_FCT − P1D_CDM` is not a constant, and
the rejection is not marginal.**

| model | fit | χ² / dof |
|---|---|---|
| **M0** `ΔP1D = c` — the white-noise prediction | `c = −0.1285 ± 0.0557` km/s | **1174.8 / 8 = 146.8** |
| **M2** `ΔP1D = −f·P1D_CDM + c` — drain + boost | `f = 0.2666 ± 0.0081`, `c = 4.088 ± 0.140` km/s | **93.2 / 7 = 13.3** |

`Δχ² = 1082` for one extra parameter.

The useful way to state the M0 rejection is as a required error inflation:
**the bootstrap errors would have to be 12× larger for a constant `ΔP1D` to
be an acceptable fit.** No plausible problem with the error estimate is that
large, so this conclusion stands independently of section 4 below.

**Consequence.** For a discrete PBH population the Poisson term is white noise
in 3D, which predicts `ΔP1D ≈ (A_P/4π)(k_max²−k²)`, i.e. a constant across
this window. That prediction alone does not describe this pair. Something
scale-dependent is present at low `k` and it is large.

## 2. The measurement

`k` in s/km, `ΔP1D` and `P1D` in km/s. Errors are paired bootstrap over
sightlines. Full table in `figures/dp1d_40.txt`; the fit window (DESI,
`0.00100`–`0.03185` s/km) is the first nine rows.

|  k [s/km] | k [Mpc⁻¹] |  P1D_CDM | ΔP1D | err | ΔP/err | ratio |
|---|---|---|---|---|---|---|
| 0.00140 | 0.107 | 49.944 | −12.611 | 0.672 | −18.8 | 0.7475 |
| 0.00279 | 0.214 | 42.497 | −9.499 | 0.584 | −16.3 | 0.7765 |
| 0.00419 | 0.321 | 35.403 | −5.709 | 0.523 | −10.9 | 0.8388 |
| 0.00559 | 0.428 | 33.293 | −4.735 | 0.490 | −9.7 | 0.8578 |
| 0.00768 | 0.588 | 26.379 | −2.384 | 0.292 | −8.2 | 0.9096 |
| 0.01048 | 0.802 | 23.593 | −1.675 | 0.273 | −6.1 | 0.9290 |
| 0.01467 | 1.123 | 18.105 | +0.108 | 0.162 | +0.7 | 1.0060 |
| 0.02026 | 1.551 | 13.056 | +0.837 | 0.134 | +6.3 | 1.0641 |
| 0.02724 | 2.086 | 9.382 | +1.260 | 0.089 | +14.1 | 1.1343 |

**The ratio crosses 1 at `k = 0.01429` s/km = `1.094` Mpc⁻¹.**

Two supporting numbers:

- **Pairing gain 1.62.** Unpaired resampling gives errors 1.62× larger. The
  two extractions really do share sightline positions and phases, and the
  correlation is worth a factor 1.62 on every error bar quoted here.
- **Effective slope** of `P1D_CDM` in the window: `n_eff = −0.554`, implying
  a 3D slope `n_3D = −2.554`. The ratio denominator behaves as `k^−0.55`,
  not as `k^−1`: the pair sits between the power-law and logarithmic cases,
  so neither limiting form should be assumed in the algebra.

## 3. M2 is much better and still not right

M2 is preferred over M0 overwhelmingly, but `χ²/dof = 13.3` is itself a
rejected fit, and the residuals are not random:

| k [s/km] | ΔP1D | M2 | (data−M2)/err |
|---|---|---|---|
| 0.00140 | −12.611 | −9.227 | **−5.03** |
| 0.00279 | −9.499 | −7.241 | **−3.87** |
| 0.00419 | −5.709 | −5.350 | −0.69 |
| 0.00559 | −4.735 | −4.788 | +0.11 |
| 0.00768 | −2.384 | −2.944 | +1.92 |
| 0.01048 | −1.675 | −2.202 | +1.93 |
| 0.01467 | +0.108 | −0.739 | **+5.24** |
| 0.02026 | +0.837 | +0.607 | +1.72 |
| 0.02724 | +1.260 | +1.587 | **−3.65** |

The sign pattern is `- - - + + + + + -`: three coherent runs, not scatter.
M2 is missing curvature, and errors would have to be 3.6× larger to rescue it.

M2's own pre-registered prediction also fails: from `c/f` it predicts the sign
change at `k = 0.01728` s/km against a measured `0.01429`, a 19% discrepancy
that the code flags as INCONSISTENT. That failure is not independent evidence
— it is the same bad fit seen from another angle — but it was registered in
advance and it did not come out right, so it is recorded as a failure.

**Reading.** The drain is not a constant fraction of `P1D_CDM`. That is
physically unsurprising: removing half the baryons changes the *shape* of the
gas distribution, not only its amplitude, so the suppression should be
`g(k)·P1D_CDM` with `g` scale-dependent. Assuming `g = f = const` was the
simplest testable version and it has now been tested and found wanting.

## 4. A caveat that is ours, not the data's

`--resample-A 20` was written into the stage as a self-check on the
approximation of freezing the flux rescaling `A` inside the bootstrap. The
docstring predicted `A` would move by about 0.2%, well below the ~few-percent
scatter of `P1D` itself. **It does not.** Measured:

| quantity | scatter |
|---|---|
| `A` (CDM) over resamples | 1.146% |
| `A` (FCT) over resamples | 0.866% |
| `P1D` bootstrap scatter (median over bins) | 1.271% |

`A` moves by about as much as `P1D` does, so the stated criterion — "well
below" — is **not met**, and the quoted error bars are plausibly
underestimated by some factor of order unity to a few.

Consequences, stated separately because they differ:

- **Section 1 is unaffected.** M0 needs a 12× error inflation to survive, and
  nothing here is close to 12×.
- **Section 3 is provisional.** M2 needs 3.6×. Until `A` is re-solved inside
  the bootstrap, `χ²/dof = 13.3` is not a quotable number, and neither is the
  crossing-prediction discrepancy. The *shape* of the residuals is still
  informative, because inflating errors uniformly cannot remove a coherent
  sign pattern — but its significance is not yet pinned down.

**Required fix before any χ² from this stage goes in a draft:** re-solve `A`
and recompute the per-sightline spectra inside each resample. Bracketing
`brentq` tightly around the full-sample `A` keeps this affordable at a few
hundred resamples, which is enough to measure the inflation factor even if
not enough for a full covariance.

## 5. What to do next, in order

1. **Fix the bootstrap** (section 4). Everything downstream of `χ²` waits on it.
2. **Replace the constant `f` with a measured `g(k)`.** Stage 09 now has
   `P_gas` for both runs at z=3 (`figures/pk_z3_CDM.txt`, `pk_z3_FCT.txt`,
   both passing the control against SWIFT's own `P_matter` at median ratio
   1.0000 and 1.0013). The gas power ratio supplies the shape of the drain
   with **no new free parameters** — only an amplitude. If that model fits
   where M2 fails, the two-regime picture is confirmed with the drain measured
   rather than assumed. If it does not, the low-`k` deficit is not simply the
   gas drain and the mechanism is still open.
3. **Repeat at every redshift** once the z = 5 → 2 outputs exist. One redshift
   cannot separate a drain from a boost; the evolution can.
4. **Then** `A_P` and `A_b` from the ICs, and `b(k,z)`.

## 6. What would falsify what

| claim | falsified by |
|---|---|
| `ΔP1D` is not constant | Error bars 12× larger than quoted. Re-solving `A` will not do that. |
| The low-`k` deficit is a gas drain | The `g(k)` model of step 2 fitting no better than the constant-`f` model. |
| The deficit is physical rather than numerical | `f` moving between legs (b) and (c) of the convergence suite — see `RUNS_2026-09-08`, §7.3. |
| The high-`k` excess is the primordial boost | `c` failing to scale as `k_max⁴` between resolutions, or moving with box size. |

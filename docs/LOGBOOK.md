# Logbook

**Append-only. Never rewrite an entry.**

This file records what was believed *at the time*, including the things that
later turned out to be wrong. Wrong turns are kept, not deleted: the reason a
test was run is often more useful later than its result, and an error that was
found once will be found again by whoever does not know it was already made.

Its companion is [`NARRATIVE.md`](NARRATIVE.md), which is **always rewritten**
to reflect current understanding, and which carries the list of what is still
open. If you catch yourself wanting to edit an old entry here, the edit
belongs in the narrative instead.

One entry per real step. Keep entries short and link to the detail: the stage
that did the work, the RESULTS document, the figure. Every stage log now
carries its own command, commit hash and input hashes (`common/prov.py`), so
an entry does not need to repeat them.

Strategy, paper scoping and authorship do **not** go in this file. It lives in
a public repository.

---

## Before 2026-09-04 — the sightlines were broken

SWIFT's own line-of-sight output for the 40 Mpc/h pair is unusable. The
parameter file carried `range_when_shooting_down = [0, 40]` in a box whose
internal length is 58.7372, so the rays sampled 46% of the transverse face
instead of all of it.

Every P1D number in this repository comes from sightlines regenerated with
`legacy/relos.py --uniform 512 --seed 12345`, which is 1536 rays at matched
positions in both runs.

**Still open:** when the pair is re-run with the LOS output fixed, the whole
chain has to be validated against SWIFT sightlines. Until then, `t12` cannot
test this pair, because both sides are `relos.py` output and the test would be
comparing a thing to itself.

The fix for the re-run is to delete the six optional `*_range_*` lines; the
bug was a unit error, `[0, 40]` in Mpc/h written where internal Mpc was
expected.

## Before 2026-09-04 — tau_eff fixes murgia and does NOT fix this pair

Recorded because it is the single thing most likely to be re-litigated.

In the murgia comparison, the large-scale deficit **was** an artefact of the
effective optical depth, and it went away once both runs were renormalised to
a common `tau_eff`.

**In the 40 Mpc/h CDM/FCT pair it does not.** Renormalised to the common target
0.42461 (giving A = 1.000011 and 1.090577), the tilt survives untouched, and no
choice of `tau_eff` can remove it: rescaling `tau` is one multiplicative
operation per pixel and cannot produce a `k`-dependent change to first order.

See section 7a of [`RESULTS_2026-09-04_qla_census.md`](RESULTS_2026-09-04_qla_census.md).

## 2026-09-04 — the gas census

Stages 06, 07 and 08. Where the gas goes under the QLA scheme.

- Converted gas becomes `PartType1`, verified by particle ID at 100.00%.
  `PartType4` is empty everywhere. Four independent confirmations: header
  arithmetic, stage 07, stage 08, and the t7 stars-back test.
- Conversion fractions of the initial baryons: murgia cdm 7.12%, M2 7.73%,
  M3 15.51%; **CDM40 14.08%, FCT40 49.83%**.
- The threshold is referenced to the mean **baryon** density, settled
  empirically: gas exists above Δ = 1000 and is exactly zero above Δ ≈ 3000 in
  all five runs. The residual population above the threshold is ~0.1% of the
  gas mass and is gas in flight between crossing the threshold and being
  converted on its next timestep.
- **The threshold does not preferentially eat forest gas in high-power
  models.** murgia z=10 → 5, fraction of initial baryons taken from Δ < 10:
  cdm 2.87%, M2 2.52%, M3 2.10%. It goes *down* with more small-scale power.

A sign bug in stage 06 inverted the murgia verdict for one decade (Δ = 1–10
reported +9.1% when it is −19.2%). Fixed in `59f2e0a`. The cause was clipping
signed differences with `dm[dm > 0]` before summing.

## 2026-09-08 — the additive difference, with error bars

Stage 10. The first falsifiable number in the analysis: until this point every
quantity had been quoted without an error.

- Bootstrap over sightlines, **paired** — both runs share LOS positions and
  Panphasia phases, so resampling the same index set from both keeps the
  correlation. Worth a factor **1.62** on every error bar.
- **ΔP1D is not a constant.** A flat fit gives χ²/dof = 168. The errors would
  have to be **12× larger** for the white-noise prediction to survive.
- A constant-fraction drain plus a constant boost gives f = 0.2696 ± 0.0076,
  c = 4.087 ± 0.134 km/s, χ²/dof = 11.5 — much better, still rejected, and the
  residual signs run `- - - + + + + + -`, which is a missing shape.
- The ratio crosses 1 at k = 0.01429 s/km = 1.094 Mpc⁻¹.
- Effective slope n_eff = −0.554, so the 3D slope is −2.554: between the
  power-law and logarithmic cases, and neither limit should be assumed.

Full record in [`RESULTS_2026-09-08_dp1d_fits.md`](RESULTS_2026-09-08_dp1d_fits.md).

## 2026-09-10 — a worry that was measured and turned out to be unfounded

The stage 10 docstring claimed that freezing the flux rescaling `A` inside the
bootstrap was harmless because `A` would move by ~0.2%. `--resample-A` measured
1.15% (CDM) and 0.87% (FCT) against a P1D scatter of 1.27%, and on that basis
the errors were declared possibly underestimated and every χ² provisional.

`--exact-boot 300`, which re-solves `A` and recomputes the spectra inside each
resample, gives an **error inflation of 1.00×**, range 0.90 to 1.09.

The criterion was wrong, not the approximation: the paired structure cancels
the effect, because a resample with more absorption moves `A` in both runs at
once and the difference does not notice. **χ²/dof = 11.5 is real.**

## 2026-09-10 — the 3D gas power, and what it explains

Stage 09 produced `P_gas`, `P_conv`, `P_dm`, `P_matter`, `P_baryon` and the
gas × converted cross spectrum for both runs at z = 3. The control against
SWIFT's own `P_matter` passes at median ratio 1.0000 (CDM) and 1.0013 (FCT).

Stage 11 then asked whether the P1D difference is just the 3D gas difference
propagated through `P1D(k) = (1/2π) ∫_k^∞ k' P3D(k') dk'`.

- **`g(k) = 1 − P_gas_FCT/P_gas_CDM` changes sign** near k ≈ 2.7 Mpc⁻¹: FCT has
  *less* gas power at large scales (ratio 0.806 at k = 0.15) and *more* at
  small scales (1.59 at k = 27).
- Replacing the assumed constant fraction with that measured shape takes
  χ²/dof from 11.5 to **2.42**, a Δχ² of 63.7 **at equal parameter count**.
  Substituting a measurement for an assumption, with no new freedom.
- The matter power ratio at those large scales is 1.0002. **The two models are
  identical in total matter and are not identical in gas.** The forest sees
  gas, which is why ΔP1D at large scale is −12.6 ± 0.7 km/s instead of zero.

## 2026-09-10 — corrections to the census, from the finer plots

Two earlier statements in this logbook's 2026-09-04 entry need context, which
is added here rather than by editing that entry.

**"Near-uniform drain" was based on coarse bins and is wrong at finer
binning.** The gas-mass ratio per Δ bin runs from 0.46 to 0.80 and rises
*above* 1 in the deepest voids (Δ < 0.2). That structure in Δ is what becomes
the scale dependence of `g(k)`.

**The census statement and the z=3 state are both true and are different
questions.** The census measured where converted particles *came from* — almost
none directly from Δ < 10. The z = 3 phase diagram shows the diffuse IGM is
depleted anyway. The link is accretion: gas leaves the diffuse phase, joins a
structure, and is converted there. Nothing converts at low Δ; the low-Δ
reservoir drains into the places that do.

The stage 08 maps make this visible. In cdm the converted particles sit only in
the knots of the cosmic web. In M3 they trace the entire filament network,
including thin filaments crossing voids. **This is where QLA's validation for
CDM stops transferring**: the scheme is justified for CDM because the gas it
removes is saturated and forest-irrelevant, and that premise fails when
conversion happens throughout the filaments.

## 2026-09-10 — three things stage 11 got wrong on its first run

Kept because each was found by a check that was written before the result.

**The fit-free bound was overstated.** P1D is a positive-kernel integral of
P3D, so the P1D ratio is a weighted average of the 3D ratio above that k and
must lie between its extremes. One bin fell outside — and the code declared the
point demonstrated. r_1D = 0.7765 against min r_3D = 0.7981 is a gap of 0.0216
with an error of 0.0118: **1.8σ**. The check never compared the violation to
the error. It now reports sigma and refuses to call anything below 3σ a
demonstration.

Also worth recording: the lowest P1D bin, which carries the deepest deficit,
**cannot be tested at all**, because the 3D band starts above it. That is a
box-size limit, and a 3D measurement reaching lower k would sharpen this test
considerably.

**The split-template fit was singular, provably.** Since g = g₊ + g₋ exactly,
span{T[g], 1} is a subspace of span{T[g₊], T[g₋], 1}, so χ²(M4) ≤ χ²(M3) must
hold as algebra. The run gave 25.03 > 16.97 with NaN errors from negative
variances. The excess template is nearly constant across the fit window and
degenerate with the free constant. Fixed with `pinv`, a reported condition
number, an explicit check of the nesting inequality, and M5 — the same split
without the redundant constant.

**The diagnostic panel could not have shown a problem.** It plotted P1D against
the round trip over eight decades of log axis, where the reconstruction looks
perfect regardless, while the result lives at the percent level. The median
closure was 1.0241 and was invisible. Now plotted as a ratio on a linear axis.

**And one proposed test that was discarded before running.** `P_death(Δ)` — the
conversion probability against overdensity — would be circular: stage 07
histograms against the SPH density, which is the very quantity the threshold
evaluates, so the curve is a step function at Δ = 1000 by construction. To be
informative it would have to use an independent density estimate, which is more
work than simply running the low-resolution leg.

## 2026-09-10 — the run suite, designed and costed

Ten runs from four sets of initial conditions, in
[`RUNS_2026-09-08_convergence_and_reionization.md`](RUNS_2026-09-08_convergence_and_reionization.md).

The design error worth remembering: the first version put the `z_reion` and
threshold variations on leg (a), the large box, which has the same particle
count as the production run. They belong on leg (c), an eighth of the
particles. Same science, roughly 145 core-hour-equivalents instead of 350.

Each leg is **one change** from the current setup: (c) is monofonIC grid
1024 → 512 with the box held at 40; (a) is box 40 → 80 with the grid held at
1024. Both give 8× the particle mass, which is what mass-matches them and makes
the splicing correction valid.

## 2026-09-11 — the suite generated, and four errors caught on the way

The ten parameter files and run scripts now come from
[`tools/make_lowres_suite.py`](../tools/make_lowres_suite.py) rather than from
copying one file ten times. All ten initial conditions resolve; the suite is
ready to submit.

The generator exists because the first parameter file, written by hand, had
four errors in it:

1. **The gravitational softening was copied unchanged from the production
   run.** This is the one that mattered. Softening must scale with the mean
   interparticle separation, and both new legs sit at exactly 2× the
   production spacing, so 0.0018 → 0.0036 and 0.0007 → 0.0014. Leaving it
   would have put a *force*-resolution change inside a run whose entire
   purpose is to isolate a *mass*-resolution change, and the two could not
   then have been separated. The generator derives it from box and grid
   rather than tabulating it.
2. An initial-conditions filename that does not exist.
3. `low-resolution` where the directory is `low_resolution` — one hyphen,
   and the run dies at startup.
4. A power-spectrum grid finer than the particle sampling supports, which is
   how the Fourier-cube corner modes got into the stage 11 analysis the first
   time.

**A bug inherited from the production `run.sh`:** `RESTART_FLAG` is computed
and then never used, because the `mpiexec` line carries a hardcoded `-r`.
Every run therefore started in restart mode, including the first one from
initial conditions. The generated scripts use the flag.

**And a bug of mine in the generator's first version.** It built the IC
filename from the gas particle count, when the filename encodes monofonIC's
`GridRes`: `masked` with `ParticleMaskType 2` puts gas on (GridRes/2)³, so
`IC_N1024B080` holds 512³ gas particles. Both facts are true and I conflated
them, so all ten runs reported `MISSING` against files that exist. Fixed, with
a glob fallback that *reports* the substitution rather than applying it
quietly, since a rename upstream should not silently produce ten parameter
files pointing at nothing.

**Three questions closed**, all of the kind that cannot be detected after the
fact: `zstart = 200` is correct and the z = 198 that stage 09 reads is the
first snapshot rather than the initial conditions; the FCT transfer function
was generated with the same convention as the production run; and all four
sets of initial conditions exist.

**Housekeeping:** the trailing "What is open" section was removed from this
file. It was a maintained status list living inside an append-only document,
which contradicted the file's own rule, and `NARRATIVE.md` already carries the
equivalent. One list instead of two stops them drifting. No entry was edited.

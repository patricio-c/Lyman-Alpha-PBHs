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

## 2026-09-16 — V2 against SpecWizard, and the meeting

### The extractor has an external number for the first time

`legacy/swift_extract.py` had never been compared to another code. It now has
been. Maria Marinichenko ran SpecWizard on `data/los_murgia_cdm_z5_20.hdf5` —
the 20 sightlines cut from `murgia/cdm/los_0003.hdf5` at z = 5 — and returned
tau. Our answer was already on disk in `cache/cache_murgia_cdm_z5_first100.npz`,
so the comparison needed no new extraction: her 20 line names are a subset of
those 100 and match exactly.

**Matched on tau_eff, the two codes agree to 3.5% pixel by pixel, with a
median bias of −0.5%**, over the 1087 pixels with tau < 1 in both.

That is the check `HANDOFF` has carried as pending since the beginning, and it
is a much harder test than any internal ablation.

### Why she returned 2042 pixels and not the 2048 that were asked for

Her v_box is 2730.681 km/s against our 2738.925. The cause is H(z):
SpecWizard recomputes it from Omega_m and Omega_Lambda alone and gets
554.9290 km/s/Mpc, while SWIFT's own `Cosmology/H [internal units]` in the
file header is **556.6044**, because it includes Omega_r and the massive
neutrinos. The gap is **0.3019%**, and 2730.681 / 1.337366 = 2041.84, which is
her 2042 exactly.

A 0.3% stretch of the whole k axis is small against a 19% effect and is not
small against the 0.76% error on `f`. The fix on her side is to read H from
the header. The fix on ours is that it is not needed: resampling her tau onto
our grid in **fraction of box** rather than in km/s cancels the stretch
exactly, leaving only 0.3% of the peculiar velocity, about 0.15 km/s or a
tenth of a pixel. The comparison above was done that way.

### The UVB difference behaves as a constant, and z = 5 cannot say more

Her tau_eff is 2.4459 against our 2.8877 (Gamma_HI = 4.30e-13 from
TREECOOL_HM12_G+Q). The factor that equates them is **A = 1.33507**, so if the
only difference were the ultraviolet background her table would imply
Gamma_HI = 5.74e-13 at z = 5. That number is the one to check against what the
Ploeckinger table actually gives.

Note that A = 1.335 while the raw tau_eff ratio is only 1.181. Saturation
decouples the two: at z = 5 the optical depth can be wrong by a third and
tau_eff moves by a fifth. **Matching tau_eff at z = 5 is a weak constraint**,
which is the opposite of what was assumed when this comparison was designed.

An attempt to test whether the residual carries a density dependence — which
is what would show self-shielding, dust or cosmic rays in her table acting as
something other than a rescaling — returned a slope of +0.0101 in
d ln(hers/ours) / d ln(tau). **That is not a null result, it is an
uninformative one.** The minimum tau anywhere in the file is 0.402: at z = 5
even the voids absorb. After cutting on tau < 1 in both fields the sample spans
tau from 0.65 to 0.98, a factor 1.5, and the first six bins trend upward while
the last two reverse. There is no lever arm. Worse, self-shielding acts at
Delta above ~100, which is tau far above 1, entirely inside the saturated
region the cut removes by construction. At z = 5 that question is inaccessible
in principle rather than for lack of statistics.

Three symptoms of one disease, all of which disappear at z = 3, where
tau_eff = 0.42 and tau spans orders of magnitude:

1. A moves 33.6% for 18% in tau_eff.
2. Only 2.65% of pixels are usable.
3. Those pixels cover a factor 1.5 in tau.

**What is established:** in the diffuse gas the forest actually measures at
z = 5, the two codes agree at 3.5% and the ultraviolet background difference is
a constant. What is not established is anything about the dense regime, and
z = 5 cannot establish it.

### A stale instruction, found while doing this

`HANDOFF` section 6 says to tell Maria the file is truncated when handing it
over. It is not. V1 (job 1524260) reproduced the murgia sightlines particle for
particle against regeneration from the full snapshot, and the ray positions in
the handed-over file cover 95.2% and 96.6% of the box in x and y with the
projection axis spanning all 29.521 Mpc. The `range_when_shooting_down = [0, 40]`
bug belongs to the 58.7372 Mpc box; in a 29.52 Mpc box that range clips
nothing. Her tau are physical, not merely comparable.

### The census reproduced itself from a different direction

`legacy/relosz.py --uniform 200 --seed 12345` was run on both 40 Mpc/h runs at
z = 3 to produce sightline files for the next round with Maria. Note that
`--uniform N` is rays **per axis** over three axes, so this is 600 rays, not
200 — the same convention as the `--uniform 512` that gives the production
1536.

The particle counts it printed are the conversion census, arrived at by walking
snapshots rather than by histogramming: 115318346 gas particles in CDM and
67341600 in FCT, against 512³ = 134217728 initially. That is **14.08%** and
**49.83%** converted, matching stage 06 to the decimal by an independent route.

### A risk to the small-scale half of the result

Particles per sightline: CDM median 5887, FCT median **3451**. FCT has 41%
fewer tracers per ray.

`HANDOFF` pending 1 already identified sparse SPH sampling as the leading
suspect for the residual at high k in the murgia comparison, where
`wsum_raw_med` fell from 0.924 to 0.771 and the conversion gap was 8.4 points.
Here the conversion gap is 35.7 points and the tracer gap is 41%.

Fewer tracers means a noisier SPH density, which is white noise added to the
field. White noise cannot touch low k, so **the large-scale deficit is not at
risk**. It does inflate high k, which is exactly where the other half of the
result lives — the gas ratio of 1.59 at k = 27 Mpc⁻¹ that has been read as
primordial excess.

This is not a claim that the excess is noise. It is a statement that it has not
been measured and that measuring it is cheap: the compensated-deletion test
already designed in `HANDOFF` pending 1, applied to this pair with
**f = 0.414** rather than the 0.084 written there for murgia. Delete that
fraction of gas particles from the CDM sightlines, multiply the surviving
masses by 1/(1−f) so the mean density and tau_eff do not move, and re-extract.
What comes out is the noise floor of the estimator at FCT's sampling.

### The meeting

Two objections were raised, and the tests already in this repository answer
both.

**"If the model converts 50% of the gas, that becomes stars, and the halo mass
function rules it out."** The argument does not separate the models, because it
applies to CDM with the same force: QLA-CDM converts 14.08% of the baryons
against an observed stellar mass density corresponding to a few per cent, so it
overshoots by more than the scatter as well. Reading QLA particles as stars
would invalidate the CDM Lyman-alpha literature that calibrated the
prescription. What survives is the ratio, 3.5×, which needs normalising so that
CDM matches observations before it constrains anything.

Colazo et al. (2025), A&A 702, A20, measures the halo and subhalo mass
functions for this family of models directly: the enhancement reaches a factor
of six at 10⁸–10⁹ M_sun and the high-mass end agrees with CDM and with observed
group and cluster abundances. The excess sits where there are no luminous
tracers. A follow-up paper, in preparation, performs the abundance-matching
step the objection proposes and finds that satellite stellar mass functions are
weak discriminants, because the galaxy–halo mapping reabsorbs the change in
counts; internal structure discriminates better than abundance does.

Both of those are dark-matter-only, 35 Mpc/h, at z = 0.5, with a different
blue index and PBH mass from the runs here. Same family, not the same model.
The honest statement is that the proposed test does not rule the model out, not
that it has already been run on these runs.

**"Compare against the gas spectrum and you will find the same deficit, so it
is ruled out."** The first half is this repository's 2026-09-10 result: the
deficit in P1D **is** the gas deficit, measured and propagated, Delta chi² =
63.7 at equal parameter count. The second half does not follow, because a
deficit in gas is only a prediction of the model if the model produces it. The
matter ratio at the same scale is 1.0002 while the gas ratio is 0.806, so what
differs is which particles the threshold removed, not how structure grew. And
the gas ratio is not a deficit at all: it changes sign near k ≈ 2.7 Mpc⁻¹ and
reaches 1.59 at k = 27. A physical suppression is monotonic; removing the
densest peaks is not.

What would settle it against us is the deficit surviving at threshold 3000 and
at one eighth the mass resolution. That is the suite, still unsubmitted.

### Note on this file's companion

`HANDOFF` section 0 is dated 2026-09-04 and predates stages 09, 10 and 11, the
run suite, and the error bars. Its pending 4 says nothing below k = 0.007 s/km
is quotable until the jackknife runs; it ran, and the result is in the
2026-09-08 and 2026-09-10 entries above. Its pending 6 (V2) and its pending 1
(the sampling-noise floor, with the wrong f) are likewise overtaken by this
entry. A pointer belongs at the top of that section; the section itself is not
rewritten.

## 2026-09-16 (third) — the headline corroborated, and a default that moves it

### The production numbers survive a clean re-derivation

Stage 10 on `cache/{cdm,fct}40_z3_200.npz` — 600 sightlines at ray positions
independent of the production 1536, extracted with full provenance and
`--exact-voigt`, against caches that carry no provenance at all.

    quantity   production (1536)      new (600)              agreement
    f          0.2696  +- 0.0076      0.24897 +- 0.01229     1.4 sigma
    c          4.087   +- 0.134       3.77001 +- 0.2127      1.3 sigma
    n_eff      -0.554                 -0.5545
    inflation  1.00x                  1.02x

**The reproducibility hole is closed.** Every published number for this pair
came from caches written before `common/prov.py` existed, with no record of
the flags used. They are now corroborated from a path that records everything.

Do not over-read the chi2/dof, which went from 11.5 to 3.11. That is not a
better fit: 600 sightlines against 1536 give error bars 1.6x larger, and
11.5/1.6^2 is about 4.5. M2 still fails its own sign-change test — predicted
0.01788 s/km against a measured crossing at 0.01531 — which is the same
missing shape that M3 supplies from the measured gas ratio.

### A default that silently moves a headline number

The first attempt at the above disagreed: c came out 2.96504 +- 0.186 against
the production 4.087 +- 0.134, which is **4.9 sigma**. f agreed at 1.2 sigma.

The cause was not sampling. Stage 10 defaults its common target to
`units.tau_eff_turner24`, giving 0.37187 at z = 3, while every production
number used the reference run's own raw tau_eff. Re-running with
`--tau-eff 0.42989` moved c to 3.77 and the disagreement to 1.3 sigma.

That f survives and c does not is exactly the expected signature: f is a
fraction and c is an additive amplitude in km/s, so rescaling tau by A moves
one and not the other.

`tau_eff_turner24` is the confirmed bug of the 2026-09-03 session, valid only
inside its calibration range. At z = 3 it is inside that range and is not
wrong — it is simply a **different** normalisation from the one the HANDOFF
prescribes ("the reference run's own raw value as the target, at every
redshift"). **A default that changes a headline number by 27% without saying
so is a bug even when both values are defensible.** Stage 10 should take the
reference cache's own `tau_eff_raw` as its default, or refuse to run without
an explicit `--tau-eff`. Recorded here rather than fixed, because fixing it
changes numbers and that belongs in its own step.

### Voigt is not the explanation for the SpecWizard residual

Re-extracting `data/regen_cdm_z5.0_n100.hdf5` with `--exact-voigt` and
comparing against the existing cache, which used the approximation:

    tau_eff  2.8510450 approximate   2.8510911 exact
    ln(exact/approximate) over 6234 unsaturated pixels:
        median +0.0000, sigma 0.0001  ->  0.01%

**Gap 5 of the V2 is closed.** Ours used an approximate Voigt and hers used no
profile at all, and the difference between exact and approximate on our side
is 0.01%. The 3.5% is the algorithm.

Incidentally `wsum_raw_med` on that file is 0.923, against the 0.924 the
HANDOFF records for cdm in murgia. Consistent.

### The 3.5% does not propagate into P1D where P1D is measured

Flux power ratio between the two tau fields, hers resampled onto our grid in
fraction of box and matched on tau_eff:

    k [s/km]   k [1/Mpc]   P_hers / P_ours
     0.0138       1.28         0.9998
     0.0275       2.55         1.0183
     0.0562       5.21         1.0345
     0.1101      10.22         1.0958   (5 sigma)
     0.2191      20.32         1.0217

The DESI window at z = 5 is 0.001 to 0.048 s/km. **Inside it the two codes'
P1D agree to 3%**, so the 3.5% pixel scatter does not propagate into the
measured quantity. **Gap 2 is closed where the forest is measured** and open
above it: the 9.6% at k = 0.11 s/km is a 5 sigma, non-monotonic feature that
white noise does not produce.

**The three highest-k bins are not usable, and the cause is on our side.**
Resampling with a linear interpolation is a low-pass filter whose power
attenuation is about cos^4(k dv/2): at k = 0.44 s/km that predicts 0.84 and
the bin measures 0.859. The decline at high k is the interpolation, not
SpecWizard.

This reverses an instruction given earlier the same day. It was argued that
Maria need not re-run the z = 5 file with H from the header, because
resampling cancels the grid difference. That is true for the pixel-by-pixel
comparison, which is robust to it, and false for P1D at high k, where the
interpolation sits in the middle of the measurement. The re-run is worth
asking for.

### Housekeeping that changes what path names mean

`lyman/cdm-box-40-1024` and `lyman/2-fct-box-40-1024` — the production pair —
were moved to `lyman/Old_sim/`, and the two directories previously named
`cdm-box-80-1024` and `fct-box-80-1024` were renamed to take their place, to
hold the new production runs with the line-of-sight output fixed.

**Any path in this repository or in earlier entries of this logbook that
reads `lyman/cdm-box-40-1024` refers to what is now under `lyman/Old_sim/`.**
The `relosz.py` commands in the 2026-09-16 entry above are among them. The
test runs were deliberately kept rather than deleted.

The ex-`box-80` directories are still the one thing in this analysis that was
never identified: they were named for an 80 Mpc/h box and contained a 40 Mpc/h
initial condition with production softening. They now carry production names.

### Suite checks before submitting

- **No `*_range_*` lines** in any of the ten parameter files; the generator
  documents their absence explicitly. The line-of-sight truncation bug will
  not recur.
- **Restarts are configured**: `delta_hours 3.0`, `max_run_time 11.0` under a
  12-hour wall, `resubmit_on_exit 1`. Whether a run fits in twelve hours
  stopped being a question.
- The SWIFT binary exists at the renamed path.

### Feedback: out of scope, with numbers

The binary is built `--with-subgrid=QLA`; the subgrid model is a compile-time
choice, so running with feedback is a second build with cooling and yield
tables, not a parameter change. It is out of scope for this work.

The statement that replaces it is quantitative rather than a hedge: the
threshold is a short lever, since raising it tenfold recovers only ~1.3% of
the gas; and for FCT to be consistent with the observed baryon budget at
z = 3 a model with a return channel would have to put back roughly 78% of
what QLA removed, a magnitude set by observation rather than by the subgrid
model. Both numbers are measured. Running the feedback model is not.

## 2026-09-16 (fourth) — stage 11 on the provenanced pair

Same re-derivation as the entry above, for the second headline.
`figures/drain_z3_200.txt`, from `dp1d_boot.npz` at the CDM-referenced target
and the stage 09 spectra, which are unchanged because they come from the
snapshots rather than the sightlines.

    quantity                production (1536)    new (600)
    chi2/dof M3             2.42                 0.595
    delta chi2 M2 -> M3     63.7                 17.63
    g sign change           2.7 Mpc^-1           2.725 Mpc^-1
    3D ratio range          0.806 to 1.59        0.7981 to 1.5943

The delta chi2 falls by more than the 1/1.6^2 that larger error bars predict —
24.5 expected, 17.6 measured — but M3 still beats M2 at equal parameter count
by about 4.2 sigma, and the conclusion is unchanged: the P1D difference is the
measured 3D gas difference, propagated.

### The 1.8 sigma fit-free violation was noise

The 2026-09-10 entry recorded one bin falling outside the positive-kernel
bound at 1.8 sigma, described there as suggestive rather than demonstrated.
**On independent ray positions every testable bin is inside, at n_sigma =
0.00.** It was noise. Correct this wherever the 1.8 sigma is quoted.

The lowest P1D bin still cannot be tested, because the 3D band starts at
0.1515 Mpc^-1 and that bin sits at 0.1070. Unchanged, and still a box-size
question.

**A caveat on the stage's own wording.** When every bin falls inside, the
stage prints "consistent with flux and gas changing by the same fraction".
That overstates what the bound shows. The band is [0.798, 1.594] — wide — so
falling inside is easy and is not a test of equality. The fitted amplitude is
1.88, which is 9.7 sigma from 1. Both statements are true and the sentence
invites misquoting them as one.

### No free additive constant is needed

    M3  alpha = 1.87927 +- 0.09081   c = 0.46864 +- 0.08502   chi2 = 4.16 / 7
    M5  a_sup = 1.87927 +- 0.09081   a_exc = 2.84289 +- 0.22720  chi2 = 4.16 / 7

M5 has no constant term and costs nothing against M3. **The two measured 3D
components account for the whole P1D difference on their own.** The boost c
that stage 10 fits as a free parameter — 3.77 km/s — is therefore not an extra
ingredient: it is what the measured small-scale gas excess produces. The
primordial excess does not need a parameter of its own in this fit; it is the
gas excess.

M4 remains degenerate by construction, cond(N) = 1e17, and the nesting
inequality holds exactly at 4.16 <= 4.16. Read M5.

### The two responses are different

a_sup = 1.879 +- 0.091 is the flux response to the large-scale gas deficit;
a_exc = 2.843 +- 0.227 is the response to the small-scale gas excess. They
differ by 0.964 against a combined error of 0.245: **3.9 sigma.** The forest
does not respond to the two halves of g(k) with one amplitude, which is the
whole reason the split templates exist.

That both are of order 2 is the order expected when the neutral fraction goes
as the square of the density, but the value the fluctuating Gunn-Peterson
approximation actually predicts has not been computed here. Suggestive, not a
result.

### Still a missing shape

Residual signs: M2 `---+++++-`, M3 and M5 `-++++-++-`. Better than a single
constant fraction, and still a run of four. chi2/dof of 0.595 on 9 bins is on
the low side (p about 0.76), so the error bars may be slightly generous. The
sign run is the thing to chase, not the chi2.
\n
## 2026-09-16 (fifth) — the sampling gap, what sets the drain, and two false alarms

Four things measured after the entries above, recorded together.

### The sampling gap, measured (this belongs chronologically before the third entry)

`stages/01_extract_los.py` on the 600-ray z = 3 files. The extractor's own
diagnostics say more than the particle counts did.

    wsum_raw_med       CDM 0.8548        FCT 0.5014
    frac_below_floor   CDM 0 (max 0.002) FCT 0.0024 (max 0.019)

`wsum_raw_med` is the SPH partition of unity before the Shepard correction:
how much of the kernel weight finds particles. **Half of it falls on nothing in
FCT.** The murgia case that `HANDOFF` pending 1 was written for showed
0.924 -> 0.771; this pair shows 0.855 -> 0.501. And the Shepard floor
intervened on up to 1.9% of FCT pixels against essentially never in CDM: those
are pixels with too little gas to interpolate at all. Shepard renormalisation
fixes the mean and not the variance.

**This makes the compensated-deletion test calibratable.** The criterion for
the test being correctly set up is that deleting a fraction f of the CDM gas
particles brings CDM's `wsum_raw_med` down to ~0.50. Start at f = 0.414, the
tracer gap, and tune against wsum rather than assuming it.

Two consistency checks passed. Raw tau_eff from 600 independent ray positions
is 0.42989 (CDM) and 0.40426 (FCT), against 0.42461 and 0.40350 from the
production 1536-ray caches at different positions -- 1.2% and 0.2%. Ray
placement is not biasing anything, and FCT is less opaque, as it must be if
the threshold took more of its gas.

Grid for those files, for the record and for the SpecWizard comparison:
2048 px, v_box = 4497.9099 km/s, dv = 2.196245 km/s,
H(z=3) = 306.3077 km/s/Mpc from `Cosmology/H`, Gamma_HI = 8.27580e-13 (HM12).

### Where the SpecWizard comparison stands

**Settled:** in the diffuse gas at z = 5 the two codes agree to 3.5% with a
-0.5% median bias over 1087 pixels; the 2042-vs-2048 grid difference is fully
accounted for by the H(z) convention; the flux power ratio is 1.00 to 1.03
inside the DESI window, so the pixel scatter does not propagate into what is
measured; and the Voigt profile is not the explanation, at 0.01%.

**Not settled, in order of weight:** FCT, since the check is CDM only and the
wsum numbers above say the extractor works in a different regime in each box;
z = 3, where every production number lives; and the dense regime, which z = 5
cannot reach at all.

### What sets the drain: sigma at 10^8 M_sun, not the Nyquist scale

Twelve runs in `lyman/more_power` -- six `NB` (broken primordial spectrum,
f_pbh = 0) and six `poisson` (Poisson term, no break) -- share Panphasia phases
with production and differ only in the transfer function. Conversion fractions
come from header arithmetic; the transfer functions give the linear field at
z = 200 exactly.

Binning the density field on a scale R and asking which R makes both families
fall on one curve gives **R = 0.067 Mpc/h, i.e. M = 1.05e8 M_sun/h**:

    conv(z=7) = 3.06% x (sigma_R / sigma_R,CDM)^3.19     rms 7.6%

Zero rank inversions across twelve models spanning a factor 7.4 in conversion
(2.98% to 22.16%). Two different mechanisms, one relation.

**The obvious label does not work.** Ranked by Delta^2 at the Nyquist scale the
rank correlation with conversion is only +0.46, and the failure is stark:
`poisson_1` has the most Nyquist power of all twelve (Delta^2 = 0.0500, at the
2LPT ceiling) and converts 3.41%, while `NB_1` has seven times less (0.0069)
and converts 10.12%. The three NB models with identical Delta^2 = 0.0291
convert 22.16%, 14.53% and 6.22% -- a factor 3.6 at fixed Nyquist power.

Power below ~10^8 M_sun does not convert gas, because those scales do not
collapse. That is the scale where Colazo et al. (2025) find the FCT subhalo
excess, which is a consistency rather than a derivation.

**Production is outside the grid and the relation cannot be extrapolated.** At
the same R, sigma_FCT/sigma_CDM for the production model is **6.17**, against a
grid spanning 1.00 to 1.94. In raw power, P_FCT/P_CDM is 8 at k = 10 h/Mpc,
285 at k = 40, and 1870 at k = 80. The power law predicts over 1000%
conversion there, which is the fit announcing it has left its regime. The grid
calibrates the mechanism; it does not predict production. A CDM run at the
grid's own `sc` particle load would anchor the family on a measured zero
instead of the `poisson_6` proxy used here.

### False alarm 1: the 2LPT criterion

Raised and withdrawn the same day, kept because the wrong number was quoted
first.

Production FCT (`IC_lyman_alpha/6_10^6/B40`) starts at z = 200 with
`ParticleLoad = masked` and `GridRes = 1024` in a 40 Mpc/h box. Computing
Delta^2 at the **monofonIC grid** Nyquist of 80.4 h/Mpc gives 1.63, which is
33 times the Delta^2 <= 0.05 criterion used to design the twelve-model grid.
That number was quoted as a possible invalidation of the production run.

**It is the wrong k.** With a masked load the grid is not what the simulation
represents. Gas sits on 512^3 and dark matter on 4 x 512^3, so:

    k_Nyq [h/Mpc]              Delta^2 FCT
    80.4   monofonIC grid         1.63     <- not the relevant scale
    63.8   dark matter particles  0.83
    40.2   gas particles          0.205

Four times over at the gas Nyquist and seventeen at the dark matter one, not
thirty-three. `IC_lyman_alpha/Validate_ic.py` makes exactly this argument in
its docstring and it is the standard one.

Two further points, and the second closes it. Reproducing the input P(k) --
which is what `Validate_ic.py` checks -- is necessary but not sufficient,
because a 2LPT field with Delta^2 ~ 1 still reproduces the linear P(k) at
z_start; what degrades is the higher-order displacement, which the power
spectrum cannot see. But **2LPT transients decay as a^-2**, so from z = 200 to
z = 3 the error is suppressed by (201/4)^2 ~ 2500, and at Delta^2 between 0.2
and 0.8 there is no shell crossing in the initial conditions. The concern does
not survive that.

Closed unless a direct measurement of Delta^2 from the initial-condition
particle distribution says otherwise.

### False alarm 2, and one real problem, from the same logs

Job `1534358`, reported as `FAILED` with exit code 1 after 27 seconds, was
**the dry run**, and it passed: *"Time integration ready to start. End of
dry-run."* The exit code is the usual MPI complaint at teardown. It was
briefly treated as a startup failure in a queued run.

The dry run did confirm the suite is sound: threshold 3000 applied where
intended, the correct initial conditions read, 16777216 gas and 67108864 dark
matter particles, softening 0.003600 / 0.001400 Plummer, all three output
lists read, restarts every three hours, and `--line-of-sight` and `--power`
both in the engine policies. No `*_range_*` lines anywhere in the ten
parameter files, so the line-of-sight truncation bug will not recur.

**The real problem is in the leg (a) initial conditions, and both are
affected.** `slurm-1531635.out` (FCT) and the corresponding CDM job both ended
`CANCELLED ... DUE TO TIME LIMIT` at 20:07 and 20:12. The FCT log shows all
eight ranks writing `PartType1` and **only two of eight writing `PartType0`**
before the kill. The two leg (c) jobs, eight times smaller, completed in 2:39
and ended `Done writing SWIFT IC file`.

So `IC_N1024B080_200_{cdm,fct}.hdf5` are truncated in their gas. Both were
regenerated with a longer wall time, and the SWIFT jobs waiting on them were
held. The check that they came out whole is the header `NumPart_Total`
against the actual dataset lengths, plus gas masses that are not zero.

### A note on the cosmology, settled

The dry run reports `Omega_k = 0.001389` and `N_nu = 0`. The production z = 3
line-of-sight file reports `Omega_k = 0.00138918`, `N_nu = 0`: identical, so
the resolution test is clean. The origin is exact -- the monofonIC config uses
`Omega_m = 0.306` with `m_nu1 = 0.06`, while the SWIFT parameter files declare
`Omega_cdm + Omega_b = 0.3046109` with no neutrino block, and
0.306 - 0.3046109 = 0.0013891. The initial conditions carry a massive neutrino
that the runs do not, and SWIFT puts the difference into curvature. It cancels
in every ratio; it matters only if absolute P1D is ever quoted against a
Planck baseline that includes 0.06 eV.

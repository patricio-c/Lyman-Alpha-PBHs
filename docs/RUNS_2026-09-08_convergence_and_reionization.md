# Planned simulation suite — convergence and reionization

**Written:** 2026-09-08. **Revised:** 2026-09-10 with the costed suite.
**Status:** configurations settled, not yet submitted.

This file exists so the design is not re-derived from scratch in a later session.
Every run below is required by the project plan; none is exploratory.

---

## 1. Why these runs exist

Two items in the plan cannot be answered by post-processing the existing cache,
and therefore need new SWIFT runs:

1. **Resolution and box convergence.** The broken-spectrum component
   (`n_b = 2`) contributes to `ΔP1D` as `k_max^{n_b+2} = k_max⁴`. It is by far
   the most fragile number in the analysis: a 10% error in the effective cutoff
   is a 46% error in that term. Nothing from it is quotable until convergence
   is demonstrated.
2. **The `z_reion` derivative.** Reionization is the one degenerate partner of
   the PBH boost that is deliberately inside the current scope, and the only
   one that cannot be obtained by repainting sightlines.

A third question turned out to need a run as well: how much of the low-`k`
deficit depends on the star-formation threshold itself. It is included below.

Everything else — mean flux, `T₀`, `γ` — is pure post-processing on the LOS
cache and must **not** consume compute.

---

## 2. Convergence: the splicing triad

Large boxes and high resolution cannot be afforded simultaneously. The standard
resolution (Borde et al. 2014; Lukić et al. 2015) is to run three boxes and
combine them:

```
P(large box, high resolution)  ≈  P(a) × P(b) / P(c)
```

The correction is only valid if **(a) and (c) share the same particle mass**.
That is the one constraint that cannot be relaxed.

| leg | box [Mpc/h] | gas grid | m_gas (rel.) | status |
|-----|-------------|----------|--------------|--------|
| **(b)** small box, high res | 40 | 512³ | 1 | **already run** |
| **(a)** large box, low res  | 80 | 512³ | 8 | to run |
| **(c)** small box, low res  | 40 | 256³ | 8 | to run |

`(80/512)³ = (40/256)³` exactly, so (a) and (c) are mass-matched and the ratio
isolates the finite-box effect cleanly. (b) sits in the same box as (c) with 8×
better mass resolution, so `P(b)/P(c)` isolates the resolution effect.

### The IC recipe: one change per leg

Since `m_gas ∝ Ω_b ρ_c L³ / N_gas` and `N_gas ∝ grid³`, each leg is a single
change from the current setup:

| leg | change | result |
|---|---|---|
| **(c)** | monofonIC grid **1024 → 512**, box stays at **40** | gas 256³, mass × 8 |
| **(a)** | box **40 → 80**, grid stays at **1024** | gas 512³, mass × 8 |

The scaling holds whatever `masked=2` does internally, so this does not depend
on knowing the masking details.

**The box must stay at 40 for leg (c).** Halving the box *and* the grid keeps
the mass resolution unchanged, which is a box-size test rather than a splicing
leg, and it breaks the mass matching that makes the correction valid.

### Pre-flight check before submitting

The existing run is named `40-1024` but its header reports
`Ngas = 134,217,728 = 512³`. The label is **not** the gas grid. Confirm
`Ngas` in each IC header before submitting: 512³ in the 80 Mpc/h box and
256³ in the small low-resolution box.

---

## 3. The suite: four sets of ICs, ten runs

All variations run on leg **(c)**, which has an eighth of the particles of leg
(a) and is being run anyway. Leg (a) has the same particle count as the current
production run, so every copy of it costs a full run.

| runs | box / grid | `H_reion_z` | `over_density` | purpose |
|---|---|---|---|---|
| (a) CDM, (a) FCT | 80 / 1024 | 7.5 | 1000 | large-box leg of the splice |
| (c) CDM, (c) FCT | 40 / 512 | **7.5** | 1000 | small-box leg **and** the central reionization point |
| (c) CDM, (c) FCT | 40 / 512 | **6** | 1000 | reionization derivative |
| (c) CDM, (c) FCT | 40 / 512 | **9** | 1000 | reionization derivative |
| (c) CDM, (c) FCT | 40 / 512 | 7.5 | **3000** | threshold sensitivity |

**Ten runs, two of them expensive, from four sets of ICs.** Every variation
reuses the leg-(a) or leg-(c) initial conditions unchanged; only the SWIFT
parameter file differs.

### Why not put the variations on leg (a)

Because it costs four times the compute for the same science. Seven copies of
leg (a) is roughly 350 core-hour-equivalents against roughly 145 for the layout
above. The variations measure *responses*, and a response does not need the
large box.

### Every variation runs in both models

Including the threshold. Varying it only in FCT would mean that when the result
moves, there is no way to tell a model-specific sensitivity from a generic
property of the scheme. The whole analysis is built on ratios between the two
models, so a one-sided change confounds the observable itself. The control
costs one cheap run.

### The values, and why these

**`z_reion` = 6, 7.5, 9.** These bracket what is observationally allowed:
Planck puts the midpoint near z ≈ 7.7 and the forest ends reionization near
z ≈ 5.3–6. A derivative wants points close enough for the response to be
linear and inside the range that matters.

A model-consistent early value for FCT is a *different* and defensible run,
because FCT forms structure — and therefore ionizing sources — much earlier.
But that number must be **measured, not chosen**: the converted-mass history
of both runs gives the collapsed fraction against redshift, and the epoch
follows from it. A value picked by hand does not survive a referee.

**`over_density` = 3000, not 10000.** QLA exists precisely to remove dense gas
before it forms cold clumps that collapse the timestep. Without feedback, gas
allowed to reach Δ ~ 10⁴ risks a run that crawls. A factor of three is a real
lever, and the census already shows that nothing survives above Δ ≈ 3000 at the
current threshold, so 3000 is where the gas will actually pile up. Go to 10⁴
only if 3000 shows a large effect, and then knowing what to expect.

Note that this is **not** made redundant by the post-hoc test on the existing
snapshot, which found that raising the threshold tenfold recovers only ~1.3% of
the deficit. That test cannot capture what a re-run does: gas that is not
converted keeps evolving, accreting and affecting its surroundings.

### Stated caveat on the derivative

`∂P1D/∂z_reion` is measured at the particle mass of the low-resolution legs.
Whether it transfers to leg (b)'s resolution is **plausible but untested** —
responses are generally more robust than absolute values, but that is an
expectation, not a measurement. It is written here rather than left implicit,
and testing it would cost one more copy of leg (a), which is not worth it now.

---

## 4. What each run answers, and what would falsify it

| run | question | falsified if |
|---|---|---|
| (a), (c) | Is the `k_max⁴` broken-spectrum term converged in the 40/512³ box? | The spliced `P(a)P(b)/P(c)` differs from `P(b)` by more than the quoted error on `ΔP1D`. Then no number from that term is quotable at current resolution. |
| (b) vs (c) | Is the drain fraction `f` physical or a sampling artefact? | `f` moves between the two. Same box, same phases, only the mass resolution differs, so a shift means the deficit is numerical. **This is the cheapest decisive test in the suite.** |
| (a) vs (b) | Does the low-k deficit survive a larger box? | The deficit shrinks with box size → it is a missing-large-scale-power artefact, not the baryon drain. |
| `z_reion` grid | Can reionization timing mimic the FCT/CDM tilt? | `∂ ln P1D/∂z_reion` reproduces the measured tilt shape within errors. Then the two effects are degenerate at `z=3` and the whole case rests on redshift evolution. |
| `z_reion` grid, CDM vs FCT | Is the reionization response model-independent? | The two derivatives agree → `z_reion` can be marginalized as an independent nuisance, which is the *easy* case. |
| threshold grid | Is the deficit an artefact of where the threshold sits? | The deficit moves substantially between Δ = 1000 and 3000 in FCT but not in CDM → the result is a property of the subgrid choice, not of the model. |

---

## 5. Output lists

Every run uses the **same output list as the production run** (z = 5 to 2 in
steps of 0.1, plus the snapshot list already settled). The splicing correction
is only needed at z = 3 to answer the convergence question, but these runs are
the cheapest that will exist and the entire analysis is about redshift
evolution. Making the correction available at every z costs nothing here and
cannot be added later without re-running.

---

## 6. Explicitly not in this suite

- **Full hydrodynamics with feedback.** Needed eventually to show the deficit is
  not a QLA artefact. The gas census bounds it partially, but the maps from
  stage 08 show that in the high-power model conversion happens throughout the
  filament network rather than only in the knots, which is exactly where QLA's
  validation for CDM stops transferring. Deferred, and this suite is what
  decides how urgent it is.
- **An extreme-FCT model.** Worth running to test whether small-scale power can
  be pushed until it is unambiguously detectable, but the model has not been
  chosen yet. Do not generate ICs until the target is pinned down.
- **Anything varying mean flux, `T₀` or `γ`.** These are repaints of the existing
  cache. Spending compute on them is a mistake.

---

## 7. Coordination

The `z_reion` variation overlaps with a Sherwood-style suite (multiple
resolutions, multiple slopes) that is wanted independently — same runs, two
purposes. Settle that overlap before submitting; it is a scheduling decision,
not a code change.

---

## 8. What to do with the output

The runs are worthless without the analysis being fixed in advance. This is
that analysis, written before the data exist so it cannot be tuned to them.

### 8.1 Order of operations — this order, not another

1. Extract sightlines on **matched positions within each pair** (same
   `--uniform N --seed`), so the CDM and FCT members of a box can be
   bootstrapped as paired samples. Pairing does not transfer across boxes;
   only within one.
2. **Renormalise every run to a common `tau_eff` before anything else.**
   Splicing multiplies ratios of power spectra; if the mean flux differs
   between the legs, that difference propagates into the spliced curve as a
   spurious scale dependence.
3. Splice, then difference. **Not** the reverse — see 8.2.
4. Run stage 10, then stage 11, on the spliced pair.

### 8.2 Two ways to splice, and the one that has to be checked

The correction can be applied either to each model separately,

```
P1D_CDM_spliced = P1D_CDM(a) x P1D_CDM(b) / P1D_CDM(c)      (then likewise FCT)
```

or directly to the ratio,

```
R_spliced = R(a) x R(b) / R(c)
```

These are algebraically identical and numerically are not: the second lets
correlated errors between CDM and FCT cancel inside each leg, and should be
the tighter estimate. **Compute both and report the difference.** If they
disagree by more than the bootstrap error, the splicing correction is not in
its linear regime and neither number is usable.

### 8.3 The sharp test the new runs make possible

Stage 10 fits `DP1D = -f P1D_CDM + c`, and stage 11 replaces the assumed shape
with the measured 3D gas suppression. Whichever parametrisation is used, the
drain amplitude and the boost amplitude behave differently under a change of
resolution:

- **The boost** carries the broken-spectrum term, which scales as `k_max^4`. It
  is *expected* to move between (c) and (b), and the size of that move is the
  convergence result.
- **The drain** is supposed to be a statement about how much absorbing gas the
  star-formation scheme removes. If that is physical, it should be **stable**
  between (b) and (c), because the two boxes are the same volume with the same
  phases.

| outcome | reading |
|---|---|
| drain stable | The drain is a property of the scheme acting on the gas, not of the sampling. It is a real effect to be reported, and the low-`k` deficit becomes a prediction of the model. |
| drain moves with resolution | The deficit is a resolution artefact of how the threshold interacts with the SPH sampling. Then the deficit is an upper bound on what a model with feedback would do, not a prediction, and the signal reduces to the small-scale excess. |

Both outcomes are publishable and they change what is claimed, not whether
there is a result.

### 8.4 The reionization derivative

For each model separately, finite-difference the three `H_reion_z` values:

```
J_reion(k) = d ln P1D(k) / d z_reion
```

Use the three-point central difference at 7.5 and report the two one-sided
differences as well; if they disagree badly the response is not linear over
this range and a wider or finer grid is needed before it can be marginalised.

Then compute the quantity the whole scope question turns on:

```
Delta J(k) = J_reion^FCT(k) - J_reion^CDM(k)
```

- `Delta J` consistent with zero → reionization is a **separable** nuisance.
  It can be marginalised independently of the signal, which is the easy case
  and should be stated plainly as such.
- `Delta J` non-zero → the reionization response depends on the model being
  tested, so it cannot be marginalised as an external systematic. Every
  constraint must then be derived on a joint grid, and saying so is itself a
  result.

Finally, compare `J_reion(k)` in shape to the measured `ln(P1D_FCT/P1D_CDM)`.
If reionization timing can reproduce that shape within errors, the two are
degenerate at a single redshift and the entire case rests on redshift
evolution — which is what the z = 5 to 2 output list exists to measure.

### 8.5 What must be recorded for every run

Non-negotiable, because past unit mismatches make it so: box in Mpc/h *and*
internal units, gas grid, particle masses, `H_reion_z`, `over_density`, the
`tau_eff` used and the `A` it implied, the sightline seed and count, and the
git hash of the extraction. Since `common/prov.py` was added, every stage log
carries its own command, commit and input hashes; the rest belongs in the run
directory next to the parameter file.

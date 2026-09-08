# Planned simulation suite — convergence and reionization

**Written:** 2026-09-08. **Status:** configurations settled, not yet submitted.

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

Each leg is run for **CDM and FCT** → 4 new runs.

### Pre-flight check before generating ICs

The existing run is named `40-1024` but its header reports
`Ngas = 134,217,728 = 512³`. The label is therefore **not** the gas grid.
Confirm the actual gas grid each configuration produces before submitting;
what matters is 512³ gas in the 80 Mpc/h box and 256³ gas in the small
low-resolution box, whatever the directory ends up being called.

### Cost

- **(a) ×2** — same particle count as the current run in a larger box. Comparable
  wall time, possibly slightly less (less collapse, longer timesteps).
- **(c) ×2** — one eighth of the particles. Substantially cheaper.

The suite is *not* four runs of equal cost: two are comparable to the current
run and two are cheap.

---

## 3. Reionization: fold it into the cheap leg

Rather than commissioning a separate small box, vary `H_reion_z` on the **(c)**
configuration, which is the cheapest available and is being run anyway.

| `H_reion_z` | CDM | FCT | note |
|---|---|---|---|
| 6.0 | ✓ | ✓ | |
| **7.5** | ✓ | ✓ | **this pair *is* leg (c)** |
| 9.0 | ✓ | ✓ | |

Six cheap runs deliver both deliverables: leg (c) of the splicing triad comes
for free as the central case, and the remaining four give
`∂P1D/∂z_reion` **separately for CDM and for FCT**.

That separation is the point. If the derivative differs between the two models,
`z_reion` is not an additive systematic that can be corrected independently — it
is coupled to the very signal being measured. A difference is physically
expected, since FCT forms structure, and therefore ionizing sources, earlier.

### Resolves a standing inconsistency

The murgia comparison runs use `H_reion_z: 9.0` while the 40 Mpc/h pair uses
`7.5`. Once the derivative is measured, that difference stops being an
uncontrolled discrepancy between suites and becomes a quantified correction.

---

## 4. What each run answers, and what would falsify it

| run | question | falsified if |
|---|---|---|
| (a), (c) | Is the `k_max⁴` broken-spectrum term converged in the 40/512³ box? | The spliced `P(a)P(b)/P(c)` differs from `P(b)` by more than the quoted error on `ΔP1D`. Then no number from that term is quotable at current resolution. |
| (a) vs (b) | Does the low-k deficit survive a larger box? | The deficit shrinks with box size → it is a missing-large-scale-power artefact, not the baryon drain. |
| `z_reion` grid | Can reionization timing mimic the FCT/CDM tilt? | `∂ ln P1D/∂z_reion` reproduces the measured tilt shape within errors. Then the two effects are degenerate at `z=3` and the whole case rests on redshift evolution. |
| `z_reion` grid, CDM vs FCT | Is the reionization response model-independent? | The two derivatives agree → `z_reion` can be marginalized as an independent nuisance, which is the *easy* case. |

---

## 5. Explicitly not in this suite

- **Full hydrodynamics with feedback.** Needed eventually to show the deficit is
  not a QLA artefact, but the gas census already bounds it: QLA drains baryons
  near-uniformly in overdensity rather than preferentially from forest gas.
  Deferred.
- **An extreme-FCT model.** Worth running to test whether small-scale power can
  be pushed until it is unambiguously detectable, but the model has not been
  chosen yet. Do not generate ICs until the target is pinned down.
- **Anything varying mean flux, `T₀` or `γ`.** These are repaints of the existing
  cache. Spending compute on them is a mistake.

---

## 6. Coordination

The `z_reion` variation overlaps with a Sherwood-style suite (multiple
resolutions, multiple slopes) that is wanted independently — same runs, two
purposes. Settle that overlap before submitting; it is a scheduling decision,
not a code change.

---

## 7. What to do with the output

The runs are worthless without the analysis being fixed in advance. This is
that analysis, written before the data exist so it cannot be tuned to them.

### 7.1 Order of operations — this order, not another

1. Extract sightlines on **matched positions within each pair** (same
   `--uniform N --seed`), so the CDM and FCT members of a box can be
   bootstrapped as paired samples. Pairing does not transfer across boxes;
   only within one.
2. **Renormalise every run to a common `tau_eff` before anything else.**
   Splicing multiplies ratios of power spectra; if the mean flux differs
   between the legs, that difference propagates into the spliced curve as a
   spurious scale dependence.
3. Splice, then difference. **Not** the reverse — see 7.2.
4. Run stage 10 on the spliced pair.

### 7.2 Two ways to splice, and the one that has to be checked

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

### 7.3 The sharp test the new runs make possible

Stage 10 fits `DP1D = -f P1D_CDM + c`. The two coefficients have different
expected behaviour under a change of resolution:

- **`c` is the boost.** It carries the broken-spectrum term, which scales as
  `k_max^4`. It is *expected* to move between (c) and (b), and the size of
  that move is the convergence result.
- **`f` is the drain fraction.** It is supposed to be a statement about how
  much absorbing gas the star-formation scheme removes. If that is physical,
  `f` should be **stable** between (b) and (c), because the two boxes are the
  same volume with the same phases.

So: **fit `f` on (b) and on (c) and compare.**

| outcome | reading |
|---|---|
| `f(b) = f(c)` within errors | The drain is a property of the scheme acting on the gas, not of the sampling. It is a real effect to be reported, and the ratio approach survives. |
| `f` moves with resolution | The low-k deficit is a resolution artefact of how the threshold interacts with the SPH sampling. Then the ~19% number is not quotable and the two-regime reading collapses. |

This is the cheapest decisive test in the whole suite and it comes free with
leg (c).

### 7.4 The reionization derivative

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

### 7.5 What must be recorded for every run

Non-negotiable, because past unit mismatches make it so: box in Mpc/h *and*
internal units, gas grid, particle masses, `H_reion_z`, the `tau_eff` used and
the `A` it implied, the sightline seed and count, and the git hash of the
extraction. Stage 10 already writes most of this into its `_boot.npz`; the
rest belongs in the run directory next to the parameter file.

# What is going on, and why it matters

**This file is always rewritten to reflect current understanding.** It carries
no history — for that, and for the things that turned out to be wrong, see
[`LOGBOOK.md`](LOGBOOK.md), which is append-only.

Last rewritten: 2026-09-10.

---

## The question

Does extra primordial power on small scales leave a signature in the Lyman-α
forest that nothing else can imitate?

The difficulty is old and well understood: at a single redshift the forest
cannot separate extra small-scale power from a colder, less pressure-smoothed
intergalactic medium. Both push the flux power spectrum the same way. The
standard answer is to impose physical priors on the thermal history. The answer
being developed here is to use the *redshift evolution* instead, because the
competing effects evolve differently — a fixed comoving excess grows
gravitationally, while thermal smoothing is progressively erased.

## The experiment

Two SWIFT simulations, 40 Mpc/h on a side, gas on a 512³ grid with four times
as many dark matter particles, run to z = 3.

They share their initial phases. The **only** difference is the primordial
power spectrum: one is CDM, the other adds a broken spectrum plus the Poisson
isocurvature of a discrete population of primordial black holes.

That construction matters. Because the phases are shared, every structure sits
in the same place in both boxes, so the difference between them is the effect
of the model rather than the noise of a different realisation. And because both
extra terms vanish at large scales — the break sits at k_t = 10 Mpc⁻¹ and the
Poisson term contributes as k³ — **the two boxes are supposed to be identical
where the forest is measured best.**

They are, in matter. The matter power ratio at k = 0.14 Mpc⁻¹ is **1.0002**.

## The result that needs explaining

The flux power spectra are not identical there. Renormalised to a common
effective optical depth, so that no part of the difference is a difference in
mean flux, the additive difference at the largest measured scale is

> **ΔP1D = −12.6 ± 0.7 km s⁻¹**, eighteen sigma from the zero it should be,

with the sign reversing at k = 0.0143 s/km and the ratio rising to well above 1
at small scales, where the primordial excess does live.

The large-scale half of that is the part that should not exist.

## Why it exists: the forest does not see matter, it sees gas

At the same large scale where the matter power ratio is 1.0002, the **gas**
power ratio is **0.806**.

The chain is this. More small-scale power means structure collapses earlier and
in more places. More gas therefore crosses the density threshold of the
star-formation scheme and is converted into collisionless particles: **49.8% of
the initial baryons in the high-power run against 14.1% in CDM**, leaving gas
mass fractions of 0.082 and 0.138 at z = 3.

Less absorbing gas, distributed differently. And the distribution matters as
much as the amount: the gas that is removed sits in the densest peaks, which
are the most strongly biased tracers of the large-scale field, so removing them
lowers the bias of everything that survives. That is a large-scale signature
produced entirely by small-scale physics, carried by a process — star formation
— that is non-linear and non-local.

The maps from stage 08 show it directly. In CDM the converted particles sit
only in the knots of the cosmic web. In the high-power model they trace the
entire filament network, including thin filaments crossing voids.

## The difference in P1D is the difference in gas, propagated

This is measured rather than assumed. Stage 09 gives the 3D gas power of both
runs; stage 11 propagates the measured difference through

```
P1D(k) = (1/2π) ∫_k^∞ k' P3D(k') dk'
```

and compares it to the flux difference. Replacing an assumed constant
suppression with the measured shape improves the fit by **Δχ² = 63.7 at equal
parameter count** — a measurement substituted for an assumption, with no new
freedom — taking χ²/dof from 11.5 to 2.42.

The measured suppression **changes sign** near k ≈ 2.7 Mpc⁻¹: less gas power at
large scales, more at small scales. Both halves of the flux signature come from
one measured curve.

There is also a constraint that needs no fit at all. Since P1D is a
positive-kernel integral of P3D, the ratio of two P1D curves is a weighted
average of the ratio of their 3D curves, and therefore lies between that
ratio's extremes above every k. Two things follow: a deficit in P1D **requires**
a deficit in 3D and cannot be manufactured by how the integral weights small
scales; and if the measured ratio falls outside the allowed band, flux power
and gas power are not changed by the same fraction. At present the violation is
1.8σ in one bin — suggestive, not a demonstration — and the most informative
bin cannot be tested because the 3D measurement does not reach low enough k.

## What is not yet known, and what decides it

**Is the drain physical, or is it the scheme?**

The star-formation scheme used here is a fast prescription designed to make the
forest cheap, not to model star formation. It is validated for CDM, and the
validation rests on a premise: the gas it removes is dense, saturated, and
irrelevant to the forest anyway. The maps show that premise failing in the
high-power model, where conversion happens throughout the filament network.

That both models convert far more mass than the observed stellar mass density
allows — by more than an order of magnitude, in *both* cases — is the direct
statement that this is a gas-removal prescription rather than star formation.
The absolute numbers therefore cannot be read as a stellar mass prediction.
What survives is the **ratio: 3.5× more converted mass** in the high-power run,
which is a consequence of the primordial spectrum and is testable against
observations once normalised so that CDM matches them.

**One cheap run decides it.** Re-running the same box with the same phases at
one eighth the mass resolution isolates the sampling: if the drain amplitude is
stable, the deficit is a property of the scheme acting on the gas and the
large-scale signature is a prediction. If it moves, the deficit is a resolution
artefact and becomes an *upper bound* on what a model with feedback would do,
leaving the small-scale excess as the signal.

Both outcomes are results. They change what is claimed, not whether there is
something to claim.

## What is being built

The pieces in place are the estimator and its error bars, the 3D fields, the
propagation from 3D to 1D, and the census that explains the mechanism.

What is missing is the part the whole idea rests on: **redshift evolution**.
Everything above is a single snapshot at z = 3. One redshift cannot separate a
drain from a boost, because at fixed z any smooth amplitude change can be
absorbed. The evolution can, because a fixed comoving excess and a thermal
effect are different functions of z — and because the observational window is
fixed in velocity units while the excess is fixed in comoving wavenumber, so
the two slide against each other by about 25% between z = 2 and z = 4. That is
a change in shape, and no bin-by-bin amplitude rescaling reproduces it.

The remaining ingredients are the amplitudes of the two primordial terms, which
are calculable from the initial conditions rather than fitted; the reionization
response, measured separately for both models so that it can be shown to be a
separable nuisance or not; and a fixed observational window common to every
redshift bin, so that measured evolution contains no instrumental component.

# Where the gas goes: the QLA conversion census

**2026-09-04.** Self-contained record of one session's work, referenced from
section 0 of `HANDOFF.md`. Everything here is measured on snapshots that
were already on disk. No simulation was re-run.

The question this session answered, in Joop's words: *where is the gas I am
losing, and in which haloes.* Underneath it sat a sharper methodological
worry — **does the QLA star-formation threshold bias models with extra
small-scale power against the forest?** If it does, the FCT/CDM P1D
comparison measures a subgrid artefact rather than physics.

---

## 0. The short version

1. **Converted gas becomes dark matter, not stars.** Verified by ID, not
   inferred: every doomed `ParticleID` turns up in `PartType1`, 100.00% in
   all three murgia runs. `PartType4` holds N = 0 in every snapshot.
2. **The threshold does NOT preferentially eat forest gas in the
   high-power models.** In murgia, M3 converts twice the gas cdm does and
   still removes *less* mass from forest densities in absolute terms:
   2.10% of the IC baryons against cdm's 2.87%. The methodological worry
   is not supported by the data; it points the other way.
3. **FCT drains 40% of the whole baryon reservoir, near-uniformly in
   density.** That does suppress the forest - but by emptying the
   reservoir, not by any threshold bias.
4. **The FCT material was marked in the initial conditions.** Its doomed
   gas starts out 2.4x more overdense at z = 198 than CDM's does. That is
   the PBH Poisson isocurvature seeding the loss, measured directly.
5. **The total matter power of the two 40 Mpc/h runs agrees to 1.0002 at
   the fundamental mode**, and FCT exceeds CDM at every smaller scale. So
   the entire P1D suppression lives in the baryon partition, not in the
   dynamics.
6. Therefore: **the 21% P1D suppression is an upper bound set by a subgrid
   model with no feedback.** The direction is physical; the magnitude is
   not yet trustworthy. `t7` brackets it from the other side.

---

## 1. A bug that inverted one published verdict

`stages/06_gas_census.py` attributed the missing gas per decade with
`dm[dm > 0]`, which silently discarded every bin where the comparison run
holds MORE gas than the reference, then normalised by the sum of the
positives. Consequences:

- the printed "total missing mass" was not the mass actually missing:
  murgia M3 reported 0.11678 against a true net loss of 0.08655, and the
  FCT pair 0.41017 against 0.40418;
- the `Delta = 1-10` decade of murgia M3 - **the forest** - is net
  POSITIVE for the model. The old table called it +9.1% of the missing
  mass. The correct answer is **-19.2%**: M3 holds about 1.7% of cdm's
  total gas mass MORE than cdm does, at exactly the densities that make
  the forest.

Fixed in commit `59f2e0a`. Shares are now signed fractions of the net
loss and the run warns when the comparison run gains mass anywhere.

The same commit removed a claim that was wrong in a subtler way. The line
`within a factor 3 above the threshold` was documented as deciding whether
raising the threshold would help. **It cannot.** Conversion is a sink, so
the standing occupancy above the threshold is near zero in every run by
construction - those bins are censored and the number comes out near zero
whatever the answer is. It was replaced by an explicit extrapolation of
the uncensored decade ratio, which self-diagnoses by printing observed
against predicted in the last uncensored decade.

**Do not repeat this class of mistake.** Any statistic that clips a
signed difference, and any statistic evaluated inside a censored region,
will produce a confident number that means nothing.

---

## 2. Conversion accounting, from the headers alone

`gas + DM` is invariant under QLA conversion. Both totals are exact in
every snapshot: murgia `268,435,456 = 2 x 512^3`, the 40 Mpc/h boxes
`671,088,640 = 512^3 x 5` (gas on a 512^3 grid, 4 DM per gas particle).
So the converted count is the gas deficit against the IC grid, exactly.

| run | z | converted | fraction of IC baryons |
|---|---|---|---|
| murgia cdm | 5 | 9,553,042 | 7.12% |
| murgia M2 | 5 | 10,379,229 | 7.73% |
| murgia M3 | 5 | 20,814,668 | 15.51% |
| CDM 40 | 3 | 18,899,382 | 14.08% |
| FCT 40 | 3 | 66,876,128 | 49.83% |

Mass tracks count to 0.4-2%, so particle splitting is a minor effect and
counts and masses may be quoted interchangeably at that precision.

**The conversion histories are opposite.** Cumulative converted fraction:

| | z=10 | z=7 | z=5 | z=3 |
|---|---|---|---|---|
| murgia cdm | 0.37% | 2.88% | 7.12% | - |
| murgia M3 | 4.92% | 11.41% | 15.51% | - |
| CDM 40 | - | 1.57% | 5.28% | 14.08% |
| FCT 40 | - | **36.62%** | 43.14% | 49.83% |

FCT converts three quarters of its total before z = 7. CDM converts most
of its own after z = 7. Any epoch-pair analysis has to account for that,
which is why `stages/07` reports its own coverage and warns below half.

---

## 3. The threshold reference density, settled empirically

`used_parameters.yml` confirms `QLAStarFormation:over_density: 1000` in
every run, but not what it is referenced against - mean baryon density or
mean total matter, a factor `Omega_m/Omega_b = 6.33` apart.

The census settles it without reading SWIFT's source. Normalising Delta by
the cosmic mean **baryon** density, gas exists above Delta = 1000
(0.00094 of the total in murgia cdm, 0.00068 in CDM at z=3) and is
**exactly zero above Delta = 3000 in all five runs**. Were the criterion
referenced to total matter, the cut-off would fall at Delta_b = 6330 and
gas would survive at 3000. It does not.

**The criterion is baryon-referenced.** Every Delta quoted in stages 06,
07 and 08 stands as printed.

---

## 4. Stage 06, corrected: where the gas is missing

Signed shares of the net loss, per decade in Delta.

**murgia, cdm - M3 at z = 5.** Net loss 0.08655 of cdm's gas mass.

| decade | signed share |
|---|---|
| < 1 | +95.7% |
| 1 - 10 | **-19.2%** |
| 10 - 100 | +20.3% |
| 100 - 1000 | +3.0% |
| > 1000 | +0.2% |

The ladder is coherent: voids drain, the forest band fills from below, and
everything above Delta ~ 10 is eaten by the sink. **The forest band gains
gas in M3.**

**FCT, CDM - FCT at z = 3.** Net loss 0.40418.

| decade | signed share | CDM's own share of its gas |
|---|---|---|
| < 1 | +37.4% | 32.0% |
| 1 - 10 | +45.6% | 47.4% |
| 10 - 100 | +14.0% | 16.0% |
| 100 - 1000 | +3.0% | 4.6% |

The deficit is distributed almost exactly as CDM distributes its own gas.
**This is a near-uniform 40% drain of the reservoir, not selective removal
at any density.** Relative to its own surviving gas FCT is slightly more
concentrated (29.3 / 48.2 / 16.8 / 5.69 against 32.0 / 47.3 / 15.9 /
4.63), so the structural push is there; it is simply swamped.

### Would a higher threshold put the gas back? No.

Extrapolating the uncensored decade ratio at z = 3: CDM would retain
0.01648 of its gas with a threshold ten times higher, FCT 0.01135. The
difference is 0.00513, **1.3% of the 40% deficit.** Raising the threshold
is not the lever. (At z = 5 the decade ratio is not geometric - observed
exceeds predicted - so the murgia extrapolation is not usable and the run
says so.)

---

## 5. Stage 07: where the doomed gas was, while it was still gas

New. Matches `ParticleIDs` between two snapshots: present in `PartType0`
early, absent late, means converted in between - exact, no statistical
matching, and splitting only ever adds IDs at the later time.

**murgia z = 10 -> z = 5** (95 / 92 / 68% of each run's conversion):

| run | median Delta | 90th pct | % > 100 | % < 10 | -> % of IC baryons taken from Delta < 10 |
|---|---|---|---|---|---|
| cdm | 11.995 | 105.19 | 10.88% | 43.67% | **2.87%** |
| M2 | 16.158 | 173.14 | 15.46% | 36.13% | **2.52%** |
| M3 | 47.213 | 524.30 | 38.77% | 20.32% | **2.10%** |

**This is the answer to the methodological worry.** M3 converts twice the
gas cdm does and takes *less* mass out of forest densities. Its extra
conversion is entirely gas that was already at Delta > 10; its doomed gas
sat four times denser than cdm's. The FCT pair at z = 7 -> z = 3 says the
same: CDM takes 5.11% of the IC baryons from Delta < 10, FCT 2.50%.

Independent confirmation of section 4: the forest band gains in M3 because
the threshold is not eating it.

**FCT z = 198 -> z = 3, 100% coverage.** At z = 198 nothing is dense, so
the question changes: does the material that later converts start out
overdense? Separation between the doomed and surviving medians:

| run | doomed | survives | separation in delta |
|---|---|---|---|
| CDM | 1.0274 | 1.0002 | **0.0272** |
| FCT | 1.0705 | 1.0064 | **0.0641** |

**Factor 2.4.** The material FCT eventually converts was already marked in
the initial conditions - the PBH Poisson isocurvature, measured directly
rather than inferred. CDM's value agrees across two independent binnings
(1.0272, 1.0273) and is solid. FCT's is a **lower bound**: its 90th
percentile is 1.2143 and the fine run clips at Delta = 1.4, so quote the
ratio as 2.4-2.7 until it is measured over the full range. Stage 07 does
not yet report mass falling outside the histogram; it should.

---

## 6. Stage 08: where the converted particles ended up

New. `PartType1` carries no `Densities`, so this builds a CIC field of all
the matter at the late snapshot and samples it at the converted particles'
positions, in units of that field's own mean. The statistic is the ratio
of the converted median Delta to a random `PartType1` sample measured
identically - shot noise affects both alike, so the ratio is robust and
the raw medians are not.

**murgia, converted between z = 10 and z = 5:**

| grid | cell | cdm ratio | M3 ratio | gap |
|---|---|---|---|---|
| 256^3 | 0.11533 Mpc/h | 27.7 | 7.27 | 3.8x |
| 512^3 | 0.05767 Mpc/h | 85.1 | 11.98 | **7.1x** |

**M3's converted particles are far less clustered than cdm's**, and
refining the grid *widens* the gap. A smoothing artefact would close it,
so the effect is physical and coarse grids were hiding it. It is **not
converged** - both medians rose by a factor 4 when the grid doubled - so
the sign and direction may be quoted, the numbers may not.

Note the resolution floor honestly: at 256^3 the cell is 28 kpc physical
at z = 5 while a 10^9 Msun halo has R_vir ~ 4 kpc. **The grid measures
environment on ~100 kpc/h scales, not haloes.** "Which haloes" needs
HBT-HERONS + SOAP on murgia; a grid cannot answer it at these masses.

**Restricting to late converters** (z = 7 -> z = 5, 256^3) closes part of
the gap: cdm 25.7, M3 11.80, gap 2.2x against 3.8x for the full set. cdm
barely moves; M3 rises 62%. So M3's early-converted population is the
dispersed one - conversion epoch explains roughly 40% of the effect.

Two live hypotheses for the rest, **not yet decided**:

- **(a) timing.** M3 converts early in small haloes that are later merged
  or tidally stripped, spreading their particles. Supported by the z=7
  test, partially.
- **(d) halo mass.** M3 converts in many low-mass haloes sitting in
  ordinary filaments, while cdm converts in the few massive haloes that
  occupy dense nodes. Untested; needs a halo catalogue.

**Ruled out: (c) resolution**, since the gap grows with grid refinement.
**Ruled out: "QLA eats the halo core, leaving it underdense."** Conversion
removes no mass - the particle stays in `PartType1` and keeps depositing
onto the grid - so it cannot lower the *total matter* density anywhere.
That mechanism would leave its signature in the gas field, which stage 06
panel (b) shows, not in the matter field stage 08 measures.

**This does not reopen the forest systematic.** Section 5 settled that at
the epoch the gas was last gas. Stage 08 describes where the material sits
afterwards as dark matter, and dark matter produces no absorption either
way.

---

## 7. The measurement that ties it together

SWIFT already wrote `power_spectra/` in every run directory;
`requested_spectra: ["matter-matter", "cdm-cdm"]`. Index `0010` is z = 3.

Use `matter-matter` and **never** `cdm-cdm` for cross-run comparison: in a
QLA run `PartType1` is DM plus the converted baryons, which is 18.9M
particles in CDM and 66.9M in FCT - a different physical quantity in each
run. Total matter is identical in composition regardless of the partition,
and its shot noise matches too, since conversion changes no particle's
mass. So only `matter-matter` divides cleanly.

**P_matter(FCT) / P_matter(CDM) at z = 3:**

| k [Mpc^-1] | ratio |
|---|---|
| 0.14 | **1.0002** |
| 0.30 | 1.0007 |
| 1.00 | 1.0104 |
| 3.00 | 1.0746 |
| 10.0 | 1.2855 |
| 20.0 | 1.6526 |
| 40.0 | 2.2846 |

Two things follow.

**The ICs are validated non-linearly.** At the fundamental mode the runs
agree to four decimals. The shared Panphasia phases did their job, and no
large-scale difference in the P1D can come from the matter field. The
excess is purely small-scale and monotonic, exactly as a broken spectrum
plus a k^3 Poisson term must behave.

**And the decisive one.** Over the P1D window, k = 0.47 to 7.3 Mpc^-1,
FCT's total matter has MORE power than CDM everywhere, up to +20%. The
flux P1D has **21% less** (0.7857 against ref-CDM). The matter says more
clustered, the flux says less power. **The entire suppression lives in the
baryon partition.**

Quote nothing above k ~ 50 Mpc^-1: Nyquist for 1024^3 in 58.7 Mpc is 55,
and the 366 bins run to 2492 only because they include the foldings.

### Why an ad-hoc k-rescaling of the P1D must not be used

It was proposed, and it should not be. Three reasons:

1. There is no tilt to correct. Against the correct tau_eff target the
   FCT/CDM ratio is 6.79-7.65% across the band, and **there are no error
   bars yet** (pending: the t9 jackknife). Fitting a free linear-in-k term
   to 0.86 points of scatter of unknown significance is fitting noise.
2. The matter ratio *rises* with k while the flux ratio is flat. There is
   no common k-filter to remove.
3. The 21% is the result, not a nuisance. A free k-dependent factor that
   makes the curves agree deletes the measurement and makes every later
   comparison circular.

The one legitimate k-dependent correction is the SPH sampling-noise floor
(pending item 1, `f = 0.084`, mass-compensated): white, flat in k, and its
amplitude is **derived from the particle sampling, not fitted**.

---

## 8. What is now safe to say, and what is not

**Established.**

- Converted gas becomes `PartType1`, verified by ID, 100.00%.
- The QLA threshold does not preferentially remove forest gas in
  high-power models; it removes less of it, absolutely.
- FCT's suppression comes from a near-uniform 40% drain of the baryon
  reservoir.
- Raising the threshold tenfold recovers ~1.3% of that deficit.
- The threshold is referenced to the mean baryon density.
- The ICs reproduce CDM at the fundamental mode to 1.0002 at z = 3.
- FCT's doomed material is 2.4-2.7x more overdense than CDM's in the ICs.

**Not established, do not quote.**

- Any converged number for the stage 08 clustering ratio.
- Which haloes hold the converted particles. The grid cannot answer it.
- Whether (a) timing or (d) halo mass dominates the residual gap.
- Whether the 0.86-point scatter in the flux ratio is significant.

**The headline for the paper.** FCT's extra small-scale power drives more
gas over the QLA threshold, and QLA, having no feedback, removes it
permanently. The 21% P1D suppression is therefore an **upper bound**; t7,
which puts the converted mass back, gives the lower one. The direction is
physical and measured three independent ways. The magnitude is set by a
subgrid model and must be presented as bracketed, not as a prediction.

---

## 9. Next, in order

1. **t9 jackknife error bars.** Nothing about k-dependence can be settled
   without them, and they gate the question above.
2. **t7**, the converted mass returned to the forest - the lower bracket
   on the 21%. Use a fixed mass cut, `--conv-mass-max 2.0e-04` in internal
   units: gas is ~1.2e-04, DM ~6.3e-04, and 2.0e-04 sits safely in the
   gap. **Do not use the value `stages/00_inspect_snapshot.py --deep`
   prints for M3.** Its verdict averages the extremes of a top-12
   population list, and M3's masses are near-continuous so the converted
   particles do not make that list; it returned 6.317e-04, which is inside
   the DM distribution and would classify half the dark matter as
   converted gas. cdm's 3.776e-04 happened to be safe. This is a bug in
   stage 00 and it is still there.
3. **P_gas,3D with Pylians on `PartType0`**, since SWIFT did not write a
   gas spectrum. The one measurement of scale dependence with no free
   normalisation at all.
4. **HBT-HERONS + SOAP on murgia**, to decide (a) against (d) and to
   answer "which haloes" properly.
5. Stage 07 should report the mass falling outside its histogram range, so
   a clipped tail cannot pass silently as it nearly did for FCT.

## Falsifiers

- **Section 5's verdict** dies if, at an epoch with better coverage, M3
  removes more mass from Delta < 10 than cdm in absolute terms.
- **Section 7's decisive claim** dies if P_gas,3D(FCT)/P_gas,3D(CDM) turns
  out to tilt across the P1D window; it predicts a near-flat ratio times
  the baryon fraction.
- **The stage 08 gap** dies if it closes at 1024^3. It widened from 256 to
  512, so this is unlikely, but it is not converged.

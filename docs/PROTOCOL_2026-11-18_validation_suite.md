# Protocol, 2026-11-18 — the validation suite, and the track for SpecWizard

This document says what is measured next, with which run, by which command,
and what would falsify each result. It also defines the **track**: the exact
inputs and steps by which Maria Marinichenko's SpecWizard produces the flux
numbers that are finally quoted.

## The division of labour

Our pipeline (stage 01 → stage 10) computes every number below, and has been
validated: Murgia+2019 reproduced, and SpecWizard agreement at 3.5% pixel by
pixel and 3% in P1D inside the DESI window at z = 5 (LOGBOOK, 2026-09-16).
**The flux numbers that end up quoted come from SpecWizard**, run on exactly
the same sightline files, and passed through exactly the same downstream
stages. Only the extractor changes. Two independent extractors on identical
inputs, every step in a stage script, every input listed with a checksum:
anyone can re-run it, and anyone can check it on GitHub.

The 3D quantities (stages 13, 14, 15) do not go through any extractor, so they
need no second code.

## Measurements

### M1 — does the new pair carry the same FCT problem?

Every number up to 2026-10-14 is the **old pair** with regenerated sightlines.
The new pair (`cdm40`/`fct40`) is the same model on the same phases, re-run
with native line-of-sight output.

| quantity | stage | command | old-pair value |
|---|---|---|---|
| gas fraction and r_gg, D = r_gg/r_mm, every z | 14 | `sbatch scripts/run_census.sh` | r_gg(0.15 Mpc⁻¹, z=3) = 0.8055 |
| IC marking of converted gas | 13 | `CDM=… FCT=… OUT=figures/qla_initial_overdensity_new40 sbatch scripts/run_initial_overdensity.sh` | 2.652 |
| flux ratio and ΔP1D, z = 3, 1500 random sightlines, seeds 1 and 2 | 01 + 10 | `run_extract_random.sh` then `run_pair_ratios.sh` | 0.7475 ± 0.011 at k = 0.0014 s/km; ΔP1D = −12.6 ± 0.7 km/s |
| the same at z = 2.2 and 4.0 | 01 + 10 | same, `Z=2.2`, `Z=4.0` | none: first measurement |

**Falsified if** the new-pair flux ratio at the largest scale is consistent
with 1, or r_gg(0.15) is consistent with r_mm. Either would mean the old
pair's deficit came from the regenerated rays and not from the gas.

### M2 — the threshold

`lr_delta_*` (Δ = 3000) against `lr_zreio75_*` (Δ = 1000), same phases,
same 512³ resolution, both models. Stage 14 gives r_gg and f_gas for
Δ3000/Δ1000 within each model, stage 15 the converted fraction, and stage 10
the flux.

**Limitation, stated up front:** at 512³ the threshold barely bites (τ_eff
equal to 0.2% on 768 sightlines, 2026-10-18). A null result here says nothing
about 1024³. Splitting pressure smoothing from the threshold at production
resolution needs a **1024³ Δ = 3000 pair**, not yet run: two runs of the size
of the production pair, 4 nodes each, about two weeks of wall clock with
resubmissions. That decision is pending.

**Falsified if** the FCT/CDM drain D at Δ = 3000 equals the one at Δ = 1000
at 1024³: then the threshold value is not what sets it.

### M3 — the response to extra power, ordered by integrated extra power

Stage 15 gives each initial spectrum one number for "how much extra", from
the `input_powerspec.txt` that monofonIC actually used, against the same CDM:

- `dI1D(k_par)`: (1/2π) ∫ k (P − P_CDM) dk from k_par to the cut, divided by the
  CDM value. This is what the P1D kernel feeds into the flux at k_par.
  Evaluated at k_par = 0.001, 0.003, 0.01 and 0.03 s/km at z = 3.
- `dA_lnk`: the extra area in Δ²(k), ∫ (Δ² − Δ²_CDM) dln k, relative.
- `sig_ratio`: σ_R/σ_R,CDM at R = 0.067 Mpc/h (10⁸ M☉/h), the scale on
  which conversion collapsed both grid families onto one relation
  (2026-09-16).

Each is computed with two cuts: the gas Nyquist of the run itself, and the
smallest gas Nyquist in the suite, which is common to every run so that
different resolutions can be compared.

The response, per model: the gas fraction at every z (stage 14, from
`statistics.txt`), D(k, z) (stage 14), and the flux ratio (stage 10).
The question is which ordering puts the twelve grid runs, the production FCT
and the two low-resolution FCTs on one monotonic curve, and which does not.

**The grid needs its own CDM.** The grid runs use `ParticleLoad = sc`, which
puts the gas on 1024³ particles. CDM b40 uses `masked`, with gas on 512³. The
grid is therefore at eight times **better** gas mass resolution than
production. The first census measured the size of that difference:
`poisson_6` has almost no extra power (P/P_CDM = 1.01 at k = 30 h/Mpc). It
converts 16.6% of its gas by z = 3 against 14.1% in `cdm40`, and sits at
D(0.15) = 0.935 against it. Until the anchor exists, every grid ratio against
`cdm40` carries that offset. The converted fraction does not: it is absolute
for each run.

The anchor, `more_power/CDM_sc`, is the grid's `conf.conf` with **one line
changed**: `transfer_file = tf_z/tf_200_v_cdm_CLASS.csv`, the TF of CDM b40.
That file and the `grilla12/*_camb.csv` files have the same 13 columns on the
same k grid. Their density and velocity columns scale together, so the
`_CLASS`/`_camb` suffixes do not mark different conventions between the two.
The SWIFT parameter file is NB_1's with the name and IC path changed, and the
binary is the same build (md5 `f9b47a5b…`, identical in all twelve grid
directories and the production pair).

    IC:     IC_lyman_alpha/more_power/CDM_sc/      job 1580446
    SWIFT:  lyman/more_power/CDM_sc/               job 1580447 (after the IC)

**Falsified if** no ordering is monotonic: then conversion and the drain do
not depend on one scalar of the linear spectrum, and the grid cannot be read
as a calibration.

### M4 — reionization redshift

Three z_reio (6, 7.5, 9) in both models at 512³. Two things are measured:

1. the **nuisance direction**: CDM(z_reio)/CDM(7.5) and FCT(z_reio)/FCT(7.5),
   in r_gg (stage 14) and in the flux (stage 10);
2. **separability**: whether FCT/CDM at fixed z_reio is the same for all
   three.

**Falsified (separability) if** FCT/CDM changes with z_reio by more than the
error bars: then the reionization response and the extra-power signal do not
simply add, and must be modelled jointly.

### M5 — after M1 to M4: the inputs to b(k,z)

`tests/t11_bias_bk.py` needs the amplitudes of the Poisson and broken-spectrum
terms and the effective cut-off, all calculable from the production transfer
function. With M1 measured at three redshifts, b(k,z) can be computed. It is
written once M1 to M4 have settled what enters it.

## The track for SpecWizard

### Inputs

Every sightline sample is written by `scripts/run_extract_random.sh`, which
appends one row to **`figures/los_manifest.tsv`** per sample:

    tag  run_dir  z_req  los_file  z_file  subset  N  seed  sha256  cache  git

`subset` is a valid SWIFT LOS file (`scripts/make_los_subset.py` copies every
root group verbatim and corrects the counts), on Clementina under
`Lyman-Alpha-PBHs/cache/los_subsets/`. `sha256` (first 16 hex) identifies the
exact file. The same seed gives the same sightline indices in both models of
a pair. `stages/02_check_los_match.py` says whether they are also the same
lines in space, and stage 10 is run paired or unpaired accordingly.

### Conventions, so that her tau and ours are on the same grid

- 2048 pixels over the box, per sightline, in the order of the subset file.
- H(z) from `Cosmology/H` in the file header, which includes Ω_r and the
  massive neutrino. Recomputing it from Ω_m and Ω_Λ stretches the k axis by
  0.30% (2026-09-16).
- The ionizing background does not need to match ours (HM12, TREECOOL). Every
  ratio is normalised to the **reference CDM's own raw τ_eff**, and a constant
  rescaling of τ drops out of it.

### What comes back

One file per subset, `tau[n_los, 2048]` plus the sightline names in subset
order, into `external/specwizard/<tag>.hdf5`. A converter into the stage 01
cache format, `tools/specwizard_to_cache.py`, **is not written yet**; it is
the next code to write. After it, the final numbers are

    python stages/10_dp1d.py external/specwizard/cache_<ref>.npz \
        external/specwizard/cache_<test>.npz --tau-eff <ref raw tau_eff> ...

which is exactly what `scripts/run_pair_ratios.sh` runs on our caches.

### Acceptance

Our P1D and hers, on the same subset, agree within 3% inside the DESI window,
as at z = 5. A larger difference on FCT and not on CDM is itself a result: it
would mean the extractors disagree where the gas has been drained.

### Order of priority

1. `cdm40`, `fct40`, z = 3, seed 1: the headline.
2. `lr_zreio75_cdm`, `lr_zreio75_fct`: the resolution test.
3. `lr_delta_*`, `lr_zreio6_*`, `lr_zreio9_*`, `lr_a_*`.
4. The twelve grid runs and `CDM_sc`, once the anchor exists.
5. `cdm40`, `fct40` at z = 2.2 and 4.0.

## Jobs launched 2026-11-18

    census (stages 14, 15)                  1580482, stage 15 re-run 1580485
    extract, new pair z=3 s1 / s2           1580426 / 1580427
    extract, new pair z=2.2 / z=4.0         1580428 / 1580429
    extract, ten low-resolution runs z=3    1580430
    extract, twelve grid runs z=3           1580431
    ratios (afterok on the extractions)     1580432 .. 1580436
    stage 13: new pair / lr75 / lr delta    1580437 / 1580438 / 1580439
    CDM_sc anchor: IC / SWIFT               1580446 / 1580447

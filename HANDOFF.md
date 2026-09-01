# HANDOFF — read this first

You are picking up a scientific analysis repository. This file is the whole
context. Read it before touching anything.

The person you are working with is **Pato Colazo**, PhD student in astronomy
at IATE/OAC (Universidad Nacional de Córdoba, CONICET). He works in
Rioplatense Spanish and writes papers in English. He is on a Linux laptop;
the cluster is remote.

---

## 0. Current status — read this first, every session

Updated 2026-09-01 (second pass, after V1). This section is the fast-changing
punch list; the rest of the file is stable background. Start every new session
here, and dip into the numbered sections below only for the *why* behind
something.

**Done:**
- Repo scaffolding pushed to `https://github.com/patricio-c/Lyman-Alpha-PBHs`
  (public). Migration of `legacy/` done 2026-08-31, `VERIFY OK`.
- Geometry check passed clean 2026-08-31. Every number in section 5 matched.
- Validation block A run 2026-09-01 (`scripts/run_validation_A.sh`, one log).
- `common/units.py` numpy-2 fix: `getattr(np, "trapezoid", np.trapz)`
  evaluated its default eagerly and died on numpy >= 2.3, where `np.trapz` no
  longer exists. `tests/test_estimator.py` now runs to the end.
- **V1 PASSED, 2026-09-01.** See below. This was the gate on t7.
- `legacy/relos.py` fixed to write `Xaxis`/`Yaxis`/`Zaxis`/`NumParts`. See
  below. Ray positions are untouched.
- The V2 hand-off file is cut: `data/los_murgia_cdm_z5_20.hdf5`, 20
  sightlines, 16.9 MB, from `murgia/cdm/los_0003.hdf5` at z=5.0. Verified to
  carry `NumParts`, `Xaxis`, `Xpos`, `Yaxis`, `Ypos`, `Zaxis` per group.

---

### V1 — do we control the sightlines? **PASSED.**

`sbatch scripts/sbatch_roundtrip.sh .../lyman/murgia/cdm 5.0 100`, job 1524260.
`figures/t12_roundtrip_cdm_z5.0_n100.txt`.

    ray positions identical                    max offset 0.000e+00
    particle sets identical                    97/100 sightlines
    particles only in SWIFT, missed by us      0
    particles only in ours                     3
    worst relative field difference            0.000e+00
    pixels where tau differs                   864 of 204800, all in LOS 26/40/52
    max |dF|                                   5.96e-06
    d<F>                                       3.73e-09
    tau_eff                                    2.851045 both
    P1D ratio                                  1.0000 at every k

t12 prints FAIL because its criterion is binary (`exact == n`). **Read the
numbers, not the verdict.** We never miss a particle SWIFT had, the field
values on shared particles are bit-identical, and the three extra particles
sit in pixels with `tau ~ 4000`, i.e. `F = 0` on both sides. The flux field
and the P1D are the same to numerical precision.

The sentence this licenses:

> Our regeneration of SWIFT's line-of-sight output reproduces the flux field
> and the 1D flux power spectrum to numerical precision. The particle sets
> agree on 97 of 100 sightlines; the three exceptions differ by one particle
> each, in fully saturated pixels, and change P1D by less than 1e-4 per cent.

**Scope of what was tested**, so nobody over-reads it: one run, one redshift,
100 sightlines, **all along axis 2** (the first 100 groups of a murgia file
are all axis 2), reusing the template's ray positions. Not tested: the
`--uniform` path, the `--skip-los` batching, and the axis detection in
`read_rays()` on axes 0 and 1. Enough for t7, which reuses positions.

**The three extra particles are not understood and it does not matter.**
`b/(gamma*h)` = 0.495, 0.879, 0.766, so not a kernel-edge rounding case.
`SplitCounts = 0` and progenitor == own ID, so not particle splitting.
`Header/Time` is bit-identical in both files, so not a different epoch. Their
`h` is in the bottom percentile of their sightline, so not a stale cell
`h_max`. Three hypotheses tried, three wrong. If it is ever worth chasing,
the answer is in SWIFT's `line_of_sight.c`, not in more numpy. It changes
nothing downstream.

---

### The `uni512_seed12345` finding — this resolves three open items at once

**The production caches were NOT made from SWIFT's own LOS files.** From the
`source` field of the caches themselves:

    cache/cache_cdm.npz  (1536, 2048)  regen/cdm40_z3.0_uni512_seed12345.hdf5
    cache/cache_fct.npz  (1536, 2048)  regen/fct40_z3.0_uni512_seed12345.hdf5

`relos.py --uniform 512` draws 512 rays per axis over three axes = **1536**,
and `make_uniform_rays` draws uniformly over the **whole** box face. Verified
bit-for-bit: regenerating with `--uniform 512 --seed 12345` on the 58.73715
box gives `LOS_0000  Xaxis=1 Xpos=13.353070342258142 Yaxis=2
Ypos=18.605482517633362 Zaxis=0`, identical to the file on COSMA.

Three things that had been logged as separate mysteries are the same fact:

1. **"6144 or 1536 sightlines?" — answered, and it was a false dichotomy.**
   They are different files. SWIFT's own files hold 6144. The analysis used a
   1536-ray `--uniform 512 --seed 12345` regeneration. Nothing was subsampled.
2. **The published numbers do NOT rest on truncated sightlines.** The 6144
   SWIFT files are truncated to 46% of the transverse face by the `h` mix-up.
   The 1536 uniform rays cover the whole face by construction. Any earlier
   note in this file implying the production results are truncated is wrong.
3. **`tau_eff` 0.42461 (cache) vs 0.26192 (100-line cut) is not a bug.** The
   cut came from the truncated SWIFT file, where each sightline is missing
   ~32% of its particles. Fewer particles, lower density, lower tau.

**And it changes what stage 02 measured.** `stages/02_check_los_match.py` was
pointed at `lyman/cdm-box-40-1024/los_0010.hdf5` and
`lyman/2-fct-box-40-1024/los_0010.hdf5`, i.e. **the 6144 SWIFT files, not the
files the caches came from.** Those two SWIFT files are genuinely unrelated to
each other. That says nothing about the 1536 uniform rays, which — same seed,
same box, same `default_rng` — are the *same positions in both runs*.

That in turn explains t9's result, which had looked like a contradiction:

    sigma(common index) / sigma(unpaired)   0.3733
    per-line correlation between runs       +0.8592

Under a common-index bootstrap draw, `corr(X_A, X_B)` is exactly the per-line
Pearson correlation in the population. +0.86 is impossible for unrelated
sightlines and expected for matched ones. **The lines are very probably
paired after all, and t8 is very probably valid.** Re-run stage 02 on the
`regen/*_uni512_seed12345.hdf5` pair to close it. Do that before requoting any
significance, and before acting on anything this file previously said about
"the pairing is broken".

---

### Correction to section 1

The table in section 1 says the `tau_eff` rescaling "makes it slightly
deeper". **It makes it shallower.** `t0_rescaling` at z=3: `r(0.003)` goes
from 0.7338 raw to 0.7797 rescaled. Fixed in the table.

Related, and worth a sentence in the paper: the rescaling is **not flat in
k**. It moves `r(k)` by +6.26% at k=0.003 s/km and +3.91% at k=0.03, a 2.35
percentage-point tilt across the window, and it flattens the slope of `r(k)`
by 2.2%. That is a systematic to quote, not an effect to deny. The mechanism:
`A` multiplies `tau`, not `delta_F`, and `|dF/dlnA| = A*tau*exp(-A*tau)`
vanishes at both `tau << 1` and `tau >> 1`. So `A` is a spatially selective
gain, saturated pixels do not respond, and since saturation traces
large-scale structure the gain is k-dependent. After rescaling, FCT has 21.1x
more transparent pixels (F>0.99) than CDM — the two runs lean on the F=1 wall
by different amounts.

---

### `legacy/relos.py` was edited — deliberately, see rule 3.2

`make_uniform_rays` wrote only `Xpos` and `Ypos`. SpectWizard needs `Xaxis`,
`Yaxis` and `Zaxis` and aborts without them, which is why a `fix.py` exists on
COSMA. Separately, `write_output` copied `NumParts` verbatim from the template
LOS file, so a regenerated file carried the *original's* particle count, not
its own.

Both fixed. **No number changes:** our extractor reads `Xpos`/`Ypos` (the log
prints `ray position source: ['attrs']`) and never reads `NumParts`, and the
RNG stream is untouched, so `--uniform N --seed S` returns exactly the same
positions as before. Verified against the COSMA file, digit for digit.

Anything regenerated from now on is directly readable by SpectWizard without
post-processing. Files made before this commit still need
`legacy/fix_los_attrs.py`, which is the script that was patching the
attributes in after the fact and which stays for exactly that reason — the
`regen/*_uni512_seed12345.hdf5` files the published numbers came from were
written before this fix.

---

### V3 — reframed onto murgia, and it is now an external check

V3 was "is the deficit real": t0, stage 04, stage 05, t9, t10 on existing
caches. Status after block A: t0 clean, stage 04 clean, t9 ran without the
block jackknife, stage 05's ratios are meaningless where `xi` crosses zero
(`xi_FCT(600) = 0.00015`, so the 0.0488 and 1.6129 in the log are noise —
plot the difference, not the ratio), and t10 contradicts itself and fits above
the window (`k < 0.0653 s/km` against a DESI top of 0.0319 and a
`common_window` top of 0.0255 — rerun with a tighter `--fit-kmax-frac`).

**A V2 failure does not invalidate all of V3, and the distinction matters.**
A purely multiplicative error in `tau` is absorbed entirely by `A` and leaves
`r(k)` untouched, because both runs go through the same extractor. An error
that depends on density or temperature does move it, because that is exactly
where FCT and CDM differ. V2 tells us which case we are in.

**New plan (Pato, 2026-09-01): run the analysis on `murgia/{cdm,M2,M3}`.**
This is not a substitute for the FCT/CDM result — murgia is a different box
(20 Mpc/h), different cosmology (`h=0.6774`, `Omega_m=0.30749`) and different
models. It is something better: **an external, published check of the whole
chain.**

The reference is **Murgia, Scelfo, Viel & Raccanelli 2019, PRL,
arXiv:1903.10509, "Lyman-alpha forest constraints on Primordial Black Holes as
Dark Matter"**. Why it lines up:

- Same box and resolution: the paper's PBH grid is **512^3 particles in a
  20 Mpc/h box**. `murgia/cdm/los_0003.hdf5` reports box 29.52465 internal =
  20 Mpc/h, and relos.py walked 124.66M gas particles, i.e. 93% of 512^3,
  consistent with ~7% converted to stars by z=5.
- Same redshifts. Our murgia LOS files sit at z = 4.4, 4.6, 4.8, 5.0, 5.2,
  5.4, 5.6. The paper's data bins are **z = 4.2, 4.6, 5.0, 5.4** (MIKE and
  HIRES/KECK). Three of four are on our grid. That is not a coincidence.
- Same physics as our Poisson term. Their `P_PBH(k) = 1/n_PBH` (their Eq. 1)
  with `n_PBH = Omega_DM rho_cr f_PBH / M_PBH` (Eq. 2) is exactly what
  pending item 6 below has to compute. Their isocurvature piece is
  `P_iso = f_PBH^2 P_PBH = (2 pi^2 / k^3) A_iso (k/k_*)^(n_iso - 1)` with
  `k_* = 0.05 /Mpc` and `n_iso = 4`, and
  `P(k,z) = D^2(z) [T_ad^2 P_ad + T_iso^2 P_iso]` (Eq. 3). **Use these
  equations for stage 03 rather than re-deriving them.**
- They ran GADGET-III, we run SWIFT. So this is a genuine code-independent
  comparison, not a re-run.

**Guess to verify, not to assume: `M2` and `M3` are probably `M_PBH = 10^2`
and `10^3` M_sun with `f_PBH = 1`**, which are exactly the two models in the
paper's Figure 1. Check `ICs_parameters` in the murgia snapshots, or ask
whoever produced them. If that is right, Figure 1 (relative difference in the
1D flux power at z=5, in k [h/Mpc]) and Figure 2 (1D flux spectra at the four
redshifts for `M_PBH f_PBH` = 30, 60, 170, 1000 M_sun) are **published curves
we should be able to hit**. That is a stronger end-to-end validation than
Sherwood: same box, same resolution, same model class, answer in print.

Two things from that paper that bear directly on our own framing:

- They fit **the full shape of the 1D flux power, not a single amplitude
  parameter**, and they marginalise over the mean flux (9 grid points,
  `{0.6...1.4} x F_REF`) and over 8 further `tau` rescalings. So "do not fix
  `A`, fit it" is not our invention, it is what the reference analysis does.
  Cite it when writing the rescaling section.
- Their Figure 1 makes the opposite-facing point to ours and both are true:
  non-linear evolution washes out the PBH signature in the **3D matter**
  power, so the **1D flux** power is the far more effective probe. Do not
  write "P1D is incomplete" in a way that reads as "P1D is weak" — the
  honest claim is that P1D is a second moment and carries no information from
  saturated pixels, so tail-sensitive statistics are *complementary*, not
  that P1D should be replaced. This is a Viel paper; the framing will be read
  carefully.

---

**Pending, in order:**

1. **Re-run stage 02 on `regen/*_uni512_seed12345.hdf5`.** Cheap, and it
   decides whether the pairing was ever broken, whether t8 is valid, and
   which significance gets quoted. Do this first.
2. **V2 with Maria.** `data/los_murgia_cdm_z5_20.hdf5` is cut and ready. Ask
   for `tau` shaped `[n_los, n_pix]` keyed by `LOS_XXXX` group name, plus her
   `n_pix`, `dv`, `Gamma_HI` and velocity convention. Our grid: 2048 px,
   `v_box = 2738.925 km/s`, `dv = 1.337366 km/s`, axis 2,
   `Gamma_HI = 4.30e-13 s^-1` from TREECOOL_HM12_G+Q at z=5. **Ask for tau,
   not flux**: this file has `tau_eff = 2.85` raw, most pixels are saturated,
   and a 30% tau disagreement inside a trough is invisible in `F`. Our own
   answer on those lines is `cache/cache_murgia_cdm_z5_first100.npz`.
3. **V3 on murgia against Murgia+2019.** Confirm what M2 and M3 are, extract
   at z = 4.6, 5.0, 5.4, and compare the ratio to the paper's Figures 1-2.
4. **Fix and rerun the two broken pieces of V3**: stage 05 must not print a
   ratio through a zero crossing; t10 must fit inside `common_window`.
5. **Extract LOS caches at more redshifts for the FCT/CDM pair.** Queue job.
   Still the biggest blocker: the whole argument is redshift evolution and
   there is one redshift on disk. 17 LOS files exist, z=5 to z=1.8.
   `t11_bias_bk.py` refuses to interpret a single point and says so.
6. **Stage 03 — `A_P` and `A_b` from the initial conditions.** Does not exist
   yet and `t11` cannot run without it. Use Murgia+2019 Eqs. 1-5 above.
   Measure the z=198 box with Pylians and check the white-noise plateau sits
   where `A_P` predicts. Watch the `h` convention: `(Mpc/h)^3` vs `Mpc^3` is
   `h^3 = 0.316`, which is why `t11` has no default for `--units`.
7. **`k_F(z)`, the filtering scale.** Fit from the P_gas/P_matter suppression
   already measured (Gnedin & Hui). The most fragile input in the programme,
   because the broken term scales as `k_max^4`.
8. `t11` for real, once 6 and 7 exist.
9. **Transferability (GOAL stage C).** Repaint at three `tau_eff` and two
   `T0`, check `R(k)` moves at the percent level. A gate: if `R(k)` is not
   stable the ratio is not the observable. `--impose-trho` exists in
   `legacy/prepatch/`.
10. **t7.** Now unblocked by V1. It is the only bound available on how much of
    `b` is QLA, and `b(QLA)` / `b(stars back)` are the two ends of the
    systematic interval on the paper's central number.
11. The 40 and 80 Mpc/h boxes are being re-run correctly by Pato (the `h`
    mix-up truncates **every** box in that family, not just the 40 pair:
    `cdm-box-80` and `fct-box-80` also stop at 80 of 117.4743, 45.6%
    coverage). Fold in the 15 new runs and the `more_power` batch afterwards.
    Their ICs were made **without** monofonIC's `masked=2`, so they are the
    "hard" case for t7 — run `stages/00_inspect_snapshot.py --deep` on each
    and read the verdict rather than assuming.

**Parked deliberately, both good ideas, neither urgent:**

- **Profile over `A` instead of fixing it.** Grid over `(A_CDM, A_FCT)`, both
  free, and ask whether *any* choice reproduces `r(k) = 1`. It cannot: the
  tilt `A` can inject over the DESI window is a few per cent while `r(k)`
  swings 43%. Two quotable statements come out — the systematic (vary `A`
  within the observed `tau_eff` error bar) and the robustness (no `A` imitates
  the effect). This is the clean answer to Viel's objection, and Murgia+2019
  already does the marginalised version.
- **Dark gaps.** Distribution of contiguous runs with `F < threshold`. Runs on
  existing caches, is a real observable (Becker+2015, Zhu+2021), is directly
  sensitive to the saturated tail that P1D cannot see, and is much less
  sensitive to continuum error than the raw flux PDF. One figure, no new runs.

**Non-obvious facts, kept from earlier sessions:**
- **Clementina has no direct outbound internet.** All git traffic goes through
  the HTTP proxy `172.28.3.3:3128`. See section 6. DESI DR1 files have to be
  downloaded on the laptop and `scp`-ed.
- **`relosz.py` is the entry point; `relos.py` is the engine.** relosz
  resolves the LOS file and the snapshot from a run directory and a redshift,
  matching on the redshift each file *reports* rather than the index in its
  name, and knowing that snapshots may live in their own subdirectory. It
  looks for snapshots as DIRECTORIES, so on a run that stores them as plain
  files beside the LOS files (the 40 Mpc/h pair) it reports "No hay snapshot a
  z=..." and you fall back to relos.py with an explicit `--snapshot` glob.
- **Snapshot and LOS indices do not correspond, in either direction**, and the
  offset is not consistent between runs. Match on the redshift read out of
  each file. relos.py refuses a snapshot more than `dz = 0.02` away.
- **murgia is the only untruncated run on disk** until the 40/80 re-runs land.
  97.6-99.3% transverse coverage, ParticleIDs present, and `frac_below_floor
  = 0` with `wsum_raw_med = 0.923`, against 0.3145 and 0.755 for the truncated
  z=3 file. A third of the pixels in the truncated files hit the Shepard
  floor; in murgia none do.
- **The only (LOS, snapshot) pair in murgia within `dz <= 0.02` is z=5.0**
  (`los_0003` + `murgia-*-lyman_0002`). The example in the header of
  `scripts/sbatch_roundtrip.sh` says 5.6, which dies on arrival — `los_0000`
  is at z=5.6 but the nearest snapshot is z=5.0.
- **Step 13 of `run_validation_A.sh` is buggy.** `CLEAN=$(head -1
  logs/v1_candidates.txt | cut -f2)` takes the *first* candidate, not a clean
  one, so it ran on the truncated `2-fct-box-40-1024/los_0015.hdf5` under an
  echo claiming the file is untruncated. `data/los_clean_sample100.hdf5` is
  misnamed and should be deleted. Filter on the ray range, not on position.
- **`desi_window()` floats with redshift** (`k_max = 0.5 pi / R_z`). Use
  `common_window(zs)` for anything integrated across redshift bins.
- **`D^2(2.2)/D^2(4.0) = 2.40`**, from `units.growth_factor`. That is the
  calculable part of the redshift lever; everything above it is filtering
  scale, thermal state and QLA.
- **`b` measured on these runs is `b` under QLA.** The scheme converts 49.8%
  of the FCT baryons against 14.1% in CDM. That makes the central measurable
  a calibration of this subgrid scheme rather than a portable transfer
  function, and it is the first thing a referee will push on. Say it out loud
  in the paper.
- **Mass-separability of the converted gas is per-run.** The original 40 Mpc/h
  pair used monofonIC's `masked=2`, so DM and QLA-converted gas have different
  particle masses in `PartType1` — the easy case `--deep` detects. The
  `more_power` batch and the 15 new runs were made without the mask; treat
  them as the hard case.
- **Lead with results, not code.** Pato is working as a supervisor on this
  repository. Give him verdicts, figures, and what would falsify them. Batch
  commands so one console dump answers several questions. When a stage passes,
  say which files he now has to read and which stay machinery.
- **The objectives document lives outside this repository, deliberately.** Ask
  Pato where. Decide the `.gitignore` question before it is ever committed: a
  file that reaches a public repo stays in the history, and forks keep it.

---

## 1. What this is about, scientifically

Two cosmological simulations, run with SWIFT, from **identical initial
conditions**. Same box (40 Mpc/h, 1024³ gas particles), same background
cosmology, same subgrid physics. The only difference is the primordial
power spectrum:

- **CDM** — the standard spectrum.
- **FCT** — a broken spectrum (Sureda et al. 2021, `n_b = 2`,
  `k_t = 10 Mpc⁻¹`) plus a Poisson isocurvature term from the discreteness
  of primordial black holes treated as dark matter.

We shoot sightlines through both boxes at z = 3, compute the Lyman-α forest
optical depth `tau` along each, and compare the 1D flux power spectrum
`P1D`. The ratio `P1D(FCT) / P1D(CDM)` sits near **0.78 at large scales**
and rises through unity at small scales.

**The whole question of this repository is whether that 0.78 is physics or a
measurement artefact.** Pato presented this at COSMO-26 in Leiden on 27
August 2026. Matteo Viel's objection there was that the normalisation to a
common effective optical depth could be manufacturing the result. The
repository exists to answer that publicly and reproducibly.

### What we already know

Five explanations were tested to completion and all failed:

| candidate | test | verdict |
|---|---|---|
| rescaling to a common `tau_eff` | A = 1 in both runs | discarded — the effect is there before any rescaling. Rescaling makes it *shallower*: `r(0.003)` goes 0.7338 → 0.7797. It is not flat in k either, and injects a 2.35 pp tilt across the window; quote it as a systematic (section 0) |
| sampling variance | paired bootstrap, split-half | discarded — 20σ |
| SPH deposition scheme | no Shepard normalisation | discarded — survives |
| dense gas removed by QLA | cut at Δ > 100 and Δ > 30 | discarded — under 1% |
| thermal history | impose T–ρ with γ = 1.40 | discarded — about 2% |

Then we measured the 3D power spectra directly, which is what settled it:

| field | FCT/CDM ratio at k = 0.107 Mpc⁻¹ |
|---|---|
| matter | **1.000** |
| baryons (gas + particles converted to stars) | **1.000** |
| gas only | **0.808** |

The matter fields of the two runs are identical to one part in a thousand.
So are the baryon fields. **The gas field alone is down by 19%.**

The reason: the quick-Lyman-α (QLA) star-formation scheme converts **49.8%
of the baryons in the FCT box** against **14.1% in CDM**. The gas left
behind traces the same matter with a lower bias. Two controls pin down that
it is *where* the gas is removed and not *how much*:

- removing 40% of the gas **at random, everywhere** → power ratio 0.9999,
  i.e. costs nothing;
- removing gas **above a density threshold** until the remaining mass
  matches FCT's → ratio drops to about 0.10, i.e. costs everything.

FCT sits at 0.808, between those extremes. Its removal is concentrated
enough to matter and diffuse enough not to destroy the field.

### The conclusion, and why it is careful

The signal is real and robust, but it is **not a clean probe of the
primordial spectrum**. It acts through the bias of the surviving gas, and
its amplitude is set by a subgrid prescription with no feedback that
converts half the baryons while the observed stellar fraction at z = 3 is a
couple of percent. A run with feedback would return much of that gas to the
IGM and the deficit would shrink by an amount we cannot currently estimate.

**Do not let the framing drift in either direction.** The honest result is
neither "the model is ruled out" nor "the signal is an artefact to correct
away". It is: we identified the channel through which the forest constrains
these models, and pinning a number on it needs a run with feedback.

---

## 2. What is in this folder

```
Lyman-Alpha-PBHs/
  README.md          per-file documentation: every option, and why you'd use it
  HANDOFF.md         this file
  common/            cosmology, the P1D estimator, cache I/O, run registry, UVB
  stages/            the pipeline, numbered
  tests/             the falsification attempts, numbered
  legacy/            EMPTY until migration — see below
  paper/audit.tex    the document, LaTeX, drops straight into Overleaf
  scripts/           migrate, verify, sbatch templates
  cache/ data/ figures/ logs/    empty, gitignored
```

**The design boundary is the array of `tau`.** Everything upstream of it is
interchangeable (our extractor, or spectWizard). Everything downstream lives
in `common/p1d.py` and is used by every script. That is what lets us swap
the extractor and know the estimator did not change underneath.

`common/p1d.py` has unit tests in `tests/test_estimator.py` that check the
arithmetic against closed-form answers (a cosine must land in one bin with
amplitude `L a²/4`, `xi(0)` must equal the variance of the flux contrast,
`solve_A` must invert `tau_eff` to machine precision). **Run them first;
they need no data and take a second.**

---

## 3. Hard rules

1. **Never fabricate, simulate, or mock scientific data.** Every number in
   this project comes from a real snapshot. If a script needs data that is
   not there, it must fail loudly, not invent a plausible array. The only
   synthetic inputs allowed anywhere are the analytic ones in
   `tests/test_estimator.py`, and that file says so in its docstring.
2. **Do not edit anything in `legacy/`.** That is the code the published
   numbers came from. The pipeline imports it, never rewrites it. If a
   behaviour needs changing, wrap it from outside — `stages/01_extract_los.py`
   shows the two ways (override the metadata object, wrap `_cgs_factor`).
   If you genuinely must edit it, say so explicitly in the commit message.
3. **Do not commit `cache/*.npz` or `data/*.hdf5`.** They are hours of
   compute and gigabytes. `.gitignore` already excludes them.
4. **Do not commit the GitHub token, or any credential, ever.** Check
   `git diff --cached` before every commit.
5. Long jobs on the cluster go through `sbatch`, **never** `srun --pty`. An
   interactive session dies with the SSH connection and takes eight hours of
   extraction with it.

---

## 4. The cluster, and where the data lives

Two clusters, both at CCAD. **Never confuse them:**

- `snmgt01` = **Clementina** ← all of this work
- `serafin` = **Serafín**

Data root: `/data/contrib/pad_140/pcolazo/lyman/`

```
lyman/
  cdm-box-40-1024/
    cdm-40-m6-lyman_0003.hdf5      snapshot at z=3 (virtual file, 0 bytes)
    cdm-40-m6-lyman_0003.{0..3}.hdf5   the actual pieces, ~11 GB each
    los_0010.hdf5                  LOS file at z=3, 6144 sightlines, 3.5 GB
  2-fct-box-40-1024/
    fct-40-m6-lyman_0003.hdf5      same layout
    los_0010.hdf5                  6144 sightlines, 2.4 GB
  more_power/NB_{1..6}/  poisson_{1..6}/   the next generation of runs
```

Facts worth having in mind:

- Snapshots exist only at **z = 198, 7, 5, 3, 2**. LOS files exist at
  **17 redshifts from z = 5 down to z = 1.8**, in steps of 0.2. z = 3 is
  `los_0010` in **both** runs — same index, which is convenient but is not
  something to rely on for other runs.
- The snapshots use SWIFT's `distributed: 1` output: a **virtual** `.hdf5`
  of zero size alongside `.0.hdf5 … .3.hdf5` pieces. Read the virtual file.
  Globbing all of them and reading every one counts every particle twice.
  `common/runs.py::_virtual_first` handles this; do not undo it.
- Both snapshots contain `PartType0`, `PartType1` **and `PartType4`**.
  If `PartType4` is non-empty this is the easy path for test t7 — see the
  open questions below.
- The FCT LOS files are 2.2–2.4 GB against CDM's 3.4–3.8 GB, with the same
  number of sightlines and the same box. That 37% difference in raw bytes
  *is* the missing gas, visible with `ls`. It is a nice sanity check and
  worth a sentence in the paper.
- Each LOS file has **6144** sightlines. Earlier analysis used 1536, so
  either `max_los` was set or only a subset was processed. Worth resolving
  before quoting error bars.
- `lyman/murgia/{cdm,M2,M3}` is a separate batch (Murgia's simulations),
  not in the `common/runs.py` registry yet. `run_dir()` accepts a raw
  directory path, so `--run /data/contrib/pad_140/pcolazo/lyman/murgia/cdm`
  works without touching the code; add registry entries only if short
  names are wanted. It has only 3 snapshot outputs and 7 LOS files, a
  different grid from the main pair — run `stages/01_extract_los.py --run
  <path>` with no `--z` first to list what redshifts actually exist before
  assuming any file is z=3.

Environment on the cluster: `conda activate astro`.

---

## 5. What to do, in order

### Step 0 — sanity, locally, no data needed

```bash
cd Lyman-Alpha-PBHs
python -m compileall -q common stages tests && echo compile ok
python tests/test_estimator.py
```

### Step 1 — GitHub

See section 6 for the token. Then:

```bash
git init -b main
git add -A
git commit -m "reproducible Lyman-alpha pipeline: FCT vs CDM"
git remote add origin https://github.com/patricio-c/Lyman-Alpha-PBHs.git
git push -u origin main
```

Public repository — it is meant to be linked from the paper.

### Step 2 — on Clementina: clone and migrate

Clementina has no direct outbound internet; git has to go through the
cluster's SSH proxy. Set up `~/.ssh/config` and the GitHub SSH key first
(section 6), then:

```bash
ssh snmgt01
cd /data/contrib/pad_140/pcolazo
git clone git@github.com:patricio-c/Lyman-Alpha-PBHs.git
cd Lyman-Alpha-PBHs && conda activate astro
bash scripts/migrate.sh /data/contrib/pad_140/pcolazo/LOS
bash scripts/verify.sh
```

(This was already done: the repo is cloned on Clementina, not yet
migrated.)

`migrate.sh` **copies and deletes nothing.** It brings his working scripts
into `legacy/`, the pre-patch `.bak_delta` / `.bak_trho` files into
`legacy/prepatch/` (the only record of the code before the `--delta-max`
and `--impose-trho` flags existed), the `cache_*.npz` files, and the
TREECOOL tables.

Pato wants to delete the old `LOS/` directory afterwards. **Only after**
`verify.sh` prints `VERIFY OK` and the migration is committed and pushed.
Do not suggest deleting it before that.

Then commit `legacy/` and push.

### Step 3 — the two checks that gate everything else

```bash
python stages/02_check_los_match.py --run-a cdm40 --run-b fct40 --z 3.0
python stages/01_extract_los.py --run fct40 --z 3.0 --geometry-only
```

The first one matters a lot and has never been run. The sightline positions
come from SWIFT's LOS output configuration, not from anything chosen at
analysis time, so "line *i* of CDM is line *i* of FCT" is an **assumption**.
Two results depend on it: `t8_single_los.py`, which compares line *i*
against line *i* pixel by pixel, and the paired bootstrap that produced the
20σ. The ensemble results never assumed a pairing and are unaffected either
way. If the lines are unmatched, the script looks for a permutation so the
caches can be reordered rather than re-extracted.

The second prints `v_box`, `dv`, `k_fund`, `k_nyq`, the unit conversions
and the DESI window without reading a single particle. Expected at z = 3
for this box: `H = 306.3 km/s/Mpc`, `v_box = 4498 km/s`,
`dv = 2.1963 km/s` at 2048 pixels, `k[s/km] → k[h/Mpc]` is `× 112.45`.
**If those numbers come out different, stop and work out why** before
running anything else.

### Step 4 — the analysis, on the caches that already exist

Minutes each, no queue:

```bash
python tests/t0_rescaling.py --cdm cache/cache_cdm.npz --fct cache/cache_fct.npz --out figures/t0_rescaling
python stages/04_p1d.py cache/cache_cdm.npz cache/cache_fct.npz --labels CDM FCT --out figures/p1d
python stages/05_xi.py  cache/cache_cdm.npz cache/cache_fct.npz --labels CDM FCT --out figures/xi
python tests/t8_single_los.py --cdm cache/cache_cdm.npz --fct cache/cache_fct.npz --pick median --out figures/t8
```

`t0_rescaling.py` is Matteo's test and it is the first section of the paper.

### Step 5 — t7, the test that closes the argument

Check first what `PartType4` actually holds:

```bash
python stages/00_inspect_snapshot.py --run fct40 --z 3.0 --deep
```

It ends with a verdict block naming the exact flag for t7. Then t7 writes an
augmented snapshot with the converted gas returned to `PartType0`. **One
step in between:** `swift_extract.py` reads SWIFT LOS files, not snapshots,
so sightlines have to be shot through the augmented snapshot before stage 01
can touch it — either by re-running SWIFT in LOS mode, or with
`legacy/relos.py` / `legacy/relosz.py`, which is how the current LOS files
were regenerated after an earlier bug.

### Step 6 — the document

`paper/audit.tex` compiles with pdflatex, no exotic packages. Pending
numbers are marked `\pending{...}` and render in red. Fill them from the
`.txt` files each script writes next to its figure.

Style for anything written into the paper: direct, concrete, lab-notebook
English. Imitate Nelson Padilla's or Federico Stasyszyn's register. **No
AI-sounding phrasing** — no "it's worth noting", no "not X but Y", no
em-dash asides, no hedging strings. Short declarative sentences. When a
result is weak, say it is weak.

---

## 6. GitHub — exactly which permissions

Two decisions first, then the table.

**Do not paste the token into a chat.** Put it in the credential helper or
the MCP config file and refer to it by name.

### Recommended path: create the empty repo in the browser first

Creating a repository through the API needs `Administration: Read and write`
on a fine-grained token, which cannot be scoped to a repository that does
not exist yet — so it ends up scoped to *all* repositories, which is much
broader than this job needs. Creating the empty repo by hand takes thirty
seconds and removes that entirely.

1. github.com → **New repository**
2. Name it, set **Public**, and **do not** initialise with a README,
   `.gitignore` or licence. An initialised repo makes the first push
   conflict and you have to `--force` or rebase.
3. Copy the HTTPS URL.

### Fine-grained personal access token

github.com → Settings → Developer settings → **Personal access tokens** →
**Fine-grained tokens** → Generate new token.

- **Token name:** something like `claude-lya-repro`
- **Resource owner:** your personal account (unless the repo lives under an
  organisation — see the open questions)
- **Expiration:** 90 days is a sensible default. It is renewable.
- **Repository access:** **Only select repositories** → pick the one you
  just created. Not "All repositories".

**Repository permissions** — set only these:

| permission | level | why |
|---|---|---|
| **Metadata** | Read-only | mandatory, GitHub enables it automatically and will not let you turn it off |
| **Contents** | **Read and write** | this is the one that matters. Push commits, create branches, read files |
| Pull requests | Read and write | *only* if you want Claude to open PRs instead of pushing straight to `main` |
| Issues | Read and write | *only* if you want Claude to file issues as a to-do list |
| Workflows | Read and write | *only* if you ever commit a file under `.github/workflows/`. There is no CI in this repo, so **leave it off**. If you add CI later, a push touching that directory will fail with a confusing error until you enable it |

**Leave everything else off.** In particular: Administration (not needed
once the repo exists), Actions, Secrets, Environments, Packages, Pages,
Webhooks, Deployments, Codespaces, Dependabot, Security advisories, Code
scanning, Secret scanning.

**Account permissions:** none. Leave every one at "No access".

Minimum viable set is exactly two lines: **Metadata: read** and
**Contents: read and write**.

### Classic token, if you use `gh repo create` instead

Classic tokens are coarser. If you want the CLI to create the repo for you:

| scope | tick it? |
|---|---|
| `public_repo` | **yes** — enough for a public repository |
| `repo` | only if the repo will be **private**; it is the full-access scope |
| `workflow` | no, unless you add CI later |
| `delete_repo` | **no.** Never. There is no reason to grant it |
| `admin:org`, `gist`, `user`, `notifications`, everything else | no |

### On Clementina specifically: git over SSH through the proxy

The cluster has no direct outbound internet — every outbound SSH connection
goes through an HTTP proxy at `172.28.3.3:3128`. Add to `~/.ssh/config`
(create the file if it does not exist; `chmod 700 ~/.ssh`, `chmod 600` the
key):

```
Host github.com
   User git
   IdentityFile ~/.ssh/id_ecdsa
   ProxyCommand ncat --proxy 172.28.3.3:3128 --proxy-type http %h %p
```

Then generate a key on Clementina if you do not already have one
(`ssh-keygen -t ecdsa`), copy `~/.ssh/id_ecdsa.pub`, and add it at
github.com -> Settings -> SSH and GPG keys -> New SSH key. Test before
cloning: `ssh -T git@github.com` should answer "Hi patricio-c! You've
successfully authenticated...". `Permission denied (publickey)` at this
step means the public key was never added on the GitHub side, not a proxy
problem. Clone with the SSH form, not HTTPS:

```bash
git clone git@github.com:patricio-c/Lyman-Alpha-PBHs.git
```

This is untested for the HTTPS+token flow above — the docs at
`docs.clementinaxxi.org.ar` only confirm the SSH route through the proxy,
so use SSH on the cluster and HTTPS+token on the laptop.

### Storing it

```bash
git config --global credential.helper store   # or 'osxkeychain', or 'libsecret'
git push -u origin main
# username: your GitHub username
# password: paste the token, NOT your GitHub password
```

Or with the CLI:

```bash
gh auth login --with-token < /path/to/tokenfile
```

Delete the token file afterwards.

### If you wire up the GitHub MCP server instead

The same fine-grained token works. It goes in the Claude Desktop config as
an environment variable for the server (usually
`GITHUB_PERSONAL_ACCESS_TOKEN`). The config file lives at
`~/.config/Claude/claude_desktop_config.json` on Linux. Restart the app
after editing it. **That file is not inside this repository and must never
be copied into it.**

Plain `git push` from a shell is simpler and needs no MCP server at all. Use
the MCP route only if you want Claude to read and edit files on GitHub
without a local clone.

---

## 7. Open questions — ask, do not guess

1. **Does `PartType4` actually contain particles at z = 3?** The group
   exists in both runs, but it also exists at z = 198 where nothing can have
   formed, so SWIFT may be writing empty groups. Run
   `stages/00_inspect_snapshot.py --run fct40 --z 3.0 --deep` and read the
   verdict block. This decides which of three routes t7 takes, and it is the
   single most valuable unknown right now.
2. **6144 or 1536 sightlines?** **Closed 2026-09-01, and it was a false
   dichotomy.** They are different files, not a subsample. SWIFT's own LOS
   files hold 6144; the analysis used `regen/*_uni512_seed12345.hdf5`, a
   `relos.py --uniform 512` regeneration = 512 rays per axis x 3 axes = 1536,
   drawn uniformly over the whole box face. Read the `source` field of the
   cache. See section 0.
3. **Are the sightlines matched between runs?** **REOPENED 2026-09-01. The
   2026-08-31 "no" was measured on the wrong files.** Stage 02 compared the
   two 6144-ray SWIFT files, which are indeed unrelated; the caches come from
   the 1536-ray uniform regenerations, which share a seed and a box and are
   therefore expected to be identical in both runs. t9's measured per-line
   correlation of +0.86 between the two runs is consistent with matched lines
   and impossible for unmatched ones. Re-run stage 02 on
   `regen/cdm40_z3.0_uni512_seed12345.hdf5` and
   `regen/fct40_z3.0_uni512_seed12345.hdf5`. Until that is done, treat both
   the "pairing is broken" claim and the corrected significance as unsettled,
   and do not send anyone a correction to the 20σ.
4. **What temperature to give the reinjected gas in t7?** There is no right
   answer, so the script sweeps `T0` by factors of 0.5, 1 and 2 and the
   result has to be shown insensitive to the choice. Do not quietly pick one.
5. **Does `swift_extract.py` agree with SpectWizard?** Two independently
   written codes computing the same physics (SPH deposition + Voigt
   integral along the sightline) are not guaranteed to agree in detail —
   `swift_extract.py`'s own bug history (a 5.6x tau overestimate from a
   missing impact parameter, catastrophic cancellation in the Voigt
   approximation, see `paper/audit.tex` appendix) shows how sensitive this
   calculation is to implementation choices, even within one codebase. Not
   yet tested. The design boundary (section 2) makes this cheap: get
   María to run SpectWizard on the same LOS file, save a `tau`/`dv`/`z`
   cache in the same format, and feed it into `stages/04_p1d.py` alongside
   the `swift_extract` cache — same pattern as the existing `--sherwood`
   flag. Worth doing as an independent, external check, not just the five
   internal ablations already in the "what we already know" table.

---

## 8. If something does not work

- `could not import legacy/swift_extract.py` → `scripts/migrate.sh` has not
  been run, or was pointed at the wrong directory.
- `no SWIFT LOS files found` → the error message lists the hdf5 files it did
  find. The run directory is probably wrong.
- `closest LOS file is z=…, you asked for z=…` → deliberate. Nothing in this
  repository uses snapshot indices, because they are not comparable between
  runs. Pass `--los-file` explicitly if you really mean that file.
- A geometry number disagrees with section 5 → **stop.** That is exactly the
  class of bug this repository exists to catch. A units mismatch between Mpc
  and Mpc/h once silently truncated 31.9% of every sightline and cost weeks.

---

## 9. Two clones, one remote — the sync flow

There are two working copies: the laptop (where code gets written and
pushed via the GitHub MCP, no local git needed there) and Clementina
(`/data/contrib/pad_140/pcolazo/Lyman-Alpha-PBHs`, where the actual extraction and
analysis run, cloned over SSH as in section 6). Both point at the same
`https://github.com/patricio-c/Lyman-Alpha-PBHs` remote. Ordinary
distributed-git discipline applies, nothing special:

- Whoever is about to edit code, `git pull` first (or, on the laptop side,
  just re-read the file before editing — the MCP push always goes straight
  to `main`).
- Small, logical commits beat one giant one, same as any repo.
- `cache/*.npz` and `data/*.hdf5` are gitignored on purpose — caches made on
  Clementina do not travel to the laptop and do not need to. Only code,
  docs and `paper/` cross the wire.
- If both sides edit the same file in the same window, it is a normal merge
  conflict, resolved normally. It has not happened yet because the split
  so far is clean: code from the laptop, `legacy/` and caches from
  Clementina.

## 10. Working across Claude sessions

This file is written so that a new Claude session can pick up the project
cold: read section 0 for what is done and what is next, and the rest for
context as needed. The convention going forward:

- **At the end of a session that changed something real** (migrated code,
  ran a gating check, resolved an open question), update section 0 before
  closing: move finished items out of "Pending," add anything genuinely
  new to "Non-obvious facts," and update the date. Keep the rest of the
  file mostly stable — it is background, not a log, and does not need to
  be rewritten each time.
- **A new session** starts by reading this file (from the GitHub repo or
  the local clone, either is fine, but push/pull first if unsure which is
  current) and treats section 0 as the actual task list. No separate
  scratch/TODO file is needed for this repo specifically — the overhead of
  keeping a second file in sync with this one is worse than the cost of
  section 0 being slightly denser than a pure TODO would be. If a task
  genuinely does not belong in a project handoff (something ephemeral,
  laptop-only, unrelated to the science) it does not need to go here at
  all.
- Do not let section 0 silently drift out of sync with what actually
  happened — if you are not sure whether the previous session finished
  something, check the repo state (`git log`, does the file exist, does
  the test pass) rather than trusting a stale bullet point.

**On model choice: switch to a stronger model for anything delicate.**
Routine work — git, running scripts, migration, formatting, following a
checklist like this one — is fine on the default model. Switch to a
stronger one (e.g. Opus) for moments where the reasoning itself is
load-bearing for the science: interpreting whether a result like the
LOS-mismatch above changes what the paper can claim, auditing `legacy/`
code for correctness rather than just running it, deciding how to frame a
borderline or weakened result for a critical reviewer (Viel), or writing
the parts of the paper that make the actual scientific argument. If it is
the kind of thing that would be embarrassing to get subtly wrong in front
of a reviewer, it is worth the switch.

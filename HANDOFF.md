# HANDOFF — read this first

You are picking up a scientific analysis repository. This file is the whole
context. Read it before touching anything.

The person you are working with is **Pato Colazo**, PhD student in astronomy
at IATE/OAC (Universidad Nacional de Córdoba, CONICET). He works in
Rioplatense Spanish and writes papers in English. He is on a Linux laptop;
the cluster is remote.

---

## 0. Current status — read this first, every session

Updated 2026-09-03. This section is the fast-changing punch list; the rest of
the file is stable background. Start every new session here, and dip into the
numbered sections below only for the *why* behind something.

**The headline, if you read nothing else.** We reproduce Murgia+2019 with an
independent code. Their Figure 1 was recovered exactly from the arXiv vector
source (not digitised by eye), and with the normalisation fixed our 1D flux
ratio sits within ~15% of their published curve across the whole fitted
window, while our non-linear matter power matches theirs to ~0.2 percentage
points. SWIFT + monofonIC + CAMB against GADGET-III + 2LPTic + CLASS. That is
now the most valuable thing in this repository, and it is what makes anything
else we claim credible.

**And the FCT/CDM pair is settled.** The z=3 deficit that this repository was
built to audit survives both normalisations: `0.7338` with `A = 1` everywhere
and `0.7857` rescaled to the CDM run's own opacity. It is not a normalisation
artefact. The headline 0.78 stands.

Two things died in the process and one survived. **The low-k deficit is dead**
— it was `units.tau_eff_turner24(5.0)` extrapolating outside its calibration
range. **The eyeballed Figure 1 column is dead** — it was wrong by more than a
factor two at high k. **The QLA conversion differential survives** — it is a
particle count, it owes nothing to any normalisation, and it now has a
sharper role than the one the previous version of this section gave it.

---

### Done this session (2026-09-03)

- **Murgia+2019 Figure 1 recovered exactly.** The arXiv e-print
  (`https://arxiv.org/e-print/1903.10509`) ships `relative_spectra.pdf` as a
  vector PDF. `pdftocairo -svg` exposes the polyline vertices matplotlib
  wrote; the calibration comes from the tick-mark paths in the same file. All
  six curves (linear / non-linear matter / 1D flux, for 10^2 and 10^3 M_sun).
  Two internal checks passed: their vertices land on exact multiples of
  `k_fund = 2 pi / 20 = 0.3142 h/Mpc`, so these are their P1D bins and not a
  smoothed trace; and their **linear** M3/M2 ratio is `10.00` at every k,
  which is exactly `P_iso ∝ M_PBH f_PBH`. Reproduce with
  `scripts/extract_murgia_fig1.py` and `notebooks/murgia_vs_swift.py`.
- **Pending 1 run** (`--norm none`) and then **the normalisation settled**
  (`--norm taueff --tau-eff 2.84299`). Numbers below.
- **The initial conditions validated** against their published curves, using
  the on-disk `p_matter_{cdm,m2,m3}z5.txt` from Pylians.
- **The FCT/CDM pair re-normalised at z=3.** Section below. The deficit is
  real, and the k-dependent tilt previously attributed to the rescaling was
  itself an artefact of the Turner target.
- **`notebooks/murgia_vs_swift.py`** reproduces the extraction, the tables
  and the figure end to end. It is in percent-cell format, so VS Code and
  Jupyter open it as a notebook while git still sees plain text. It reads
  `figures/p1d_murgia_z5_*.txt` when they exist and says so out loud when it
  falls back to transcribed console values.

---

### THE CROSS-VALIDATION — the new central result

**The normalisation is a real choice and it was being made wrong.** Three
normalisations of the same three caches, z = 5.0, 1536 sightlines. Excess over
`cdm` in per cent. The Murgia column is the recovered curve, interpolated in
log k to our anchors — quotable, not eyeballed.

    M3 = 10^3 M_sun
    k [s/km]  k [h/Mpc]   none   Turner   ref CDM   Murgia   refCDM/Murgia
      0.0030      0.411   +9.91    -3.47     -1.22    -2.55        -
      0.0050      0.685  +10.61    -5.73     -1.18    +0.62        -
      0.0100      1.369  +21.96    +4.97     +8.51    +7.19      1.18
      0.0200      2.739  +32.24    +9.69    +16.18   +18.55      0.87
      0.0300      4.108  +53.84   +18.77    +34.36   +33.84      1.02
      0.0600      8.217 +120.10   +43.25    +89.26   +78.22      1.14

    M2 = 10^2 M_sun
      0.0030      0.411   +4.80    +1.89     +1.25    -0.15        -
      0.0050      0.685   +5.95    +2.39     +2.24    +0.46        -
      0.0100      1.369   +6.85    +1.35     +2.82    +1.44      1.96
      0.0200      2.739   +8.57    +2.96     +4.23    +3.18      1.33
      0.0300      4.108  +11.17    +3.25     +6.64    +5.78      1.15
      0.0600      8.217  +21.74   +10.90    +16.70   +14.18      1.18

`none` overshoots by x1.54 (M3) and x1.53 (M2) at k = 0.06. `Turner`
undershoots by x0.55 and x0.77. **`ref CDM` lands at x1.14 and x1.18**, and
both models go to ~0 at large scales, as theirs do. A factor 2.8 spread that
depended on which flag you passed is now a ~15% residual.

**What `ref CDM` means, exactly.** `--norm taueff --tau-eff 2.84299`: all
three runs are rescaled to the *same* target, and that target is the `cdm`
run's own raw `tau_eff`. So `A(cdm) = 1.00000` by construction, `A(M2) =
0.96644`, `A(M3) = 0.89491`. No external observational target is imported.
This is the normalisation to use from here on, at every redshift, with the
reference run's own raw value as the target. It is defensible on its own
terms: the ratio is between models, so the only mean-flux difference that
needs removing is the one *between* them.

**Why the old default broke.** `units.tau_eff_turner24(5.0)` implied
`A(cdm) = 0.37877` — a factor 2.6 on the optical depth, which is not a
calibration, it is a rebuild of the forest. Everything else in this repository
works at z = 2.2 to 4.0. Treat that function as valid inside its calibration
range and nowhere else. It is the single confirmed bug found this session and
it is in `common/units.py`.

**The initial conditions are validated.** Non-linear matter power at z = 5,
excess over `cdm` in per cent, ours from Pylians on the snapshots against
their published dotted curves. `k[Mpc^-1] -> k[h/Mpc]` divides by
`h = 20/29.52465 = 0.6774`; checked that nothing below depends on that at
`h = 0.6774 / 0.681 / 0.702`.

    k [h/Mpc]   M2 ours   M2 Murgia    M3 ours   M3 Murgia
        4.11      -0.02       +0.07      +0.22       +0.33
        5.00      -0.02       +0.08      +0.38       +0.51
        8.22      +0.01       +0.14      +1.39       +1.58
       10.00      +0.04       +0.22      +2.09       +2.09
       15.00      +0.21       +0.39      +4.33       +4.08
       19.00      +0.43       +0.65      +6.65       +6.36

M3 agrees within ~0.2 points everywhere and is exact at k = 10. **The ICs and
the non-linear evolution are not the problem and are not a free parameter any
more.** M2 sits ~0.15 points low, but their M2 curve is pressed against the
zero of their axis and our extraction precision there is ~0.05 points — read
nothing into it.

Note in passing what panel b says on its own: +6.6% in matter at k = 19
against +89% in flux at k = 8. Their point that non-linear evolution washes
out the PBH signature in 3D matter is correct, and we reproduce it.

**What the residual ~15% at high k can still be.** It is the *same* factor for
M2 and M3, which means one common cause downstream of the density field. Ruled
out: the ICs and the non-linear evolution (panel b); the resolution (their PBH
grid is 512^3 in a 20 Mpc/h box, stated in their methods — the same as ours);
anything primordial at large scales (`n_iso = 4` makes `P_iso` white noise,
contributing as k^3, and their own linear curve gives +0.08% at k = 1).
Still open, in the order worth testing:

  1. **SPH sampling noise in our extractor.** `wsum_raw_med` degrades with the
     model: cdm 0.924, M2 0.888, M3 0.771 at z=5. Fewer gas particles per
     sightline means a noisier SPH density. Shot noise is white, so it cannot
     touch low k and **does** inflate high k — which is where the residual is,
     in the model that lost the most particles. This is pending 1 and it is
     the leading suspect.
  2. **Thermal state and UV background.** Their reference is
     `T_0(z=4.2) = 9200 K`, `gamma = 1.47`, Planck-2015 parameters
     (H0 = 70.2, Omega_m = 0.301, sigma_8 = 0.829, n_s = 0.961),
     `z_reio = 9`, mean flux anchored on BOSS. Ours is HM12 TREECOOL at
     `Gamma_HI = 4.30e-13 s^-1`, `h0 = 68.1`, `Omega_m = 0.3053`. Measure T_0
     and gamma in our z=5 snapshot and re-extract with `--impose-trho` at
     theirs.
  3. **The codes themselves.** GADGET-III SPH against SWIFT SPH, and their
     `T < 10^5 K` cut against our overdensity-only QLA plus the SWIFT entropy
     floor. This is the hardest to isolate and the least likely to be worth
     15% on its own.

---

### THE FCT/CDM PAIR AT z = 3 — the deficit is real

Run 2026-09-03 on `cache/cache_{cdm,fct}.npz`, 1536 sightlines, z = 3.0.

    run   tau_eff raw    A (target = CDM's own 0.42461)
    CDM      0.42461        1.00001
    FCT      0.40350        1.09058

    k [s/km]   norm=none   ref CDM   ref/none
      0.0030      0.7338    0.7857     1.0707
      0.0050      0.7939    0.8498     1.0704
      0.0100      0.8749    0.9347     1.0683
      0.0200      0.9945    1.0620     1.0679
      0.0300      1.0731    1.1479     1.0697
      0.0600      1.1144    1.1997     1.0765

**1. The deficit is not a normalisation artefact.** It is there with `A = 1`
everywhere (`0.7338`) and it is there rescaled to the reference run's own
opacity (`0.7857`). Quote the pair, not one of them: *the large-scale ratio is
0.73 unnormalised and 0.79 normalised*. Both are far from unity and the
conclusion does not depend on the choice. This is the opposite of what happened
to the murgia low-k deficit, and the contrast is the point — the same check
killed one result and confirmed the other.

**2. The k-dependent tilt was also a target artefact.** The previous version of
this file recorded that the rescaling is "not flat in k", moving `r(k)` by
+6.26% at k=0.003 against +3.91% at k=0.03, a 2.35-point tilt. That was
measured against the Turner target. Against the reference run's own value the
gain is **flat**: 7.07%, 7.04%, 6.83%, 6.79%, 6.97%, 7.65% across the whole
window — a spread of 0.9 points with no monotonic trend. A normalisation that
is a scale-independent gain is what one wants; the tilt was the mismatch, not
the method. Correct this wherever it is quoted, `paper/audit.tex` included.

**3. `A(FCT) = 1.09058` is a 9% correction — modest and defensible**, unlike
the factor 2.6 the Turner target demanded at z=5. And it points the other way
than murgia: FCT is *less* opaque than CDM (`tau_eff` 0.40350 vs 0.42461),
where M3 was *more* opaque than cdm. That is consistent with QLA converting
49.8% of the FCT baryons against 14.1% in CDM — far less gas, so less
absorption. The direction of the mean-flux shift is set by how much gas the
threshold removes, which is the QLA lever again, now visible in two pairs
pointing opposite ways.

**Still to do on this pair:** the same normalisation at every other redshift
once the caches exist, and error bars. Nothing at k < 0.007 s/km is quotable
until the jackknife runs.

---

### THE QLA THRESHOLD RESULT — what survives, and it is sharper now

**The measurement, exact, from the z=5 snapshot headers.** Initial gas is
`512^3 = 134217728` particles:

    run   gas at z=5     converted    (M_PBH, f_PBH = 1)
    cdm   124664686       7.12%       LCDM
    M2    123838499       7.73%       10^2 M_sun
    M3    113403060      15.51%       10^3 M_sun

`M3` converts **more than twice** what `cdm` does — an 8.4 percentage-point
gap. All three share the Panphasia descriptor
`[Panph6,L20,(235287,445214,422255),S1,KK1025,CH-999,COLIBRE050]` and the same
`Omega_cdm`: same realisation, different input transfer function, so sample
variance cancels in the ratio at the level of the ICs. Counted near the 1536
sightlines instead of globally the gap is 14.2% rather than 9.0%, which is the
bias change showing up directly.

**`PartType4` is empty in all three runs**, yet gas is missing. SWIFT's QLA
sends converted gas to `PartType1` as collisionless particles, exactly as
Viel+2004 describes. Do not count stars to get the conversion fraction; count
the missing gas, or split `PartType1` by mass. This answers open question 1
for murgia.

**The claim, restated so it survives the corrections.** The previous version
of this section said the threshold produces a low-k deficit in P1D that their
model cannot make. That is now known to be false — the deficit was the
normalisation. The claim that survives is narrower and aimed better:

> The star-formation threshold `Delta > 1000` is applied identically to both
> models, but it removes more than twice as much gas in the PBH run. That
> changes the run's **mean flux** — raw `tau_eff` is 2.843 in `cdm` against
> 3.008 in M3, a 5.8% shift. The mean flux is not a spectator: it is the
> nuisance parameter every Lyman-alpha constraint marginalises over, with a
> Gaussian prior of `sigma = 0.04` on `F_bar` in Murgia+2019. So a numerical
> choice in the subgrid model produces a model-dependent shift in a parameter
> that is degenerate with the one being constrained.

Viel+2004 validated QLA against a full multi-phase model **inside LCDM only**
— it says nothing about `QLA(model A) - QLA(model B)`, and the whole mechanism
of these models is to send more gas over the threshold. In their check both
sides had the same amount of gas crossing. In a PBH constraint they do not.

**This is a mechanism with a measured size, not a demonstrated bias.** We have
not shown that it moves `f_PBH M_PBH`. Saying so would be overclaiming and a
referee would catch it in one line. What closes the gap is pending 2: move the
threshold for real and watch whether `R(k)` moves. State it that way in the
paper.

**Do not repeat this mistake.** The QLA differential has now been credited
with a low-k deficit (wrong, it was normalisation) and a factor-two high-k
deficit (wrong, also normalisation). Each time the mechanism was plausible and
the evidence was an artefact. Before attributing anything else to it, check
what the same plot looks like under `--norm none` and under the cdm-referenced
target.

---

### CODE STATUS — what is safe to read and trust

**Validated, read these first, they are the argument:**

- `common/p1d.py` — the estimator, and the design boundary of the whole
  repository. Checked against closed-form answers in `tests/test_estimator.py`
  and now, indirectly, by reproducing a published P1D ratio. This is the file
  to read to understand what the project measures.
- `stages/04_p1d.py` — thin driver over `p1d.py`. Its `--norm` switch is the
  most consequential flag in the repository; read the three modes and the
  `p1d_from_tau(..., target=)` call.
- `stages/01_extract_los.py` — the extraction driver. Read the docstring for
  the convention switches (`--hubble-mode`, `--vel-convention`); the geometry
  block it prints is the first thing to check on any new run.
- `legacy/relos.py` and `legacy/relosz.py` — **V1 passed**: they reproduce
  SWIFT's own LOS output particle for particle and the flux field to numerical
  precision. Machinery; read only when regenerating sightlines.
- The IC chain (monofonIC + CAMB, FCT transfer function, Poisson term) —
  validated at z=5 against Murgia's published linear and non-linear curves.

**Trusted but not independently verified:**

- `legacy/swift_extract.py` — the `tau` engine, and the one with a bug history
  (a 5.6x overestimate from a missing impact parameter, catastrophic
  cancellation in the Voigt approximation; see `paper/audit.tex`). It now has
  strong *external* corroboration — reproducing a published P1D at the 15%
  level is a much harder test than any internal ablation — but V2 against
  SpectWizard is still the only true code-vs-code check and is still open.

**Known broken, do not read them for understanding until fixed:**

- `common/units.py::tau_eff_turner24` — wrong outside its calibration range.
  This is the confirmed bug of this session.
- `stages/05_xi.py` — prints ratios through a zero crossing
  (`xi_FCT(600) = 0.00015`), so the 0.0488 and 1.6129 in the log are noise.
  Plot the difference, not the ratio.
- `tests/t10_delta_p1d.py` — contradicts itself and fits above the window
  (`k < 0.0653 s/km` against a `common_window` top of 0.0255).
- `tests/t8_single_los.py` — its docstring still states the sightline pairing
  as fact. Fix the docstring whatever else happens.

---

### V1 — do we control the sightlines? **PASSED.**

`sbatch scripts/sbatch_roundtrip.sh .../lyman/murgia/cdm 5.0 100`, job 1524260,
`figures/t12_roundtrip_cdm_z5.0_n100.txt`. Ray positions identical to
`0.000e+00`; particle sets identical on 97 of 100 sightlines; `max |dF| =
5.96e-06`; `tau_eff = 2.851045` both; P1D ratio 1.0000 at every k. t12 prints
FAIL because its criterion is binary (`exact == n`) — **read the numbers, not
the verdict**. The three extra particles sit in pixels with `tau ~ 4000`.

The sentence this licenses: *our regeneration of SWIFT's line-of-sight output
reproduces the flux field and the 1D flux power spectrum to numerical
precision.*

Scope: one run, one redshift, 100 sightlines, all along axis 2, reusing the
template's ray positions. Not tested: the `--uniform` path, `--skip-los`
batching, axis detection on axes 0 and 1.

---

### The `uni512_seed12345` finding

**The production caches were NOT made from SWIFT's own LOS files.** From the
`source` field: `cache/cache_{cdm,fct}.npz` are `(1536, 2048)` from
`regen/{cdm,fct}40_z3.0_uni512_seed12345.hdf5`. `relos.py --uniform 512` draws
512 rays per axis over three axes = 1536, uniformly over the whole box face.
Verified bit-for-bit and confirmed by Pato.

1. "6144 or 1536 sightlines?" was a false dichotomy — different files, not a
   subsample.
2. The published numbers do **not** rest on truncated sightlines.
3. `tau_eff` 0.42461 (cache) vs 0.26192 (100-line cut) is not a bug: the cut
   came from the truncated SWIFT file.

**Stage 02 measured the wrong files** — it was pointed at the 6144 SWIFT
files, which are genuinely unrelated to each other. The caches come from the
1536 uniform regenerations, which share seed, box and `default_rng`, so they
are the same positions in both runs. That explains t9's
`sigma(common)/sigma(unpaired) = 0.3733` and per-line correlation `+0.8592`.
Re-run stage 02 on the `regen/*_uni512_seed12345.hdf5` pair before requoting
any significance.

---

### Untapped, and worth more than another extraction

Every murgia run directory holds files nobody has opened:

- `power_spectra/` — SWIFT computed **3D power spectra on the fly**. If
  `P_matter` and `P_gas` are there for all three runs, the measurement that
  settled the FCT argument (matter 1.000, baryons 1.000, gas 0.808) is
  available for free, without Pylians and without a queue. Do this first.
- `statistics.txt` — mass per particle type against time, i.e. the full
  conversion history rather than one number at z=5.
- `SFR.txt` — the same as a rate.
- `used_parameters.yml` — what SWIFT **actually applied**, defaults included.
  Cite this, not `M3.yml`. `unused_parameters.yml` shows what it ignored.
  This is where to confirm whether we apply a `T < 10^5 K` cut and what the
  entropy floor (`QLAEntropyFloor:over_density_threshold = 10`,
  `temperature_norm_K = 8000`) actually did.

---

**Pending, in order:**

1. **The sampling-noise floor.** Delete a random fraction of gas particles from
   the `cdm` sightlines, **multiply the surviving particle masses by
   `1/(1-f)`** so the mean gas density and `tau_eff` are untouched, and
   re-extract. The only thing that changes is the number of tracers, so what
   comes out is the white noise floor of the estimator. Use `f = 0.084`, the
   measured cdm-to-M3 gap. If that alone adds ~15% at k = 0.06, the residual is
   fully explained and there is nothing left to attribute. The uncompensated
   version (delete and stop) is a *different* experiment — it measures mass
   removal, not noise — and the 3D control already answered it (random removal
   gives 0.9999). Run the compensated one first; run the other only if the
   first does not close the gap.
2. **Move the threshold for real.** Re-run `cdm` and `M3` with
   `QLAStarFormation:over_density` at 1000 (baseline), 10^4, and with star
   formation off. If `R(k)` moves, the threshold is a numerical choice with an
   observable consequence and that is the paper. If it does not, the QLA line
   of work is dead and the repository still has the cross-validation. The box
   is 20 Mpc/h at 512^3 — cheap.
3. **Read `power_spectra/`, `statistics.txt`, `used_parameters.yml`** in the
   three murgia runs. Free, and one of them may replace a queue job.
4. **Error bars.** The two lowest k sit at 1.3 and 2.2 times the fundamental
   mode of a 20 Mpc/h box. 1536 sightlines are not 1536 independent samples of
   that mode. `t9` does bootstrap plus block jackknife and takes any two
   caches. Nothing at k < 0.007 s/km should be quoted until this runs.
5. **Re-run stage 02 on the `regen/*_uni512_seed12345.hdf5` pair.** Cheap.
   Decides whether the pairing was ever broken and what significance is
   quotable.
6. **V2 with Maria.** File is cut: `data/los_murgia_cdm_z5_20.hdf5`, 20
   sightlines, 16.9 MB, from `murgia/cdm/los_0003.hdf5` at z=5.0. Ask for
   `tau` shaped `[n_los, n_pix]` keyed by `LOS_XXXX`, plus her `n_pix`, `dv`,
   `Gamma_HI` and velocity convention. Our grid: 2048 px,
   `v_box = 2738.925 km/s`, `dv = 1.337366 km/s`, axis 2,
   `Gamma_HI = 4.30e-13 s^-1`. **Ask for tau, not flux.** Our answer is
   `cache/cache_murgia_cdm_z5_first100.npz`. Say the file is truncated when
   handing it over — it is a code-vs-code comparison, both sides get the same
   input.
7. **Their Figure 2.** `flux_spectra_NEW.pdf` is also vector and also
   extractable, and it carries the **MIKE/HIRES data points** at z = 4.2, 4.6,
   5.0, 5.4 together with their model curves in *absolute* P1D. Worth doing
   precisely because it tests the absolute normalisation, which is what just
   bit us — a ratio is blind to a common mean-flux error and Figure 2 is not.
   Do it after pendings 1-2, and expect it to be harder: 589 paths, error bars,
   and four redshifts overlaid.
8. **Fix `stages/05_xi.py` and `tests/t10_delta_p1d.py`.**
9. **Extract FCT/CDM caches at more redshifts.** Queue job. 17 LOS files
    exist, z = 5 to z = 1.8. `t11_bias_bk.py` refuses to interpret a single
    point.
10. **Stage 03** — `A_P` and `A_b` from the ICs, using Murgia+2019 Eqs. 1-5.
    Lower priority than it was: the ICs are already validated against their
    published curves, so this is now bookkeeping rather than a check. Watch
    the `h` convention: `(Mpc/h)^3` vs `Mpc^3` is `h^3 = 0.316`.
11. **`k_F(z)`, the filtering scale**, from the P_gas/P_matter suppression
    (Hui & Gnedin 1997), then `t11`.
12. **Transferability (GOAL stage C).** Repaint at three `tau_eff` and two
    `T0`. Note that pending 1 is already a weak version of this and the answer
    was *not* reassuring: `R(k)` moved by tens of per cent between
    normalisations. Quote the sensitivity honestly.
13. **t7**, on murgia as well as FCT. Cheap, brackets, gives a sign — but read
    Pato's critique below before quoting it as a number.
14. The 40 and 80 Mpc/h boxes are being re-run correctly by Pato. Fold in the
    15 new runs and `more_power` afterwards; their ICs were made **without**
    monofonIC's `masked=2`, so they are the "hard" case for t7 — run
    `stages/00_inspect_snapshot.py --deep` and read the verdict.

**Parked deliberately, both good, neither urgent:**

- **Profile over `A` instead of fixing it.** Grid over `(A_CDM, A_FCT)`, both
  free, and ask whether *any* choice reproduces `r(k) = 1`. Murgia+2019
  already does the marginalised version and it is what a referee will expect.
- **Dark gaps.** Distribution of contiguous runs with `F < threshold`. Runs on
  existing caches, is a real observable (Becker+2015, Zhu+2021), sensitive to
  the saturated tail P1D cannot see, and much less sensitive to continuum
  error than the raw flux PDF.

---

**Pato's critique of t7, and it is correct.** Putting the converted particles
back gives you particles that have been *collisionless* since they converted:
no pressure forces, so they over-collapsed, and they carry no temperature.
`t7` sweeps `T0` (open question 4), which addresses the temperature but not
the position or the velocity. **t7 gives a sign and an order of magnitude, not
a number.** The number needs the threshold moved for real.

---

### Correction to section 1

The section 1 table said the `tau_eff` rescaling "makes it slightly deeper".
**It makes it shallower.** `t0_rescaling` at z=3: `r(0.003)` goes 0.7338 raw
-> 0.7797 rescaled against the Turner target, or 0.7857 against the CDM run's
own opacity. Fixed in the table. The "not flat in k" claim (+6.26% at k=0.003
against +3.91% at k=0.03) holds **only for the Turner target** — see the
FCT/CDM section above, where the reference-run target gives a flat ~7% gain.
Mechanism: `A` multiplies `tau`, not `delta_F`, and
`|dF/dlnA| = A*tau*exp(-A*tau)` vanishes at both `tau << 1` and `tau >> 1`, so
`A` is a spatially selective gain, saturated pixels do not respond, and since
saturation traces large-scale structure the gain is k-dependent. This is the
mechanism that produced the phantom low-k deficit — it is real, it is just
much smaller when the target is sane.

---

### `legacy/relos.py` was edited — deliberately, see rule 3.2

`make_uniform_rays` wrote only `Xpos`/`Ypos`; SpectWizard needs `Xaxis`,
`Yaxis`, `Zaxis`. Separately, `write_output` copied `NumParts` verbatim from
the template. Both fixed. **No number changes:** our extractor reads
`Xpos`/`Ypos` and never reads `NumParts`, and the RNG stream is untouched.
Files written before that commit still need `legacy/fix_los_attrs.py` — the
production `regen/*_uni512_seed12345.hdf5` pair among them.

---

### The Murgia+2019 reference

**Murgia, Scelfo, Viel & Raccanelli 2019, PRL, arXiv:1903.10509.** Their PBH
grid is 512^3 in a 20 Mpc/h box — the same as our murgia runs, which are
29.52465 internal = 20 Mpc/h and walked 124.66M gas particles. Their data bins
are z = 4.2, 4.6, 5.0, 5.4; our LOS grid is 4.4 to 5.6 in steps of 0.2, so
three of four line up. GADGET-III + 2LPTic + CLASS against SWIFT + monofonIC +
CAMB. Their result: `f_PBH M_PBH < 170 M_sun` (2 sigma, flat `z_reio` prior),
`< 60` with a Gaussian prior.

Their equations, for stage 03: `P_PBH(k) = 1/n_PBH` (Eq. 1) with
`n_PBH = Omega_DM rho_cr f_PBH / M_PBH` (Eq. 2);
`P_iso = f_PBH^2 P_PBH = (2 pi^2 / k^3) A_iso (k/k_*)^(n_iso - 1)` (Eq. 4)
with `k_* = 0.05 /Mpc` and `n_iso = 4`; and
`P(k,z) = D^2(z) [T_ad^2 P_ad + T_iso^2 P_iso]` (Eq. 3).

Their reference simulation is 2x768^3 in 20 Mpc/h; the **PBH grid is 512^3**.
Reference thermal history `T_0(z=4.2) = 9200 K`, `gamma = 1.47`; cosmology
Planck-2015, `z_reio = 9` reference with `{7,9,15}` sampled; mean flux from
SDSS-III/BOSS with 9 rescalings of `F_bar` and 8 of `tau`.

Two things from it that bear on our framing:

- They **fit the full shape of the 1D flux power and marginalise over the mean
  flux**, with a Gaussian prior `sigma = 0.04` on `F_bar`. So "do not fix `A`,
  fit it" is what the reference analysis already does — cite it. It is also
  why the QLA claim has to be phrased in terms of the mean flux: that is the
  parameter it actually moves.
- Their Figure 1 makes the opposite-facing point to ours and both are true:
  non-linear evolution washes out the PBH signature in the 3D matter power, so
  the 1D flux power is the far more effective probe. Do not write "P1D is
  incomplete" in a way that reads as "P1D is weak". This is a Viel paper; the
  framing will be read carefully.

References worth having: **Irsic et al. 2017 (arXiv:1702.01764)** the
reference suite, **Murgia, Irsic & Viel 2018 (arXiv:1806.08371)** the method
paper, **Bolton et al. 2017** Sherwood, **Hui & Gnedin 1997** the filtering
scale.

---

**Non-obvious facts, kept:**
- **Normalise to the reference run's own raw `tau_eff`**, not to an external
  target, and say which value you used. This is the single most important
  operational lesson of the project so far.
- **Clementina is `pcolazo@ssh.clementinaxxi.org.ar` from the laptop**, not
  `snmgt01` (that is the internal node name and it does not resolve outside).
  Password auth only, so an assistant cannot run anything there — hand Pato
  the commands. Clementina also has no direct outbound internet; git goes
  through the HTTP proxy `172.28.3.3:3128`. DESI DR1 has to be downloaded on
  the laptop and `scp`-ed.
- **Figures in papers on arXiv are often vector and therefore exact.**
  `curl https://arxiv.org/e-print/<id>` then `pdftocairo -svg`. This replaced a
  "+-5-10 percentage points, do not quote" caveat with quotable numbers, and it
  is worth trying on any figure this project needs to compare against.
- **`relosz.py` is the entry point; `relos.py` is the engine.** relosz matches
  on the redshift each file *reports*, not the index in its name, and knows
  snapshots may live in their own subdirectory. It looks for snapshots as
  DIRECTORIES, so on a run that stores them as plain files beside the LOS files
  (the 40 Mpc/h pair) it reports "No hay snapshot a z=..." and you fall back to
  relos.py with an explicit `--snapshot` glob.
- **Snapshot and LOS indices do not correspond, in either direction**, and the
  offset is not consistent between runs. relos.py refuses a snapshot more than
  `dz = 0.02` away.
- **The murgia LOS grid**: z = 5.6, 5.4, 5.2, 5.0, 4.8, 4.6, 4.4 for
  `los_0000` through `los_0006`. The only (LOS, snapshot) pair within
  `dz <= 0.02` is **z=5.0** (`los_0003` + `murgia-*-lyman_0002`). The example
  in the header of `scripts/sbatch_roundtrip.sh` says 5.6, which dies on
  arrival. Stage 01 needs no snapshot, so all seven redshifts are extractable.
- **murgia is the only untruncated run on disk** until the 40/80 re-runs land.
  97.6-99.3% transverse coverage, ParticleIDs present, `frac_below_floor = 0`
  and `wsum_raw_med = 0.92` against 0.3145 and 0.755 for the truncated z=3
  file. A third of the pixels in the truncated files hit the Shepard floor.
- **Step 13 of `run_validation_A.sh` is buggy.** `CLEAN=$(head -1
  logs/v1_candidates.txt | cut -f2)` takes the *first* candidate, not a clean
  one. `data/los_clean_sample100.hdf5` is misnamed and should be deleted.
  Filter on the ray range, not on position.
- **`desi_window()` floats with redshift** (`k_max = 0.5 pi / R_z`). Use
  `common_window(zs)` for anything integrated across redshift bins.
- **`D^2(2.2)/D^2(4.0) = 2.40`**, from `units.growth_factor`.
- **`b` measured on these runs is `b` under QLA.** Say it out loud in the paper.
- **Mass-separability of the converted gas is per-run.** The 40 Mpc/h pair used
  monofonIC's `masked=2`, so DM and QLA-converted gas have different particle
  masses in `PartType1` — the easy case `--deep` detects. `more_power` and the
  15 new runs were made without the mask; treat them as the hard case.
- **Lead with results, not code.** Pato works as a supervisor on this
  repository. Give him verdicts, figures, and what would falsify them. Batch
  commands so one console dump answers several questions.
- **State confidence honestly, and check artefacts before mechanisms.** Two
  separate "results" this month turned out to be the same normalisation
  artefact wearing different clothes. Label a hypothesis as a hypothesis and
  name the command that would settle it.
- **The objectives document lives outside this repository, deliberately.** Ask
  Pato where. Decide the `.gitignore` question before it is ever committed.

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
and rises through unity at small scales. Measured 2026-09-03: **0.734 with no
normalisation and 0.786 rescaled to the CDM run's own opacity.** Quote the
pair; the conclusion does not depend on the choice.

**The whole question of this repository is whether that 0.78 is physics or a
measurement artefact.** Pato presented this at COSMO-26 in Leiden on 27
August 2026. Matteo Viel's objection there was that the normalisation to a
common effective optical depth could be manufacturing the result. The
repository exists to answer that publicly and reproducibly.

### What we already know

Five explanations were tested to completion and all failed:

| candidate | test | verdict |
|---|---|---|
| rescaling to a common `tau_eff` | A = 1 in both runs | discarded, and **re-confirmed 2026-09-03 against the reference run's own opacity**: `r(0.003)` is 0.7338 raw, 0.7857 rescaled. The deficit survives both. Rescaling makes it *shallower*, and against the correct target the gain is flat in k (~7% across the window), not the 2.35 pp tilt the Turner target produced. See section 0 |
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
   single most valuable unknown right now. **Answered for murgia, 2026-09-02:
   `PartType4` is empty and the converted gas is in `PartType1`.** See
   section 0.
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
   **And see Pato's critique in section 0**: temperature is not the only
   thing wrong with a reinjected star particle.
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
- `no snapshots found under ...` from `stages/00_inspect_snapshot.py` → it
  does not recurse. Point `--run` at the snapshot subdirectory itself.
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

- **The laptop has no git credential helper configured, and no `gh`.** This
  has come up repeatedly. `git ls-remote` works (read is anonymous), but
  `git push` from the laptop would prompt for a token interactively, which
  is not something to do in a Claude session. **From the laptop, push with
  the GitHub MCP** (`create_or_update_file`, one call per file, full file
  content, needs the current blob SHA). It goes straight to `main`. Verify
  afterwards by `curl`-ing the raw file back and diffing against the local
  copy — the MCP call takes the whole file, so a transcription slip is
  silent otherwise. **This has already destroyed this file once**: on
  2026-09-03 a push landed the literal string `PLACEHOLDER` as the entire
  HANDOFF. Diff after every push of a large file, without exception. Note
  that `raw.githubusercontent.com` is CDN-cached for a few minutes, so
  verify through
  `https://api.github.com/repos/.../contents/<path>?ref=main` with the
  `Accept: application/vnd.github.raw` header instead. From Clementina,
  plain `git push` over SSH through the proxy works; that is the route for
  anything large, and it is the one to prefer.
- Whoever is about to edit code, `git pull` first (or, on the laptop side,
  pull the current file with
  `curl -sfL https://raw.githubusercontent.com/patricio-c/Lyman-Alpha-PBHs/main/<path>`
  and edit that, rather than editing a copy from earlier in the session —
  the MCP push always goes straight to `main`).
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

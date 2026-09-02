# HANDOFF — read this first

You are picking up a scientific analysis repository. This file is the whole
context. Read it before touching anything.

The person you are working with is **Pato Colazo**, PhD student in astronomy
at IATE/OAC (Universidad Nacional de Córdoba, CONICET). He works in
Rioplatense Spanish and writes papers in English. He is on a Linux laptop;
the cluster is remote.

---

## 0. Current status — read this first, every session

Updated 2026-09-02. This section is the fast-changing punch list; the rest of
the file is stable background. Start every new session here, and dip into the
numbered sections below only for the *why* behind something.

**The headline, if you read nothing else.** The validation phase is over: V1
passed and the chain is under control. What came out of it is a second,
larger result that was not the plan. The star-formation threshold used by
every Lyman-alpha constraint in this literature removes a *different amount
of gas* in a model with enhanced small-scale power than it does in LCDM, and
that differential has never been reported. We measured it from the snapshot
headers: at `Delta > 1000`, the LCDM run converts **7.12%** of its baryons and
the `M_PBH = 10^3 M_sun` run converts **15.51%** — more than double, a gap of
8.4 percentage points. The approximation was validated inside LCDM only, and a
constraint needs the differential, not the approximation. See "The QLA
threshold result" below. That is now the most valuable thing in this
repository.

---

### Done this session (2026-09-01/02)

- **V1 PASSED.** `relos.py` reproduces SWIFT's own LOS output to numerical
  precision. Details below. t7 is unblocked.
- **The `uni512_seed12345` finding.** Closed three separately-logged
  mysteries, reopened one. Details below.
- **`legacy/relos.py` fixed** to write `Xaxis`/`Yaxis`/`Zaxis`/`NumParts`.
  No number changes; ray positions verified bit-identical.
- **`common/units.py` numpy-2 fix** (Pato). `getattr(np, "trapezoid",
  np.trapz)` evaluated its default eagerly and died on numpy >= 2.3.
  `tests/test_estimator.py` now runs to the end.
- **Nine murgia caches extracted** — `cdm`, `M2`, `M3` at z = 4.6, 5.0, 5.4,
  1536 sightlines each, `scripts/sbatch_murgia.sh`. First V3 numbers below.
- **The ratio measured under both normalisations** (`--norm taueff` and
  `--norm none`), and the **exact converted-gas fractions** read out of the
  z=5 snapshot headers. Those two together are what produced the QLA result
  and both corrections to it.
- **The V2 hand-off file is cut and verified:**
  `data/los_murgia_cdm_z5_20.hdf5`, 20 sightlines, 16.9 MB, from
  `murgia/cdm/los_0003.hdf5` at z=5.0. Carries `NumParts`, `Xaxis`, `Xpos`,
  `Yaxis`, `Ypos`, `Zaxis` per group. Not sent yet.

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
numbers, not the verdict.** The three extra particles sit in pixels with
`tau ~ 4000`, i.e. `F = 0` on both sides.

The sentence this licenses:

> Our regeneration of SWIFT's line-of-sight output reproduces the flux field
> and the 1D flux power spectrum to numerical precision. The particle sets
> agree on 97 of 100 sightlines; the three exceptions differ by one particle
> each, in fully saturated pixels, and change P1D by less than 1e-4 per cent.

**Scope**: one run, one redshift, 100 sightlines, all along axis 2, reusing
the template's ray positions. Not tested: the `--uniform` path, `--skip-los`
batching, axis detection on axes 0 and 1. Enough for t7.

**The three extra particles are not understood and it does not matter.**
`b/(gamma*h)` = 0.495, 0.879, 0.766, so not a kernel-edge rounding case.
`SplitCounts = 0`, progenitor == own ID, so not particle splitting.
`Header/Time` bit-identical, so not a different epoch. Their `h` is in the
bottom percentile of their sightline, so not a stale cell `h_max`. Four
hypotheses, four wrong. If it is ever worth chasing, the answer is in SWIFT's
`line_of_sight.c`, not in more numpy.

---

### The `uni512_seed12345` finding

**The production caches were NOT made from SWIFT's own LOS files.** From the
`source` field of the caches:

    cache/cache_cdm.npz  (1536, 2048)  regen/cdm40_z3.0_uni512_seed12345.hdf5
    cache/cache_fct.npz  (1536, 2048)  regen/fct40_z3.0_uni512_seed12345.hdf5

`relos.py --uniform 512` draws 512 rays per axis over three axes = **1536**,
uniformly over the **whole** box face. Verified bit-for-bit: regenerating with
`--uniform 512 --seed 12345` on the 58.73715 box gives `LOS_0000  Xaxis=1
Xpos=13.353070342258142 Yaxis=2 Ypos=18.605482517633362 Zaxis=0`, identical to
the file on COSMA. Confirmed by Pato: that is how they were made.

1. **"6144 or 1536 sightlines?" was a false dichotomy.** Different files, not
   a subsample.
2. **The published numbers do NOT rest on truncated sightlines.** The 6144
   SWIFT files are truncated to 46% of the transverse face by the `h` mix-up;
   the 1536 uniform rays cover the whole face by construction.
3. **`tau_eff` 0.42461 (cache) vs 0.26192 (100-line cut) is not a bug.** The
   cut came from the truncated SWIFT file, missing ~32% of particles per line.

**And stage 02 measured the wrong files.** It was pointed at
`lyman/{cdm,2-fct}-box-40-1024/los_0010.hdf5`, the 6144 SWIFT files, which are
genuinely unrelated to each other. The caches come from the 1536 uniform
regenerations, which — same seed, same box, same `default_rng` — are the same
positions in both runs. That explains t9's `sigma(common)/sigma(unpaired) =
0.3733` and per-line correlation `+0.8592`, which is impossible for unpaired
lines. **Re-run stage 02 on the `regen/*_uni512_seed12345.hdf5` pair.** Until
then, do not requote any significance and do not send anyone a correction to
the 20 sigma.

---

### THE QLA THRESHOLD RESULT — the new line of work

**Read the two corrections at the end of this section before using anything
here.** The first reading of these numbers was wrong twice in a row, and the
version that survived is better than the one that did not.

**The measurement, exact, from the z=5 snapshot headers.** Initial gas is
`512^3 = 134217728` particles:

    run   gas at z=5     converted    (M_PBH, f_PBH = 1)
    cdm   124664686       7.12%       LCDM
    M2    123838499       7.73%       10^2 M_sun
    M3    113403060      15.51%       10^3 M_sun

`M3` converts **more than twice** what `cdm` does. That 8.4 percentage-point
gap is the number this line of work is about. All three share the Panphasia
descriptor `[Panph6,L20,(235287,445214,422255),S1,KK1025,CH-999,COLIBRE050]`
and the same `Omega_cdm`: same realisation, different input transfer function,
so sample variance cancels in the ratio at the level of the ICs.

Counted near the 1536 sightlines instead of globally, M3/cdm is 14.2% rather
than 9.0% at z=5. That difference is not noise — the surviving gas in M3 is
distributed differently, which is the bias change showing up directly.

**`PartType4` is empty in all three runs**, yet gas is missing. SWIFT's QLA
sends the converted gas to `PartType1` as collisionless particles, not to
`PartType4` — which is exactly what Viel+2004 says it does ("converted to
collisionless particles"), and consistent with what this file already records
about `masked=2` giving DM and converted gas different masses inside
`PartType1`. So do not count stars to get the conversion fraction; count the
missing gas, or split `PartType1` by mass. This also answers open question 1
for murgia.

**The prescription is nominally identical to the literature's.** From the
SWIFT parameters stored in the LOS file:

    Parameters/QLAStarFormation:over_density = 1000

And the chain: Murgia+2019 -> Irsic+2017 -> Bolton+2016 (Sherwood) ->
Viel+2004, all `Delta > 1000`. **But Viel adds `T < 10^5 K` and we appear not
to**, and SWIFT also carries an entropy floor
(`QLAEntropyFloor:over_density_threshold = 10`, `temperature_norm_K = 8000`)
that GADGET does not. So "same threshold" is not "same converted fraction",
and that asymmetry is probably why our differential is large. **Confirm it in
`used_parameters.yml` in the run directory, not in the LOS attributes** — see
the untapped files below.

**Where the justification fails.** Viel, Haehnelt & Springel 2004, which the
whole chain rests on:

> "In a pixel-to-pixel comparison with a simulation which adopted the full
> multi-phase star formation model of Springel & Hernquist (2003) we
> explicitly checked for any differences introduced by this approximation. We
> found that the differences in the flux probability distribution function
> were smaller than 2%, while the differences in the flux-power spectrum were
> smaller than 0.2%."

That test is **QLA versus multi-phase, inside LCDM**. It validates the
approximation for a given cosmology and says nothing about
`QLA(model A) - QLA(model B)`, which is what a constraint needs — and the
entire mechanism of these models is to send more gas over the threshold. In
their check both sides had the same amount of gas crossing. In a PBH
constraint they do not.

**Why the "that gas makes no P1D" defence does not hold.** True for its
*direct* contribution: it is self-shielded and saturated, `F = 0`. False for
its *indirect* effect, and this repository already has the controls (section
1): removing 40% of the gas **at random** gives a power ratio of 0.9999;
removing gas **above a density threshold** until the masses match drops it to
about 0.10. Same mass, two orders of magnitude apart. Threshold removal is
spatially correlated with the density field, so it changes the gas-matter bias
on *all* scales, including the unsaturated ones that do produce P1D.

**The measured ratio, both normalisations.** `stages/04_p1d.py`, ratio to
`cdm`, z=5.0, and the same k converted to h/Mpc (x136.9463):

    k [s/km]  k [h/Mpc]   M3 norm=taueff   M3 norm=none   Murgia+19 Fig.1
      0.0030      0.41         -3.5%           +9.9%          ~0
      0.0050      0.68         -5.7%          +10.6%          ~0
      0.0100      1.37         +5.0%          +22.0%          ~3
      0.0200      2.74         +9.7%          +32.2%          ~9
      0.0300      4.11        +18.8%          +53.8%         ~16
      0.0600      8.22        +43.3%         +120.1%         ~33

The Figure 1 column is **eyeballed off a rendered page**, worth +-2 points at
the bottom and +-8 at the top. Do not quote it. Get the real curve, or compute
the *linear* matter curve exactly from their Eqs. 1-5 — that one is analytic
and needs no digitising, and it is also stage 03.

**CORRECTION 1: the low-k deficit is produced by the rescaling.** Under
`--norm none` the ratio is positive at every k and every redshift, rising
monotonically, which is what a positive-definite Poisson term must give:

    k [s/km]   M3 z=4.6   M3 z=5.0   M3 z=5.4
      0.0030    1.1070     1.0991     1.0387
      0.0050    1.0855     1.1061     1.0485
      0.0100    1.1947     1.2196     1.1492
      0.0600    1.7873     2.2010     2.4223

So the deficit is not an independent signature of QLA in the flux field. It is
what the common-`tau_eff` rescaling does when the two runs sit far apart in
raw `tau_eff`.

**CORRECTION 2: and yet Murgia+2019 rescale too, so the deficit is still the
anomaly.** An earlier draft of this section argued their Figure 1 is at fixed
UVB. The numbers say otherwise: the **normalised** ratio tracks their curve
within the reading error from 1 to 8 h/Mpc, while the **unnormalised** one
overshoots by a factor 3-4. Rescaling to a common mean flux is also the
standard in this literature, because the mean flux is what is observed.

**So the chain is this, and it is stronger than the version it replaces:**

> more small-scale power -> more gas over `Delta = 1000` -> more gas removed
> -> the two runs' raw `tau_eff` move apart (2.843 vs 3.008 at z=5) -> their
> fitted `A` move apart (0.37877 vs 0.30997, an 18% gap) -> **the
> common-mean-flux rescaling that everyone performs then imprints a
> k-dependent tilt whose size is set by the star-formation differential.**

Nobody made an error. The standard procedure is correct and still imprints a
model-dependent distortion. Why it bites us and not them is the size of the
lever, and the lever is the converted-gas differential.

**And it makes a falsifiable prediction**: raise the threshold, the converted
fractions converge, the `A` converge, and **the low-k deficit must shrink**.
That is pending item 2, and it now has a number to hit instead of being an
open-ended check.

**Two things still missing before any of this is a result:**

1. **There are no error bars.** The two lowest k sit at 1.3 and 2.2 times the
   fundamental mode of a 20 Mpc/h box (`k_fund = 0.002294 s/km` at z=5, i.e.
   a wavelength 0.76 of the box). The 1536 sightlines are not 1536 independent
   samples of that mode. Shared phases between the runs help a lot; unmatched
   sightline positions between the runs do not. Only a jackknife settles it,
   and `t9` already does bootstrap plus block jackknife and takes any two
   caches.
2. **Sampling noise, and which way it pushes.** `wsum_raw_med` degrades with
   the model: cdm 0.930/0.924/0.918, M2 0.900/0.888/0.874, M3
   0.800/0.771/0.731. Fewer gas particles means a noisier SPH density
   estimate; shot noise is uncorrelated between pixels, so it is **white, flat
   in k**, and P1D falls steeply with k. It therefore **inflates the high-k
   enhancement** — the part that validates against the paper — and **cannot
   create the low-k deficit**. Measure it by deleting a random 9% of gas
   particles from `cdm` and re-extracting; `legacy/removal_curve_pk.py` is the
   3D version of the same control.

**Pato's critique of t7, and it is correct.** Putting the converted particles
back gives you particles that have been *collisionless* since they converted:
no pressure forces, so they over-collapsed, and they carry no temperature.
`t7` sweeps `T0` (open question 4), which addresses the temperature but not
the position or the velocity. **t7 gives a sign and an order of magnitude, not
a number.** The number needs the threshold moved for real. Keep t7 — it is
cheap and it brackets — but do not quote it as the systematic.

**Untapped, and worth more than another extraction.** Every murgia run
directory holds files nobody has opened yet:

- `power_spectra/` — SWIFT computed **3D power spectra on the fly**. If
  `P_matter` and `P_gas` are in there for all three runs, that is the
  measurement that settled the FCT argument (matter 1.000, baryons 1.000, gas
  0.808) available for free, without Pylians and without a queue. Do this
  first.
- `statistics.txt` — mass per particle type against time, i.e. the **full
  conversion history**, not one number at z=5. Turns the differential from a
  point into a curve against redshift.
- `SFR.txt` — the same thing as a rate.
- `used_parameters.yml` — what SWIFT **actually applied**, defaults included.
  Cite this one, not `M3.yml`, which is only the input.
  `unused_parameters.yml` shows anything SWIFT silently ignored.

---

### Correction to section 1

The section 1 table said the `tau_eff` rescaling "makes it slightly deeper".
**It makes it shallower.** `t0_rescaling` at z=3: `r(0.003)` goes 0.7338 raw
-> 0.7797 rescaled. Fixed in the table.

And the rescaling is **not flat in k**: `r(k)` moves +6.26% at k=0.003 s/km
and +3.91% at k=0.03, a 2.35 percentage-point tilt across the window, and the
slope of `r(k)` flattens by 2.2%. Quote it as a systematic, do not deny it.
Mechanism: `A` multiplies `tau`, not `delta_F`, and
`|dF/dlnA| = A*tau*exp(-A*tau)` vanishes at both `tau << 1` and `tau >> 1`.
So `A` is a spatially selective gain; saturated pixels do not respond; and
since saturation traces large-scale structure, the gain is k-dependent. After
rescaling, FCT has 21.1x more transparent pixels (F>0.99) than CDM — the two
runs lean on the F=1 wall by different amounts.

---

### `legacy/relos.py` was edited — deliberately, see rule 3.2

`make_uniform_rays` wrote only `Xpos`/`Ypos`; SpectWizard needs `Xaxis`,
`Yaxis`, `Zaxis` and aborts without them, which is why `legacy/fix_los_attrs.py`
exists. Separately, `write_output` copied `NumParts` verbatim from the template,
so a regenerated file carried the *original's* particle count.

Both fixed. **No number changes:** our extractor reads `Xpos`/`Ypos` (the log
prints `ray position source: ['attrs']`) and never reads `NumParts`, and the
RNG stream is untouched. Files written before this commit still need
`legacy/fix_los_attrs.py` — the production `regen/*_uni512_seed12345.hdf5`
pair among them.

---

### V3 status

t0 clean, stage 04 clean, t9 ran without the block jackknife. Two pieces are
broken and must be fixed before they are read: **stage 05** prints ratios
through a zero crossing (`xi_FCT(600) = 0.00015`, so the 0.0488 and 1.6129 in
the log are noise — plot the difference, not the ratio), and **t10**
contradicts itself and fits above the window (`k < 0.0653 s/km` against a DESI
top of 0.0319 and a `common_window` top of 0.0255).

**A V2 failure does not invalidate all of V3.** A purely multiplicative error
in `tau` is absorbed entirely by `A` and leaves `r(k)` untouched, because both
runs go through the same extractor. An error that depends on density or
temperature does move it, because that is exactly where the models differ.
V2 tells us which case we are in.

---

### The Murgia+2019 reference

**Murgia, Scelfo, Viel & Raccanelli 2019, PRL, arXiv:1903.10509**, "Lyman-alpha
forest constraints on Primordial Black Holes as Dark Matter". Their PBH grid is
512^3 in a 20 Mpc/h box; our murgia files are 29.52465 internal = 20 Mpc/h and
walked 124.66M gas particles = 93% of 512^3. Their data bins are z = 4.2, 4.6,
5.0, 5.4; our LOS grid is 4.4 to 5.6 in steps of 0.2, so three of four line up.
They ran GADGET-III + 2LPTic + CLASS where we run SWIFT + monofonIC + CAMB, so
this is a genuine code-independent comparison. Their result:
`f_PBH M_PBH < 170 M_sun` (2 sigma, flat `z_reio` prior), `< 60` with a
Gaussian prior.

**Use their equations for stage 03 rather than re-deriving them:**
`P_PBH(k) = 1/n_PBH` (Eq. 1) with `n_PBH = Omega_DM rho_cr f_PBH / M_PBH`
(Eq. 2); `P_iso = f_PBH^2 P_PBH = (2 pi^2 / k^3) A_iso (k/k_*)^(n_iso - 1)`
(Eq. 4) with `k_* = 0.05 /Mpc` and `n_iso = 4`; and
`P(k,z) = D^2(z) [T_ad^2 P_ad + T_iso^2 P_iso]` (Eq. 3).

Two things from it that bear on our framing:

- They fit **the full shape of the 1D flux power, not a single amplitude
  parameter**, and marginalise over the mean flux (9 grid points,
  `{0.6...1.4} x F_REF`) plus 8 further `tau` rescalings. So "do not fix `A`,
  fit it" is what the reference analysis already does. Cite it.
- Their Figure 1 makes the opposite-facing point to ours and both are true:
  non-linear evolution washes out the PBH signature in the **3D matter**
  power, so the **1D flux** power is the far more effective probe. Do not
  write "P1D is incomplete" in a way that reads as "P1D is weak". The honest
  claim is that P1D is a second moment of a bounded, strongly non-Gaussian
  field and carries no information from saturated pixels, so tail-sensitive
  statistics are *complementary*. This is a Viel paper; the framing will be
  read carefully.

Their references worth having: **[68] Irsic et al. 2017 (arXiv:1702.01764)**
is the reference simulation suite, **[63] Murgia, Irsic & Viel 2018
(arXiv:1806.08371)** the method paper, **[74] Bolton et al. 2017** is Sherwood,
**[73] Hui & Gnedin 1997** is the filtering scale for pending item 9.

---

**Pending, in order:**

1. **Read `power_spectra/`, `statistics.txt` and `used_parameters.yml` out of
   the murgia run directories.** No compute at all, and it may hand over the
   3D gas-field measurement and the conversion history that the argument
   needs. See the end of the QLA section.
2. **Move the threshold for real.** Re-run `cdm` and `M3` with
   `QLAStarFormation:over_density` at 1000 (baseline), 10^4, and with star
   formation off. If `R(k)` moves with the threshold then `R(k)` is not an
   observable, it is a function of a numerical choice. It now has a specific
   prediction to hit: raising the threshold must shrink the low-k deficit.
   This is the decisive test and it is the paper. The box is 20 Mpc/h at
   512^3 — cheap next to the 1024^3 runs.
3. **Error bars on the murgia ratio**, via `t9` on those caches, and the
   random-deletion noise floor (drop a random 9% of gas from `cdm`,
   re-extract). Neither number above is interpretable without them.
4. **Re-run stage 02 on the `regen/*_uni512_seed12345.hdf5` pair.** Cheap.
   Decides whether the pairing was ever broken, whether t8 is valid, and what
   significance gets quoted.
5. **V2 with Maria.** File is cut. Ask for `tau` shaped `[n_los, n_pix]` keyed
   by `LOS_XXXX`, plus her `n_pix`, `dv`, `Gamma_HI` and velocity convention.
   Our grid: 2048 px, `v_box = 2738.925 km/s`, `dv = 1.337366 km/s`, axis 2,
   `Gamma_HI = 4.30e-13 s^-1`. **Ask for tau, not flux** — `tau_eff = 2.85`
   raw, most pixels saturated, a 30% tau disagreement inside a trough is
   invisible in `F`. Our answer is `cache/cache_murgia_cdm_z5_first100.npz`.
6. **Fix stage 05 and t10** (see V3 status).
7. **Extract FCT/CDM caches at more redshifts.** Queue job. The whole argument
   is redshift evolution and there is one redshift on disk. 17 LOS files
   exist, z=5 to z=1.8. `t11_bias_bk.py` refuses to interpret a single point.
8. **Stage 03** — `A_P` and `A_b` from the ICs, using Murgia+2019 Eqs. 1-5.
   Measure the z=198 box with Pylians and check the white-noise plateau sits
   where `A_P` predicts. This is also an exact check on the murgia ICs: the
   *linear* matter curve in their Figure 1 is analytic, so it can be compared
   without digitising anything. Watch the `h` convention: `(Mpc/h)^3` vs
   `Mpc^3` is `h^3 = 0.316`.
9. **`k_F(z)`, the filtering scale**, from the P_gas/P_matter suppression
   already measured (Hui & Gnedin 1997). The most fragile input in the
   programme, because the broken term scales as `k_max^4`.
10. `t11`, once 8 and 9 exist.
11. **Transferability (GOAL stage C).** Repaint at three `tau_eff` and two
    `T0`; `R(k)` must move at the percent level or the ratio is not the
    observable. `--impose-trho` exists in `legacy/prepatch/`.
12. **t7**, on murgia as well as FCT. Cheap, brackets the systematic, gives a
    sign — but read Pato's critique above before quoting it as a number.
13. The 40 and 80 Mpc/h boxes are being re-run correctly by Pato (the `h`
    mix-up truncates **every** box in that family: `cdm-box-80` and
    `fct-box-80` also stop at 80 of 117.4743, 45.6% coverage). Fold in the 15
    new runs and `more_power` afterwards; their ICs were made **without**
    monofonIC's `masked=2`, so they are the "hard" case for t7 — run
    `stages/00_inspect_snapshot.py --deep` and read the verdict.

**Parked deliberately, both good, neither urgent:**

- **Profile over `A` instead of fixing it.** Grid over `(A_CDM, A_FCT)`, both
  free, and ask whether *any* choice reproduces `r(k) = 1`. It cannot: the tilt
  `A` can inject over the window is a few per cent while `r(k)` swings 43%.
  Two quotable statements — the systematic (vary `A` within the observed
  `tau_eff` error bar) and the robustness (no `A` imitates the effect).
  Murgia+2019 already does the marginalised version.
- **Dark gaps.** Distribution of contiguous runs with `F < threshold`. Runs on
  existing caches, is a real observable (Becker+2015, Zhu+2021), is directly
  sensitive to the saturated tail P1D cannot see, and is much less sensitive
  to continuum error than the raw flux PDF.

**Non-obvious facts, kept:**
- **`stages/00_inspect_snapshot.py` does not recurse into subdirectories.**
  It reports "no snapshots found" on a run that keeps them in their own
  directory, which is most of them. `common/runs.py::list_hdf5` does recurse,
  so the two discovery paths disagree. Workaround: point `--run` at the
  snapshot subdirectory itself
  (`.../murgia/M3/murgia-M3-lyman_0002`). Worth fixing to use `list_hdf5`.
- **SWIFT's QLA converts gas to `PartType1`, not `PartType4`.** `PartType4` is
  empty in the murgia runs even though 7-16% of the gas is gone. Count the
  missing gas or split `PartType1` by particle mass; do not count stars.
- **Clementina has no direct outbound internet.** git goes through the HTTP
  proxy `172.28.3.3:3128`. See section 6. DESI DR1 has to be downloaded on the
  laptop and `scp`-ed.
- **`relosz.py` is the entry point; `relos.py` is the engine.** relosz matches
  on the redshift each file *reports*, not the index in its name, and knows
  snapshots may live in their own subdirectory. It looks for snapshots as
  DIRECTORIES, so on a run that stores them as plain files beside the LOS
  files (the 40 Mpc/h pair) it reports "No hay snapshot a z=..." and you fall
  back to relos.py with an explicit `--snapshot` glob.
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
  one, so it ran on the truncated `2-fct-box-40-1024/los_0015.hdf5` under an
  echo claiming otherwise. `data/los_clean_sample100.hdf5` is misnamed and
  should be deleted. Filter on the ray range, not on position.
- **`desi_window()` floats with redshift** (`k_max = 0.5 pi / R_z`). Use
  `common_window(zs)` for anything integrated across redshift bins.
- **`D^2(2.2)/D^2(4.0) = 2.40`**, from `units.growth_factor`. The calculable
  part of the redshift lever; everything above it is filtering scale, thermal
  state and QLA.
- **`b` measured on these runs is `b` under QLA.** 49.8% of FCT baryons
  converted against 14.1% in CDM at z=3. Say it out loud in the paper.
- **Mass-separability of the converted gas is per-run.** The 40 Mpc/h pair used
  monofonIC's `masked=2`, so DM and QLA-converted gas have different particle
  masses in `PartType1` — the easy case `--deep` detects. `more_power` and the
  15 new runs were made without the mask; treat them as the hard case.
- **Lead with results, not code.** Pato works as a supervisor on this
  repository. Give him verdicts, figures, and what would falsify them. Batch
  commands so one console dump answers several questions. When a stage passes,
  say which files he now has to read and which stay machinery.
- **State confidence honestly.** Four mechanism guesses were wrong in one
  session (the three extra particles), and the QLA section above had to be
  corrected twice. Label a hypothesis as a hypothesis and name the command
  that would settle it. Pato reads the reasoning, not just the conclusion.
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
  silent otherwise. Note that `raw.githubusercontent.com` is CDN-cached for
  a few minutes, so verify through
  `https://api.github.com/repos/.../contents/<path>?ref=main` with the
  `Accept: application/vnd.github.raw` header instead. From Clementina,
  plain `git push` over SSH through the proxy works; that is the route for
  anything large.
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

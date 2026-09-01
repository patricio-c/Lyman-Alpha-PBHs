# HANDOFF — read this first

You are picking up a scientific analysis repository. This file is the whole
context. Read it before touching anything.

The person you are working with is **Pato Colazo**, PhD student in astronomy
at IATE/OAC (Universidad Nacional de Córdoba, CONICET). He works in
Rioplatense Spanish and writes papers in English. He is on a Linux laptop;
the cluster is remote.

---

## 0. Current status — read this first, every session

Updated 2026-09-01. This section is the fast-changing punch list; the rest
of the file is stable background. Start every new session here, and dip
into the numbered sections below only for the *why* behind something.

**Done:**
- Repo scaffolding (`common/`, `stages/`, `tests/`, `scripts/`, `paper/`,
  docs) pushed to `https://github.com/patricio-c/Lyman-Alpha-PBHs` (public).
- Migration done 2026-08-31: `legacy/` populated on Clementina,
  `VERIFY OK`, committed, pushed. The old `LOS/` directory can be deleted
  whenever — both conditions are met, no rush.
- Geometry check passed clean 2026-08-31. Every number in section 5 matched
  exactly. No unit bug in this pair.
- LOS-match check run 2026-08-31 and FAILED. The sightlines are not paired.
  See "Non-obvious facts".
- Resampling machinery and three new tests added 2026-08-31:
  `common/boot.py`, `tests/t9_unpaired_significance.py`,
  `tests/t10_delta_p1d.py`, `tests/t11_bias_bk.py`. `common/units.py`
  gained `growth_factor(z)` and `common_window(zs)`.
  `tests/test_estimator.py` now covers all of it and passes. It needs no
  data and takes a second — run it before anything else.
- Validation tooling added 2026-08-31, corrected against the real disk
  layout 2026-09-01: `tests/t12_relos_roundtrip.py`,
  `scripts/make_los_subset.py`, `scripts/run_validation_A.sh` (block A, one
  log) and `scripts/sbatch_roundtrip.sh` (V1, queue). See "the validation
  phase" below, which is what the project is actually doing right now.

**THE VALIDATION PHASE — this is what is happening now.**

Pato's call, 2026-08-31, and it governs the order of everything below:
nothing new gets built until the existing chain is checked against
something independent. He is deliberately not reading the analysis code
yet; he is asking for figures and verdicts, and will read the code once
the results come out right. So the job is to produce results that are
checkable without reading the code, and to say plainly what each one
would look like if it were broken.

Three legs, run in parallel:

  V1. **Do we control the sightlines?** Our regeneration must reproduce a
      LOS file SWIFT wrote itself, particle for particle, from the snapshot
      at the same redshift. `tests/t12_relos_roundtrip.py` is the check.
      Needs a run whose LOS output is not truncated, which rules out the
      40 Mpc/h pair and points at murgia — see the facts below. Step 4 of
      `scripts/run_validation_A.sh` finds the candidates and step 4b prints
      each one's ray range so a truncated candidate is visible before the
      queue time is spent on it. Run it through
      `scripts/sbatch_roundtrip.sh <run dir> <z>`, which goes through
      `legacy/relosz.py`. Gates t7.
  V2. **Is tau right?** Our extractor against SpectWizard on the same 100
      sightlines. `scripts/make_los_subset.py` cuts the file to send Maria.
      This is HANDOFF open question 5 and it is the only external check in
      the whole project. Gates every quoted number.
  V3. **Is the deficit real?** Finish the tests already designed: t0,
      stage 04, stage 05, t9, t10. These run on existing caches in minutes
      and need nothing from V1 or V2.

**Pending, in order:**
1. **Extract LOS caches at more redshifts.** Queue job. This blocks
   everything else and it is not close: the whole argument is redshift
   evolution, and there is exactly one redshift on disk. 17 LOS files
   exist from z=5 to z=1.8. `t11_bias_bk.py` refuses to interpret a single
   point and says so in its own output.
2. Run `tests/t9_unpaired_significance.py` on the existing z=3 caches. Pass
   `--los-a`/`--los-b`: without the sightline positions it skips the block
   jackknife, which is the part worth having. Cheap, no queue.
3. Run `tests/t10_delta_p1d.py` on the same caches. Also cheap.
4. **Stage 03 — `A_P` and `A_b` from the initial conditions.** Does not
   exist yet and `t11` cannot run without it. The Poisson amplitude is
   analytic (`P_Poisson = 1/n_PBH` in Mpc^3, from `f_PBH` and the mass
   function); `A_b` comes from the Sureda et al. parametrisation. Computing
   them is not enough — measure the z=198 box with Pylians and check the
   white-noise plateau sits where `A_P` predicts. Watch the `h` convention:
   `(Mpc/h)^3` vs `Mpc^3` is a factor `h^3 = 0.316`, which is why `t11`
   has no default for `--units`.
5. **`k_F(z)`, the filtering scale, for `k_max`.** No new runs needed: fit
   it from the P_gas/P_matter suppression already measured (Gnedin & Hui).
   It is the most fragile input in the programme because the broken term
   scales as `k_max^4`.
6. `t11` for real, once 4 and 5 exist.
7. Transferability test (GOAL stage C): repaint at three `tau_eff` and two
   `T0`, check `R(k)` moves at the percent level. Hours of cache, not of
   compute, and it is a gate — if `R(k)` is not stable the ratio is not the
   observable and the strategy changes. `--impose-trho` already exists in
   `legacy/prepatch/`.
8. `t7`. Its priority went up, see the QLA fact below.
9. Steps 4-6 of section 5 (the analysis on existing caches, the paper) are
   unchanged. Step 4 can proceed at any time — it never assumed pairing.
10. Two more datasets to fold in once the CDM/FCT pair is done — see
    section 4:
    - 15 additional runs, extracted with correct on-the-fly LOS output.
      Their ICs were made *without* monofonIC's `masked=2`, so treat them
      like the `more_power` case for `t7`, not like the original pair.
    - `lyman/murgia/{cdm,M2,M3}` — a different snapshot/LOS grid (3
      snapshots, 7 LOS files, not the usual 5/17), and a different box:
      29.52465 internal, not 58.7372. `murgia/cdm/los_0000.hdf5` is z=5.600,
      measured. Do not assume any other `los_000N` lines up with a
      particular z, and do not assume a LOS index matches a snapshot index —
      step 4 of `run_validation_A.sh` reads them all.

**Decisions taken 2026-08-31, and why:**
- **Re-quote the significance before deciding whether to re-shoot the LOS.**
  Two options were on the table: re-run the SWIFT LOS output for both boxes
  with a shared seed/position list, or re-derive the 20 sigma as an unpaired
  test. They are not alternatives. The re-quote is hours on caches that
  already exist and it is *diagnostic*: it says whether the pairing was
  load-bearing at all. Re-shooting is days of queue plus the risk of a new
  discrepancy against the caches everything else is built on, and `t7` needs
  sightlines re-shot anyway, so bundle it there rather than running a
  separate campaign. Order: re-quote first, re-shoot with `t7`.
- **Do not assume the unpaired number comes out weaker.** A bootstrap
  *measures* the covariance between the two samples, it does not assume one.
  Resampling a common index set and differencing gives
  `var_A + var_B - 2cov`, and when the lines are unrelated the `cov` in the
  data is already zero. A paired bootstrap over a fake pairing should
  therefore degrade to the unpaired answer by itself. That is a claim about
  the data, so `t9` runs both and prints them side by side instead of
  arguing it. If they agree, the correction to make is to the *description*
  of the test, not to the number.
- **Quote the block jackknife, not the bootstrap.** 6144 sightlines through
  a 40 Mpc/h box sit ~0.5 Mpc/h apart, far below the correlation length.
  They are not 6144 independent measurements, and every estimator that
  resamples individual lines assumes they are. `t9` reports the effective
  number of independent sightlines. This is a bigger threat to the quoted
  significance than the pairing ever was, and neither of the two original
  options addressed it.
- **The tiling recovers part of what the broken pairing cost.** Line *i* of
  CDM is not line *i* of FCT, but tile *b* of CDM and tile *b* of FCT are
  the same region of the same initial conditions. The large-scale variance
  a pairing was supposed to cancel is a property of the region, not of the
  line, so deleting the same tile from both runs cancels it without any
  sightline being matched.
- **Validate before producing.** The chain from snapshot to P1D has three
  links: SWIFT's LOS output, our regeneration of it, and our tau. Only the
  first is somebody else's code. V1 and V2 above pin the other two against
  an independent implementation, and until they are done a new result only
  adds to the pile of things that would have to be redone. V3 runs in
  parallel because it costs minutes and depends on neither.
- **Lead with results, not code.** Pato is working as a supervisor on this
  repository: he has not read the analysis code and will not until the
  results validate, at which point he will read it properly to check the
  reasoning. Give him verdicts, figures, and what would falsify them.
  Batch the commands so one console dump answers several questions. When a
  stage passes, tell him which files he now has to read and which he can
  keep treating as machinery.
- **The objectives document lives outside this repository, deliberately.**
  Ask Pato where. It is not in `.gitignore` yet — decide that before it is
  ever committed, because a file that reaches a public repo stays in the
  git history even after it is deleted, and forks keep it.

**Non-obvious facts learned this session:**
- **The paired-sightline assumption is false — checked for real on
  2026-08-31.** `stages/02_check_los_match.py --run-a cdm40 --run-b fct40
  --z 3.0` on Clementina: 0% of lines match by index (median offset 20.4
  internal units against a tolerance of 0.006 — off by ~35% of the box),
  and no valid one-to-one permutation exists either (0.07% within
  tolerance). An offset that size is not a shift, it is independent random
  positions: the two runs drew their sightlines separately. Before planning
  any re-shoot, compare the `LineOfSight` block of the two parameter files.
  If it is a different seed, the fix is editing a `.yml`; if SWIFT seeds
  from something run-dependent, it has to go through `legacy/relos.py` with
  an explicit position list, which only works at the redshifts where
  snapshots exist (z = 198, 7, 5, 3, 2).
  `t8_single_los.py` is currently comparing unrelated sightlines — do not
  trust its output, and note that its docstring still states the pairing as
  fact in the header. Fix that docstring whatever else happens, or the next
  session reads it and believes it again.
  **The ensemble results — P1D, the FCT/CDM ratio, the 3D power — never
  assumed pairing and are unaffected; the main finding stands.** If the 20
  sigma is anywhere in circulating talk material, flag it for correction —
  but send one mail with the corrected number after `t9`, not a "we might
  be wrong" mail now.
- **`relos.py`'s own docstring probably explains the pairing failure, and
  it is one command to check.** It says the original `los_*.hdf5` were
  written with `range_when_shooting_down_* = [0, 40]` in a 58.7372 box, so
  each sightline is missing ~32% of its particles AND the rays only sample
  46% of the transverse face. It has a `--uniform N --seed S` flag that
  throws away the old rays and draws new ones over the whole box. If one
  run's file was regenerated with `--uniform` and the other was not, or
  with a different seed, the positions are independent draws — exactly what
  stage 02 measured. Run `tests/t12_relos_roundtrip.py --a FILE` on both
  files and read the ray-position range: if one spans [0, 40] and the other
  [0, 58.74], that is the answer, and it also means the two runs did not
  sample the same sub-volume.
- **Do not validate `relos.py` against the 40 Mpc/h originals.** Those files
  are themselves relos.py output (they were regenerated after the
  truncation bug), so the comparison is circular; and the pre-regeneration
  originals fail by construction because they are truncated. V1 needs a run
  with correct on-the-fly SWIFT LOS output. The 15 new runs are described
  as having exactly that, and `lyman/murgia/*` is the other candidate.
  Check the ray-position range of any candidate first — t12 prints it.
- **The 40 Mpc/h production LOS files ARE truncated. Confirmed by Pato,
  2026-09-01.** `los_0010.hdf5` in both `cdm-box-40-1024` and
  `2-fct-box-40-1024` only holds data out to 40 Mpc instead of the full
  58.7372, from the `h` mix-up. So they cannot serve as the reference for
  V1 - relos.py would be asked to reproduce a file that is itself wrong.
  They are still fine for V2: a code-vs-code comparison feeds both codes
  the same input, and a truncated input is truncated identically for both.
  Say so when handing the file over, so the other side does not read the
  physics as broken when it is the input that is.
  Note what this does NOT settle: if BOTH files are truncated the same way,
  truncation does not explain the pairing failure and a different seed
  does. Step 3 of `scripts/run_validation_A.sh` prints both ranges.
- **Most runs keep their snapshots in a subdirectory, not beside the LOS
  files.** `murgia/cdm` holds `los_000{0..6}.hdf5` next to directories
  `murgia-cdm-lyman_000{0,1,2}/`, each containing that snapshot's virtual
  file and pieces. Not all runs do this. Any glob written by hand will get
  it wrong half the time, so do not write one: step 4 of
  `run_validation_A.sh` reads the real layout and prints the exact,
  ready-to-paste `sbatch` line for every candidate.
- **`relosz.py` is the entry point; `relos.py` is the engine. Keep both.**
  Pato does not call `relos.py` directly - he uses `relosz.py`, and so
  should anything written from here on. They are not alternatives: relosz
  is about forty lines that resolve the LOS file and the snapshot from a
  run directory and a redshift, and then `subprocess` into relos.py, which
  is the ~470 lines that actually walk the snapshot and select particles.
  Deleting relos.py would delete the work. What relosz adds is exactly the
  part that keeps getting written wrong by hand: it matches on the redshift
  each file reports rather than on the index in its name, and it knows that
  snapshots may live in their own subdirectory. Its one limitation is that
  it looks for snapshots as DIRECTORIES, so on a run that stores them as
  plain files beside the LOS files (the 40 Mpc/h pair) it reports "No hay
  snapshot a z=..." and you fall back to relos.py with an explicit
  `--snapshot` glob. `scripts/sbatch_roundtrip.sh` takes a run directory
  and a redshift and goes through relosz.
- **Snapshot and LOS indices do not correspond, in either direction.**
  `los_0003` is not the epoch of `snap_0003`, and the offset is not even
  consistent between runs. Nothing may match on an index. Match on the
  redshift read out of each file, with a tolerance: relos.py refuses a
  snapshot more than dz = 0.02 from the LOS file, so that is the threshold
  worth proposing pairs at. Step 4 of `run_validation_A.sh` does this and
  also prints the near misses, so "this run has no candidate" comes with
  the reason attached.
- **murgia is the clean run. Measured 2026-09-01, not assumed:**
  `murgia/cdm/los_0000.hdf5` holds 1536 sightlines at z = 5.600 in a
  29.52465 box, rays along all three axes, transverse coverage 99.8%, no
  truncation warning, ParticleIDs present, 5779-8504 particles per
  sightline. That is the V1 reference. Two things follow. The
  ParticleIDs make t12's set comparison the strong version rather than the
  coordinate fallback. And murgia uses 1536 sightlines, which is the same
  number the published results quote - worth keeping in mind while
  resolving the 6144-vs-1536 question, though it does not settle it, since
  the published numbers are for the 40 Mpc/h pair and not for murgia.
- **6144 or 1536 sightlines is still open, and it matters more than the
  pairing.** The files hold 6144, the published numbers say 1536. That is a
  factor 2 in any quoted sigma, larger than anything the pairing does.
  Resolve it before requoting an error bar. `t9` aborts if the sightline
  count of the cache disagrees with the LOS file, which is one way to find
  out.
- **`desi_window()` floats with redshift.** It returns
  `k_max = 0.5 pi / R_z`, which grows with z. Integrating an observable
  over it bin by bin means integrating over a different window in each bin,
  and the "redshift evolution" that comes out contains a purely
  instrumental component. `common_window(zs)` returns the fixed
  intersection and is what the new tests use. Enforce it anywhere an
  integrated quantity is computed.
- **The `k << k_max` approximation is not automatically safe at the top of
  the window.** With the fixed common window over z = 2.2 to 4.0
  (`k2 = 0.02548 s/km`), the comoving edge is 1.77 Mpc^-1 at z = 2.2 and
  2.16 Mpc^-1 at z = 4.0. The correction `(k/k_max)^2` is 0.35% for
  `k_F = 30 Mpc^-1`, 1.2% for `k_F = 20`, and 4.7% for `k_F = 10` — and it
  grows with z, which is the axis the whole signal lives on. So the first
  fit in `t11` (constant versus `k_max^2 - k^2` at equal dof) is not a
  formality; it tests a premise that may not hold. This is another reason
  `k_F` has to be measured rather than assumed.
- **`D^2(2.2)/D^2(4.0) = 2.40`.** That is the calculable part of the
  redshift lever, from `units.growth_factor`. Everything above that ratio
  in the measured evolution is filtering scale, thermal state and QLA — the
  residual `t11` reports, and the thing that actually has to be modelled.
- **`b` measured on these runs is `b` under QLA, not `b`.** The QLA scheme
  converts 49.8% of the FCT baryons against 14.1% in CDM, and a run with
  feedback would return much of that gas. That makes the central measurable
  a calibration of this subgrid scheme rather than a portable transfer
  function, and it is the first thing a referee will push on. It does not
  block publication, but it has to be said out loud in the paper — and it
  changes what `t7` is for. `t7` is no longer only "the test that closes the
  argument": it is the only bound available on how much of `b` is QLA.
  `b(QLA)` and `b(stars put back)` are the two ends of the systematic
  interval on the paper's central number.
- **Mass-separability of the converted gas is per-run, not fixed.** The
  original 40 Mpc/h CDM/FCT pair used monofonIC's `masked=2` for the ICs,
  so DM and the QLA-converted gas end up with different particle masses in
  `PartType1` — the "easy" case `stages/00_inspect_snapshot.py --deep`
  detects automatically. The `more_power` batch was generated *without*
  the mask, so DM and baryons share a mass there — the "hard" case that
  needs `PartType4` or `--conv-from-ids`. The 15 new runs were also made
  without the mask, so treat them as the hard case by default; run stage
  00 `--deep` on each one and read the verdict rather than assuming.
- **Clementina has no direct outbound internet.** All SSH traffic,
  including git, has to go through the cluster's HTTP proxy
  (`172.28.3.3:3128`). See section 6 for the exact `~/.ssh/config` entry
  and the GitHub SSH-key setup. This is the flow that actually works on
  Clementina; the plain HTTPS+token clone in Step 2 below is for the
  laptop, not the cluster. It also means the DESI DR1 files cannot be
  fetched from the cluster: download them on the laptop and `scp`.

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
| rescaling to a common `tau_eff` | A = 1 in both runs | discarded — the effect is there before any rescaling, and rescaling makes it slightly deeper |
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
2. **6144 or 1536 sightlines?** The LOS files hold 6144. The published
   numbers say 1536. Find out which was used before requoting any error bar.
3. **Are the sightlines matched between runs?** **Answered 2026-08-31: no.**
   See section 0 for the numbers. `t8_single_los.py` and the 20σ from the
   paired bootstrap are not just provisional now, they are known wrong as
   currently computed — fix by re-shooting the LOS output with a shared
   seed/position list, or re-derive the significance unpaired.
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

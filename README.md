# Lyman-α forest: is the FCT/CDM suppression physical or a bias artefact?

Two simulations, identical initial conditions, one with a standard CDM
primordial spectrum and one with the FCT broken spectrum plus the Poisson
term from PBH discreteness. The 1D flux power spectrum of the Lyman-α
forest is suppressed in FCT by about 19% on large scales, and this
repository is the record of everything we did to try to make that
suppression go away.

It did not go away. What we found instead is that it lives in the gas field
and not in the matter field: the 3D matter power ratio between the two runs
is 1.000, the gas power ratio is 0.808. The quick-Lyman-α star formation
scheme converts 49.8% of the baryons in the FCT box against 14.1% in CDM,
and the gas left behind traces the matter with a lower bias. That is a real
effect on the forest, but its amplitude is set by a subgrid prescription
with no feedback, so it is a constraint on the star-formation treatment as
much as on the cosmology.

Everything below reproduces that chain, one script at a time.

---

## Layout

```
common/      cosmology, the P1D estimator, cache I/O, the run registry
stages/      the pipeline: snapshot -> LOS -> tau -> flux -> P1D -> xi
tests/       the falsification attempts, numbered in the order of the paper
legacy/      the extraction and deposition code the results came from
paper/       the LaTeX source of the document
scripts/     migration, verification, sbatch templates
cache/       tau caches (.npz), not tracked by git
data/        TREECOOL tables and external reference spectra
figures/     output, not tracked by git
```

**The design boundary is the array of `tau`.** Everything upstream of it is
interchangeable — our extractor or spectWizard, either is fine. Everything
downstream is written once, in `common/p1d.py`, and used by every script.
That is what makes it possible to swap the extractor and know that the
estimator did not change underneath.

---

## The data, concretely

Root: `/data/contrib/pad_140/pcolazo/lyman/`

| run | registry name | z=3 snapshot | z=3 LOS file |
|---|---|---|---|
| CDM | `cdm40` | `cdm-box-40-1024/cdm-40-m6-lyman_0003.hdf5` | `cdm-box-40-1024/los_0010.hdf5` |
| FCT | `fct40` | `2-fct-box-40-1024/fct-40-m6-lyman_0003.hdf5` | `2-fct-box-40-1024/los_0010.hdf5` |

- Snapshots exist at z = 198, 7, 5, 3, 2 only. LOS files exist at 17
  redshifts from z = 5 to z = 1.8 in steps of 0.2, which makes a
  redshift-evolution test cheap — the extraction is already done.
- Each LOS file holds **6144** sightlines. Earlier analysis quoted 1536, so
  either `--max-los` was in play or only a subset was processed. Resolve
  this before requoting an error bar.
- Snapshots use SWIFT's `distributed: 1` output: a virtual `.hdf5` of zero
  size next to `.0.hdf5 … .3.hdf5`. Read the virtual file; globbing all of
  them counts every particle twice. `common/runs.py` handles it.
- Both snapshots have `PartType0`, `PartType1` **and `PartType4`**. Whether
  `PartType4` is populated at z=3 decides which route `t7` takes — run
  stage 00 with `--deep`.
- The FCT LOS files are 2.2–2.4 GB against CDM's 3.4–3.8 GB, same number of
  sightlines, same box. That 37% difference in raw bytes is the missing gas,
  measured with `ls`.

## Getting it running on Clementina

```bash
ssh snmgt01
cd /data/contrib/pad_140/pcolazo
git clone <repo-url> lya-repro
cd lya-repro
conda activate astro
pip install -r requirements.txt        # only if something is missing

bash scripts/migrate.sh /data/contrib/pad_140/pcolazo/LOS
bash scripts/verify.sh
```

`migrate.sh` **copies**. It does not delete anything. Once `verify.sh`
prints `VERIFY OK` and you have pushed, you can remove the old `LOS/`
directory by hand — not before.

Two things the migration deliberately keeps:

- `legacy/` holds `swift_extract.py`, `forest_tools.py`, `ionization.py`
  and friends unmodified. The published numbers came out of that code and
  it is not rewritten here, only imported. Editing it means the results
  change; if you edit it, say so in the commit.
- `legacy/prepatch/` holds the `.bak_delta` and `.bak_trho` files. They are
  the only surviving record of what the code looked like before the
  `--delta-max` and `--impose-trho` flags existed.

Long jobs go through `sbatch`, never `srun --pty`. An interactive session
dies with the SSH connection and takes eight hours of extraction with it.

---

## The pipeline

### `stages/00_inspect_snapshot.py` — what is in this file?

**Run it when** you touch a run you have not used before, and **always**
before `t7`. It answers one question: is the gas the star-formation scheme
removed still recoverable, and how do you tell it apart from dark matter?

```bash
python stages/00_inspect_snapshot.py --run fct40 --z 3.0 --deep
python stages/00_inspect_snapshot.py --run NB_1 --z 3.0 --deep
```

| option | what it does |
|---|---|
| `--run NAME` | registry name (`cdm40`, `fct40`, `NB_1`, `poisson_3`, …) or a directory |
| `--z Z` | pick the snapshot closest to this redshift, by reading every header |
| `--snap PATH` | explicit file or index, overrides `--z` |
| `--deep` | also read the `Masses` arrays and cluster them into populations |
| `--max-datasets N` | how many dataset names to print per PartType |

It ends with a verdict block naming the flag you need to pass to `t7`. The
important case is the third one: in the current 40 Mpc/h pair the converted
gas sits in `PartType1` and is separable because it carries the gas particle
mass while the dark matter does not. **In the `more_power` runs the two
species share a mass and that trick dies.** The script says so and tells you
to fall back on matching particle IDs against the initial conditions.

### `stages/01_extract_los.py` — LOS file to `tau`

**Run it when** you have a new LOS file, or when you want to change a
convention and see what it costs.

`swift_extract.py` reads SWIFT **LOS files** (groups `LOS_0000`, `LOS_0001`,
…), not snapshots. It is a library with no `main`, which is why
`python swift_extract.py --help` prints nothing. This stage is the CLI for
it and the only file in the repo that imports from `legacy/`.

```bash
# which LOS files does this run have, and at what redshift?
python stages/01_extract_los.py --run fct40

# geometry only, without reading a particle
python stages/01_extract_los.py --run fct40 --z 3.0 --geometry-only

# production
python stages/01_extract_los.py --run fct40 --z 3.0 --npix 2048 \
    --treecool data/TREECOOL_HM12_G+Q --out cache/cache_fct.npz
```

| option | what it does |
|---|---|
| `--run NAME` / `--los-file PATH` / `--z Z` | resolve the LOS file by redshift, or name it |
| `--npix N` | pixels per sightline |
| `--max-los N` | stop after N lines, for a quick check |
| `--treecool PATH` | TREECOOL table; `Γ_HI` interpolated at the file's redshift |
| `--gamma-hi X` | explicit rate, overrides the table |
| `--hubble-mode file\|fixed\|ref` | `H(z)` from `Cosmology/H`, from `--h0`/`--omega-m`, or from `--reference-run` so both boxes share a *k* axis |
| `--vel-convention swift\|gadget\|none` | peculiar velocities as the file's own conversion gives them, times `sqrt(a)`, or dropped |
| `--delta-max X` | drop gas above this overdensity |
| `--impose-trho T0 GAMMA` | force a T–ρ relation |
| `--no-normalize` | skip the Shepard correction to the SPH deposition |
| `--w-floor X`, `--xh X`, `--he-state S`, `--n-sigma X`, `--exact-voigt`, `--max-b-violation X` | passed straight through to `extract_tau` |
| `--geometry-only` | print `v_box`, `dv`, `k_fund`, `k_nyq`, the conversions and the DESI window, then exit |

The two convention flags exist because they are the two places a factor can
hide. `H` sets the length of the box in km/s and therefore the whole *k*
axis; get it wrong in one run and not the other and you have manufactured a
ratio out of a stretch. The velocity convention is worse: SWIFT and Gadget
differ by `sqrt(a)`, a factor of 2 at z=3. Neither switch edits
`legacy/swift_extract.py` — the Hubble one overrides the metadata object,
the velocity one wraps `_cgs_factor`.

At the end it aggregates the extractor's own diagnostics over all
sightlines: `frac_b_gt_H` (particles outside their own kernel support, must
be ~0), `wsum_raw_med` (the SPH partition of unity before the Shepard
correction; 1.0 is perfect sampling and anything well below it is where the
undersampling of density peaks shows up), and the ray-position source.

**`--geometry-only` is the cheap check.** Run it on any new simulation the
moment the LOS file lands.

### `stages/02_check_los_match.py` — are the two runs sampling the same lines?

**Run this before trusting anything paired.** The sightline positions come
from SWIFT's LOS output configuration, not from anything we choose at
analysis time, so "line *i* of CDM is line *i* of FCT" is an assumption and
it was never checked.

```bash
python stages/02_check_los_match.py --run-a cdm40 --run-b fct40 --z 3.0
```

It reads the ray position of every line from both files using the
extractor's own logic, compares them index by index with periodic wrapping,
and if that fails looks for a permutation that does match — so you can
reorder instead of re-extracting.

| option | what it does |
|---|---|
| `--run-a` / `--run-b` / `--z` | resolve both LOS files by redshift |
| `--a` / `--b` | explicit paths |
| `--tol-frac F` | matched means the offset is under `F` × box (default `1e-4`) |
| `--save-permutation PATH` | write the permutation as `.npy` if one is found |

Two results depend on the answer: `t8_single_los.py`, which compares line
*i* against line *i* pixel by pixel, and the paired bootstrap that produced
the 20σ. The ensemble results — the P1D, its ratio, the 3D power — never
assumed a pairing and are unaffected either way.

### `stages/04_p1d.py` — the main figure

**Run it when** you want the answer.

```bash
python stages/04_p1d.py cache/cache_cdm.npz cache/cache_fct.npz \
    --labels CDM FCT --out figures/p1d

python stages/04_p1d.py cache/cache_cdm.npz cache/cache_fct.npz \
    --labels CDM FCT --norm none --out figures/p1d_raw
```

| option | what it does |
|---|---|
| `--norm taueff\|none\|own` | rescale every run to a common effective optical depth, use the field as the simulation made it, or rescale each to its own Turner+24 value |
| `--tau-eff X` | override the common target (default: Turner+24 at the cache redshift, 0.3719 at z=3) |
| `--kmax-frac F` | cut at `F` × Nyquist. **Leave it at 0.5.** Above half-Nyquist the estimator aliases and the ratio between two runs picks up pure grid structure |
| `--nbins N` | log-*k* bins |
| `--sherwood PATH` | overlay a Sherwood cache as an external reference |
| `--ylim-ratio A B` | y range of the lower panel |

The first cache listed is the denominator of the ratio. Writes a `.png` and
a `.txt` with the ratio tabulated at the DESI anchor scales, which is what
goes into the document.

### `stages/05_xi.py` — the configuration-space control

**Run it when** someone suggests the result is an artefact of the binning.

```bash
python stages/05_xi.py cache/cache_cdm.npz cache/cache_fct.npz \
    --labels CDM FCT --out figures/xi
```

| option | what it does |
|---|---|
| `--norm taueff\|none` | as above |
| `--rmax KMS` | largest separation plotted (default: half the box) |
| `--nboot N` | bootstrap resamples over sightlines for the error band; `0` skips it |

This is not an independent measurement of the field — it is the Fourier
transform of the same estimator. It is an independent check of the *binning
and the interpolation*, which is where artefacts are born. A 19%
suppression below k = 0.01 s/km has to appear on separations above roughly
600 km/s. If it did not, the Fourier result would be a binning artefact.

---

## The tests

Numbered in the order they appear in the document, which is the order of
how cheap they are to kill.

### `tests/t0_rescaling.py` — is it the normalisation?

**This is the test Matteo asked for.** The objection is that CDM and FCT
start from different raw `tau_eff`, we rescale both to the observed value,
and that rescaling is nonlinear in the flux, so it could move power between
scales and manufacture the ratio.

```bash
python tests/t0_rescaling.py --cdm cache/cache_cdm.npz \
    --fct cache/cache_fct.npz --out figures/t0_rescaling
```

Produces six panels: the flux PDF with no rescaling and with it, the P1D
both ways, the ratio computed both ways on one set of axes, and the
difference between the two ratios. If the rescaling were the cause, the raw
ratio would sit at 1 and only the rescaled one would dip.

| option | what it does |
|---|---|
| `--tau-eff X` | common target |
| `--own-taueff` | extra curve with each run at its own `tau_eff` |
| `--kmax-frac F` | half-Nyquist cut |
| `--nbins N` | log-*k* bins for the ratio |

The `.txt` output also carries the saturation statistics — the fraction of
pixels above F = 0.99 and below F = 0.01 in each run and each state. That
number matters on its own: after rescaling, FCT has roughly 21 times more
fully transparent pixels than CDM, which is the opposite of what the first
version of the talk claimed.

### `tests/t7_stars_back.py` — put the converted gas back

**This is the test that closes the argument.** Everything upstream says the
suppression is the difference between the matter field (ratio 1.000) and the
gas field (0.808), and that difference is exactly the baryons QLA removed.
So put them back and re-extract.

The script does **not** compute optical depths. It writes an augmented
**snapshot** whose `PartType0` contains the surviving gas plus the converted
particles, with a smoothing length and a temperature assigned to the latter.

**One step in between.** `swift_extract.py` reads SWIFT LOS files, not
snapshots, so the augmented snapshot has to have sightlines shot through it
before stage 01 can touch it. Two ways:

- re-run SWIFT in LOS-output mode over the augmented snapshot, with the same
  sightline positions as the original (this is the clean route and it is
  what `scripts/sbatch_t7.sh` assumes, expecting a `*_los.hdf5`);
- or use `legacy/relos.py` / `legacy/relosz.py`, which is how the current
  LOS files were regenerated after the truncation bug and already knows how
  to go from a snapshot to LOS groups.

Whichever you pick, the deposition and the Voigt integral stay in
`swift_extract.py`, so `t7` is validated by the same extractor as everything
else and there is no second code path to trust.

```bash
# always the dry run first
python tests/t7_stars_back.py --run fct40 --z 3.0 --dry-run \
    --conv-parttype 1 --conv-mass-max 2.0e-3

# then for real, with the T0 sweep that shows the answer does not
# depend on the temperature we invented for the reinjected gas
python tests/t7_stars_back.py --run fct40 --z 3.0 \
    --conv-parttype 1 --conv-mass-max 2.0e-3 \
    --t0-sweep 0.5 1.0 2.0 --out data/aug_fct.hdf5

python stages/01_extract_los.py --snap data/aug_fct_T0x1.hdf5 \
    --nlos 1536 --npix 2048 --seed 42 --out cache/cache_fct_t7.npz
python stages/04_p1d.py cache/cache_cdm.npz cache/cache_fct.npz \
    cache/cache_fct_t7.npz --labels CDM FCT "FCT + converted gas" \
    --out figures/t7
```

| option | what it does |
|---|---|
| `--conv-parttype N` | where the converted particles live, 1 or 4 |
| `--conv-mass-max X` | internal units; particles lighter than this in `--conv-parttype` are the converted gas. Get the value from stage 00 `--deep` |
| `--conv-from-ids PATH` | the fallback for the equal-mass runs: anything whose ID was `PartType0` in the ICs and is not gas now was converted |
| `--h-mode grid\|knn` | smoothing length from a CIC density field (fast, O(N)) or from the distance to the *N*th neighbour (slower, needs a tree) |
| `--ngrid N` / `--nngb N` | resolution of each |
| `--t-mode trho\|fixed` | put the reinjected gas on the T–ρ relation at its own overdensity, or at a fixed temperature |
| `--t0 K` / `--gamma G` | the T–ρ relation |
| `--t0-sweep A B C` | write one augmented snapshot per factor on T0. **Use it.** The temperature of the reinjected gas is the one thing in this test we invented, and the sweep is what shows the answer does not turn on it |
| `--frac F` | reinject only a random fraction, so you get a curve instead of a point |
| `--dry-run` | count everything, write nothing |

The `[.]` lines in the output mark datasets filled with the gas median
because there was nothing sensible to copy. Read them; if the extractor
depends on one of those fields, the fill matters.

### `tests/t8_single_los.py` — look at one sightline

**Run it when** you want to see the mechanism instead of averaging over it.
The two boxes have identical initial conditions and the sightlines are drawn
from the same seed, so line *i* passes through the same structure in both.

```bash
python tests/t8_single_los.py --cdm cache/cache_cdm.npz \
    --fct cache/cache_fct.npz --los 12 340 900 --out figures/t8

python tests/t8_single_los.py --cdm cache/cache_cdm.npz \
    --fct cache/cache_fct.npz --pick extreme --n 3 --out figures/t8_extreme
```

Three outputs: the `tau` and `F` profiles along each line with the strong
CDM absorbers shaded, the P1D of each individual line against the ensemble,
and a peak census. The census is the mechanism in one table — for every
pixel where CDM has `tau > 2`, what does FCT have at that same pixel? A
mean ratio well below 1 means the absorber is in the same place in both
runs but shallower in FCT, because the gas that made it has been converted.

| option | what it does |
|---|---|
| `--los I J K` | explicit indices |
| `--pick median\|extreme\|strongest\|random` | pick automatically: lines whose own ratio is closest to the ensemble (the typical case), the most suppressed, the deepest absorbers, or a seeded random draw |
| `--n N` | how many, when using `--pick` |
| `--no-rescale` | leave A = 1 in both |
| `--peak-tau X` | what counts as a strong CDM absorber |
| `--smooth-kms X` | boxcar smoothing, display only |

### `tests/test_estimator.py` — arithmetic

```bash
python tests/test_estimator.py
```

Unit tests on `common/p1d.py` with analytic inputs: a cosine must land in
one bin with amplitude `L a²/4`, `xi(0)` must equal the variance of the
flux contrast, `solve_A` must invert `tau_eff` to machine precision, a
constant field must have no power. No simulation output is synthesised
anywhere; every scientific number in this repository comes from a real
snapshot.

---

## Adding a new simulation

The `more_power` batch (`NB_1`…`NB_6`, `poisson_1`…`poisson_6`) is already
in the registry in `common/runs.py`. Anything else, add a line there or
pass a directory path to `--run`.

Snapshot indices are **not** comparable between runs — they come from
different output lists, so `los_0010` is z=3 in one box and `_0003` is z=3
in another. Nothing in this repository uses an index. `--z 3.0` opens every
snapshot in the directory, reads `Header/Redshift`, takes the closest, and
refuses if the closest is more than 0.05 away. That is what makes a run
that did not exist when this was written work without editing code.

When a new run lands:

```bash
python stages/01_extract_los.py --run NB_1 --z 3.0 --geometry-only
python stages/00_inspect_snapshot.py --run NB_1 --z 3.0 --deep
python stages/01_extract_los.py --run NB_1 --z 3.0 --seed 42 \
    --out cache/cache_NB1.npz
python stages/04_p1d.py cache/cache_cdm.npz cache/cache_NB1.npz \
    --labels CDM NB_1 --out figures/p1d_NB1
```

The first two commands cost seconds and catch the two mistakes that
actually happen: the wrong redshift, and converted gas you cannot separate.

---

## Provenance

Every cache written by this repository carries a `prov` field with the git
hash (marked `-dirty` if the tree was not clean), the full command line,
the hostname, the timestamp, and a short hash of the input snapshot. Read
it with:

```python
import numpy as np, json
print(json.loads(str(np.load("cache/cache_fct.npz")["prov"])))
```

Caches made before this repository existed only have `tau`, `dv` and `z`.
They still load; the missing fields come back as `unknown`.

---

## What we already know, so you can skip re-deriving it

Numbers from the audit, at z = 3, box 40 Mpc/h, 1536 × 2048:

| quantity | CDM | FCT |
|---|---|---|
| gas particles in the ICs | 134,217,728 | 134,217,728 |
| gas particles at z = 3 | 115,318,346 | 67,341,600 |
| converted by QLA | 14.1% | 49.8% |
| `tau_eff` rescaling A | 0.78961 | 0.87198 |

| test | result | verdict |
|---|---|---|
| rescale both to a common `tau_eff` | ratio moves the wrong way | not the cause |
| sampling variance, paired bootstrap and split-half | 20σ | not the cause |
| deposition scheme, no normalisation | survives | not the cause |
| remove dense gas from CDM, Δ > 100 and Δ > 30 | < 1% | not the cause |
| impose the same T–ρ relation, γ = 1.40 | ~2% | not the cause |
| 3D power, matter | 1.000 | — |
| 3D power, baryons (gas + converted) | 1.000 | — |
| **3D power, gas only** | **0.808** | **this is it** |
| random removal of 40% of the gas anywhere | 0.9999 | mass budget is not the cause |
| cut CDM by density to FCT's gas mass | 0.10 | it is *where* the gas is removed |

The last two rows together are the argument. Removing 40% of the gas at
random costs nothing. Removing the same mass by cutting the filaments costs
90% of the large-scale power. FCT loses only 19%, so its removal is neither
random nor a clean filament cut — it is concentrated in peaks, and what
survives is a less biased tracer of the same matter field.

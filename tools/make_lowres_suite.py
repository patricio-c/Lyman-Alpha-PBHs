#!/usr/bin/env python3
"""
Generate the ten low-resolution SWIFT parameter files and run scripts.

Ten hand-edited copies of a yml is where mistakes get in, and the first one
written by hand already had four of them. This builds all ten from one
template plus a table, so the only thing that can be wrong is the template.

What the suite is
-----------------
Four sets of initial conditions, ten runs. Legs (a) and (c) have the SAME
particle mass, which is what makes the splicing correction valid; leg (b) is
the existing production run at eight times better mass resolution in the same
box as (c).

    (a)  80 Mpc/h, monofonIC grid 1024 -> gas 512^3
    (c)  40 Mpc/h, monofonIC grid  512 -> gas 256^3

    leg  model  H_reion_z  over_density   purpose
    a    CDM      7.5        1000         large-box leg of the splice
    a    FCT      7.5        1000
    c    CDM      6          1000         reionization derivative
    c    FCT      6          1000
    c    CDM      7.5        1000         small-box leg AND central z_reion
    c    FCT      7.5        1000
    c    CDM      9          1000         reionization derivative
    c    FCT      9          1000
    c    CDM      7.5        3000         threshold sensitivity
    c    FCT      7.5        3000

The threshold variation runs in BOTH models on purpose. Varying it only in
the high-power run would confound the ratio, which is the observable.

What is derived rather than copied
----------------------------------
**The gravitational softening.** It must scale with the mean interparticle
separation, and copying the production value into a run with eight times the
particle mass puts a force-resolution change inside the very measurement that
is supposed to isolate a mass-resolution change. Both new legs sit at exactly
2x the production spacing:

    production  box 40, gas 512^3  ->  0.114721 Mpc
    leg (a)     box 80, gas 512^3  ->  0.229442 Mpc   (2.0000x)
    leg (c)     box 40, gas 256^3  ->  0.229442 Mpc   (2.0000x)

so 0.0018 -> 0.0036 and 0.0007 -> 0.0014. Computed here from box and grid so
that it cannot drift if the configuration changes.

**The power-spectrum grid.** Leg (c) has gas on 256^3, so a 512 grid reports
pure shot noise above its own sampling limit. Set to the gas grid.

**The mesh and the node count**, both scaled to the particle count.

What is preserved verbatim from the production file
---------------------------------------------------
Cosmology, unit system, time integration, SPH, entropy floor, cooling tables,
the three output lists, `Neutrino: use_model_none`, and `a_end = 0.4`. That
last one is not a typo: the LOS and power-spectrum lists run down to z = 1.8,
so the run has to continue past it.

The line-of-sight block carries NO `*_range_*` lines. Those are what broke the
production pair: `[0, 40]` in Mpc/h written where internal Mpc was expected,
so the rays sampled 46% of the transverse face.

Usage
-----
    python tools/make_lowres_suite.py --check      # report, write nothing
    python tools/make_lowres_suite.py              # write the tree
    python tools/make_lowres_suite.py --force      # overwrite existing files
"""

from __future__ import annotations

import argparse
import glob
import os

H = 0.681

# --- paths, all absolute because SWIFT is run from the job directory -------
ROOT = "/data/contrib/pad_140/pcolazo/lyman/low_resolution"
IC_ROOT = "/data/contrib/pad_140/pcolazo/IC_lyman_alpha/low_resolution"
SWIFT_BIN = "/data/contrib/pad_140/pcolazo/lyman/cdm-box-40-1024/swift_mpi"
COOLING = "/data/contrib/pad_140/pcolazo/lyman/UV_dust1_CR1_G1_shield1.hdf5"
OUT_SNAP = f"{ROOT}/output.txt"
OUT_LOS = f"{ROOT}/los_output_list.txt"
OUT_PS = f"{ROOT}/output_list_ps.txt"

# --- the production reference the softening is scaled from ----------------
PROD_BOX, PROD_GAS = 40.0, 512
PROD_SOFT, PROD_SOFT_MAX = 0.0018, 0.0007

# --- the two legs ---------------------------------------------------------
LEGS = {
    #      box [Mpc/h]  monofonIC grid  gas side  nodes  hours
    "a": dict(box=80.0, grid=1024, gas=512, nodes=4, hours=48),
    "c": dict(box=40.0, grid=512, gas=256, nodes=1, hours=12),
}

MODELS = {"CDM": "cdm", "FCT": "fct"}

# (leg, model, H_reion_z, over_density, subdirectory under ROOT)
RUNS = [
    ("a", "CDM", 7.5, 1000, "a/CDM"),
    ("a", "FCT", 7.5, 1000, "a/FCT"),
    ("c", "CDM", 6.0, 1000, "c/z_reio/6/CDM"),
    ("c", "FCT", 6.0, 1000, "c/z_reio/6/FCT"),
    ("c", "CDM", 7.5, 1000, "c/z_reio/7.5/CDM"),
    ("c", "FCT", 7.5, 1000, "c/z_reio/7.5/FCT"),
    ("c", "CDM", 9.0, 1000, "c/z_reio/9/CDM"),
    ("c", "FCT", 9.0, 1000, "c/z_reio/9/FCT"),
    ("c", "CDM", 7.5, 3000, "c/delta_dif/CDM"),
    ("c", "FCT", 7.5, 3000, "c/delta_dif/FCT"),
]


def spacing(box_hmpc, gas_side):
    """Mean comoving separation between gas particles, in internal Mpc."""
    return (box_hmpc / H) / gas_side


def softening(box_hmpc, gas_side):
    """
    Softening scaled from the production run by the interparticle spacing.

    Derived rather than tabulated: the whole point of these runs is a
    controlled change in mass resolution, and a softening that does not
    follow it silently adds a second change.
    """
    r = spacing(box_hmpc, gas_side) / spacing(PROD_BOX, PROD_GAS)
    return PROD_SOFT * r, PROD_SOFT_MAX * r, r


def ic_path(leg, model):
    """
    (path, how) for this leg's initial conditions.

    The filename encodes monofonIC's GridRes, NOT the gas particle count.
    `masked` with ParticleMaskType 2 puts gas on (GridRes/2)^3, so
    IC_N1024B080 holds 512^3 gas particles. Both facts are true and an
    earlier version of this function conflated them, encoding the particle
    count and then looking for files that do not exist.

    If the constructed name is absent but exactly one IC for this model
    sits in the directory, that file is used and the substitution is
    reported. A rename upstream should not silently produce ten parameter
    files pointing at nothing.
    """
    tag = MODELS[model]
    d = f"{IC_ROOT}/{model}/{leg}"
    name = (f"IC_N{LEGS[leg]['grid']:04d}"
            f"B{int(LEGS[leg]['box']):03d}_200_{tag}.hdf5")
    path = os.path.join(d, name)
    if os.path.exists(path):
        return path, "ok"
    hits = sorted(glob.glob(os.path.join(d, f"IC_*_{tag}.hdf5")))
    if len(hits) == 1:
        return hits[0], "glob"
    if len(hits) > 1:
        return path, f"AMBIGUOUS ({len(hits)} candidates)"
    return path, "MISSING"


def yml(leg, model, zreion, overdens, name, ic):
    L = LEGS[leg]
    soft, soft_max, _ = softening(L["box"], L["gas"])
    mesh = L["gas"]
    ps_grid = L["gas"]
    return f"""# Generated by tools/make_lowres_suite.py -- do not hand-edit.
# Regenerate instead, so that the ten files cannot drift apart.

MetaData:
  run_name:   {name}

InternalUnitSystem:
  UnitMass_in_cgs:     1.98841e43    # 10^10 M_sun in grams
  UnitLength_in_cgs:   3.08567758e24 # Mpc in centimeters
  UnitVelocity_in_cgs: 1e5           # km/s in centimeters per second
  UnitCurrent_in_cgs:  1             # Amperes
  UnitTemp_in_cgs:     1             # Kelvin

Cosmology:
  Omega_cdm:      0.2560109
  Omega_lambda:   0.693922
  Omega_b:        0.0486000000000
  Omega_r:        7.79165e-05
  h:              {H}
  a_begin:        0.004975124   # z_start = 200
  a_end:          0.4           # z_end = 1.5; the LOS/PS lists reach z=1.8

TimeIntegration:
  dt_min:     1e-10
  dt_max:     1e-2
  time_begin: 0
  time_end:   10

Snapshots:
  distributed:         1
  basename:            {name}
  output_list_on:      1
  output_list:         {OUT_SNAP}

Statistics:
  scale_factor_first:  0.016
  delta_time:          1.01

# Softening scaled from the production run by the interparticle spacing:
# box {L['box']:.0f} Mpc/h with gas {L['gas']}^3 gives {spacing(L['box'], L['gas']):.6f} Mpc
# against {spacing(PROD_BOX, PROD_GAS):.6f} Mpc in the 40/512^3 production run.
Gravity:
  eta:                         0.025
  MAC:                         adaptive
  epsilon_fmm:                 0.001
  theta_cr:                    0.7
  use_tree_below_softening:    1
  mesh_side_length:            {mesh}
  comoving_DM_softening:         {soft:.6f}
  max_physical_DM_softening:     {soft_max:.6f}
  comoving_baryon_softening:     {soft:.6f}
  max_physical_baryon_softening: {soft_max:.6f}

SPH:
  resolution_eta:        1.2348
  h_min_ratio:           1e-8
  h_max:                 0.5
  CFL_condition:         0.2
  minimal_temperature:   10
  initial_temperature:   268.7
  H_mass_fraction:       0.756

Scheduler:
  max_top_level_cells:   16
  cell_split_size:       200

InitialConditions:
  file_name: {ic}
  periodic:   1
  cleanup_h_factors: 0
  cleanup_velocity_factors: 0

# No *_range_* lines. Those are what broke the production pair: [0, 40] in
# Mpc/h written where internal Mpc was expected, so the rays covered 46% of
# the transverse face. Omitting them makes SWIFT use the whole box.
LineOfSight:
  basename:                   los
  output_list_on:             1
  output_list:                {OUT_LOS}
  num_along_x:                2048
  num_along_y:                2048
  num_along_z:                2048

QLACooling:
  dir_name:                {COOLING}
  H_reion_z:               {zreion}
  H_reion_eV_p_H:          2.0
  He_reion_z_centre:       3.0
  He_reion_z_sigma:        0.5
  He_reion_eV_p_H:         2.0
  rapid_cooling_threshold: 0.333333

QLAStarFormation:
  over_density:              {overdens}

QLAEntropyFloor:
  density_threshold_H_p_cm3: 0.1
  over_density_threshold:    10.
  temperature_norm_K:        8000

Restarts:
  enable:             1
  onexit:             0
  delta_hours:        3.0
  max_run_time:       {LEGS[leg]['hours'] - 1}.0
  resubmit_on_exit:   1
  resubmit_command:   sbatch ./run.sh 1

# grid_side_length is the gas grid: a finer grid reports shot noise above the
# sampling limit, which is how the cube-corner modes got into the analysis
# the first time round.
PowerSpectrum:
  grid_side_length:  {ps_grid}
  num_folds:         6
  requested_spectra: ["matter-matter", "gas-gas", "cdm-cdm", "gas-matter"]
  output_list_on:    1
  output_list:       {OUT_PS}

Neutrino:
  use_model_none: 1
"""


def runsh(leg, model, name, workdir, ymlname):
    L = LEGS[leg]
    return f"""#!/bin/bash
# Generated by tools/make_lowres_suite.py -- do not hand-edit.
#SBATCH --account=pad_140
#SBATCH --job-name={name}
#SBATCH --output=aa_swift_%j.out
#SBATCH --error=aa_swift_%j.err
#SBATCH --nodes={L['nodes']}
#SBATCH --ntasks={L['nodes']}
#SBATCH --cpus-per-task=64
#SBATCH --partition=cpunode
#SBATCH --time={L['hours']}:00:00

module purge
module load intel/2024.0.0
module load impi/2021.11
module load fftw/3.3.10
module load phdf5/1.14.0
module load metis/5.1.0
module load gsl/2.8

export I_MPI_FABRICS=shm:ofi
export I_MPI_OFI_PROVIDER=psm3
export FI_PROVIDER=psm3
export I_MPI_PIN=off
export I_MPI_EXTRA_FILESYSTEM=1
export I_MPI_EXTRA_FILESYSTEM_FORCE=gpfs

SWIFT_BIN="{SWIFT_BIN}"
cd {workdir} || exit 1
echo "romio_cb_write disable" > {workdir}/romio_hints
export ROMIO_HINTS={workdir}/romio_hints

# The production script computed this flag and then never used it, because
# its mpiexec line carried a hardcoded -r. Every run therefore started in
# restart mode, including the first one from initial conditions.
if [ "$1" == "1" ]; then
    RESTART_FLAG="-r"
    echo ">>> RESTART"
else
    RESTART_FLAG=""
    echo ">>> clean start from ICs"
fi

mpiexec -np $SLURM_NTASKS "$SWIFT_BIN" $RESTART_FLAG -v 1 --cosmology \\
        --line-of-sight --self-gravity --power --quick-lyman-alpha \\
        --threads=$SLURM_CPUS_PER_TASK {ymlname}
"""


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", default=ROOT)
    ap.add_argument("--check", action="store_true",
                    help="report and write nothing")
    ap.add_argument("--force", action="store_true",
                    help="overwrite files that already exist")
    args = ap.parse_args()

    print("softening, derived from the interparticle spacing")
    s0 = spacing(PROD_BOX, PROD_GAS)
    print(f"  production  box {PROD_BOX:.0f}, gas {PROD_GAS}^3 -> "
          f"{s0:.6f} Mpc   softening {PROD_SOFT} / {PROD_SOFT_MAX}")
    for leg, L in LEGS.items():
        soft, soft_max, r = softening(L["box"], L["gas"])
        print(f"  leg {leg}       box {L['box']:.0f}, gas {L['gas']}^3 -> "
              f"{spacing(L['box'], L['gas']):.6f} Mpc   ratio {r:.4f}   "
              f"softening {soft:.6f} / {soft_max:.6f}")
    print()

    missing, written, skipped = [], [], []
    print(f"{'run':<26s} {'leg':>3s} {'z_reion':>8s} {'Delta':>6s} "
          f"{'nodes':>5s}  {'IC':<8s} file")
    for leg, model, zre, od, sub in RUNS:
        tag = MODELS[model]
        name = f"{tag}-{int(LEGS[leg]['box'])}-lyman-{leg}-z{zre}-d{od}"
        workdir = os.path.join(args.root, sub)
        ic, how = ic_path(leg, model)
        if how not in ("ok", "glob"):
            missing.append(f"{ic}   [{how}]")
        print(f"{name:<26s} {leg:>3s} {zre:>8} {od:>6} "
              f"{LEGS[leg]['nodes']:>5d}  {how:<8s} {os.path.basename(ic)}")

        if args.check:
            continue
        os.makedirs(workdir, exist_ok=True)
        ymlname = f"{tag}-{leg}.yml"
        for fname, body in ((ymlname, yml(leg, model, zre, od, name, ic)),
                            ("run.sh", runsh(leg, model, name, workdir,
                                             ymlname))):
            path = os.path.join(workdir, fname)
            if os.path.exists(path) and not args.force:
                skipped.append(path)
                continue
            with open(path, "w") as fh:
                fh.write(body)
            if fname == "run.sh":
                os.chmod(path, 0o755)
            written.append(path)

    print()
    if missing:
        print(f"{len(missing)} initial condition file(s) unusable:")
        for m in sorted(set(missing)):
            print("  ", m)
        print("   generate them before submitting; the parameter files point "
              "at these paths.")
    if skipped:
        print(f"{len(skipped)} file(s) already existed and were left alone "
              f"(use --force to overwrite):")
        for s in skipped:
            print("  ", s)
    if written:
        print(f"{len(written)} file(s) written.")
    if args.check:
        print("--check: nothing was written.")


if __name__ == "__main__":
    main()

# %% [markdown]
# # Murgia+2019 vs SWIFT at z = 5
#
# Percent-cell format: VS Code and Jupyter open this as a notebook, and it
# stays a plain-text diff in git.  Three pieces:
#
# 1. **Their curves** come out of the vector PDF of the arXiv e-print, via
#    `scripts/extract_murgia_fig1.py`.  Not digitised by eye.
# 2. **Our P1D** is what `stages/04_p1d.py` writes.  This file reads the
#    `.txt` it leaves next to each figure; when it cannot find one it falls
#    back to values transcribed from the console and says so out loud.
# 3. **Our non-linear matter power** is the Pylians `p_matter_*z5.txt`.
#
# The one thing worth changing is the normalisation.  See the last cell.

# %%
import os
import subprocess
import sys

import numpy as np
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath("__file__")))
WORK = os.path.join(ROOT, "data", "murgia_fig1")   # pdf / svg / json
PK = os.path.join(ROOT, "data")                    # p_matter_*z5.txt
FIGS = os.path.join(ROOT, "figures")               # stage 04 outputs

BOX_MPC, BOX_HMPC = 29.52465, 20.0                 # the murgia box
H = BOX_HMPC / BOX_MPC                             # 0.677400
VBOX = 2738.925                                    # km/s at z=5, from the LOS header
KFAC = VBOX / BOX_HMPC                             # k[s/km] -> k[h/Mpc]
print(f"h = {H:.6f}    k[s/km] -> h/Mpc  x {KFAC:.3f}")

# %% [markdown]
# ## 1. The published curves

# %%
import json

js = os.path.join(WORK, "murgia2019_fig1.json")
if not os.path.exists(js):
    subprocess.run([sys.executable, os.path.join(ROOT, "scripts", "extract_murgia_fig1.py"),
                    "--workdir", WORK], check=True)
MUR = {n: (np.array(d["k_hMpc"]), np.array(d["pct"])) for n, d in json.load(open(js)).items()}


def at(name, q):
    k, v = MUR[name]
    return np.interp(np.log10(q), np.log10(k), v)


print("curves:", ", ".join(sorted(MUR)))
print("linear M3/M2 (must be 10.00):",
      " ".join(f"{at('linear_M3', q) / at('linear_M2', q):.2f}" for q in (5, 8.22, 10, 15, 19)))

# %% [markdown]
# ## 2. Our P1D
#
# `stages/04_p1d.py` prints a table of ratios at the end of its `.txt`.  Three
# normalisations matter:
#
# - `none`   : `A = 1` everywhere, the field as the simulation made it;
# - `turner` : rescaled to `units.tau_eff_turner24(z)` -- **wrong at z = 5**,
#              that function is extrapolating outside its calibration range;
# - `cdmref` : rescaled to the reference run's own raw `tau_eff`, 2.84299.
#              This is the one to use.

# %%
ANCHORS = np.array([0.0030, 0.0050, 0.0100, 0.0200, 0.0300, 0.0600])

# transcribed from the Clementina console, 2026-09-02 and 2026-09-03
FALLBACK = {
    "none":   {"M2": [1.0480, 1.0595, 1.0685, 1.0857, 1.1117, 1.2174],
               "M3": [1.0991, 1.1061, 1.2196, 1.3224, 1.5384, 2.2010]},
    "turner": {"M2": [1.0189, 1.0239, 1.0135, 1.0296, 1.0325, 1.1090],
               "M3": [0.9653, 0.9427, 1.0497, 1.0969, 1.1877, 1.4325]},
    "cdmref": {"M2": [1.0125, 1.0224, 1.0282, 1.0423, 1.0664, 1.1670],
               "M3": [0.9878, 0.9882, 1.0851, 1.1618, 1.3436, 1.8926]},
}
FILES = {"none": "p1d_murgia_z5_raw.txt", "cdmref": "p1d_murgia_z5_cdmref.txt"}


def read_stage04(path):
    """Return {label: [ratio at each anchor]} from the table stage 04 writes."""
    lines = open(path).read().splitlines()
    head = next(i for i, l in enumerate(lines) if l.strip().startswith("k [s/km]"))
    labels = lines[head].split()[2:]
    rows = {}
    for l in lines[head + 1:]:
        f = l.split()
        if len(f) != 1 + len(labels):
            break
        rows[float(f[0])] = [float(x) for x in f[1:]]
    ks = sorted(rows)
    return {lab: [rows[k][j] for k in ks] for j, lab in enumerate(labels)}


OURS = {}
for key, fn in FILES.items():
    p = os.path.join(FIGS, fn)
    if os.path.exists(p):
        OURS[key] = read_stage04(p)
        print(f"{key:7s} read from {p}")
    else:
        OURS[key] = FALLBACK[key]
        print(f"{key:7s} FALLBACK to transcribed values ({p} not found)")
OURS.setdefault("turner", FALLBACK["turner"])

# %% [markdown]
# ## 3. The table that decides

# %%
kh = ANCHORS * KFAC
for tag in ("M3", "M2"):
    print(f"\n--- {tag}: excess over CDM, per cent, z = 5.0 ---")
    print(f"{'k[s/km]':>8} {'k[h/Mpc]':>9} | {'none':>8} {'Turner':>8} "
          f"{'ref CDM':>8} {'Murgia':>8} | {'refCDM/Mur':>10}")
    for i, (ks, kq) in enumerate(zip(ANCHORS, kh)):
        m = at(f"flux_{tag}", kq)
        a, b, c = (100 * (OURS[k][tag][i] - 1) for k in ("none", "turner", "cdmref"))
        r = f"{c / m:10.2f}" if abs(m) > 1.0 else f"{'-':>10}"
        print(f"{ks:8.4f} {kq:9.3f} | {a:8.2f} {b:8.2f} {c:8.2f} {m:8.2f} | {r}")

# %% [markdown]
# ## 4. Non-linear matter -- the test of the initial conditions
#
# This goes through neither the extractor nor the normalisation.  If it
# agrees, the ICs and the non-linear evolution are fine and any residual in
# the flux panel belongs to the gas physics or to the estimator.

# %%
def read_pk(name):
    a = np.loadtxt(os.path.join(PK, name))
    return a[:, 1], a[:, 2]            # k [Mpc^-1], P [Mpc^3], shot noise subtracted


kc, Pc = read_pk("p_matter_cdmz5.txt")
_, P2 = read_pk("p_matter_m2z5.txt")
_, P3 = read_pk("p_matter_m3z5.txt")


def mat(P, q):
    return 100 * (np.interp(q * H, kc, P) / np.interp(q * H, kc, Pc) - 1)


print(f"{'k[h/Mpc]':>9} | {'M2 ours':>9} {'M2 Murgia':>10} | {'M3 ours':>9} {'M3 Murgia':>10}")
for q in (4.108, 5, 8.217, 10, 15, 19):
    print(f"{q:9.2f} | {mat(P2, q):9.2f} {at('nonlinear_M2', q):10.2f} |"
          f" {mat(P3, q):9.2f} {at('nonlinear_M3', q):10.2f}")

# %% [markdown]
# ## 5. The figure

# %%
from matplotlib.lines import Line2D

BAND = (0.6834, 10.7462)               # MIKE/HIRES window, measured off the same PDF
C2, C3 = "#1f4fbf", "#c0271f"

fig, (ax, bx) = plt.subplots(2, 1, figsize=(9.2, 8.6), sharex=True,
                             gridspec_kw={"height_ratios": [1.55, 1], "hspace": 0.08})
for a in (ax, bx):
    a.axvspan(*BAND, color="0.92", lw=0, zorder=0)
    a.axhline(0, color="0.55", lw=1.1)
    a.set_xscale("log")
    a.set_xlim(0.30, 22)
    a.grid(alpha=0.25, which="both")

q = np.logspace(np.log10(0.30), np.log10(22), 400)
for tag, col, P in (("M2", C2, P2), ("M3", C3, P3)):
    k, v = MUR[f"flux_{tag}"]
    m = (k >= 0.30) & (k <= 22)
    ax.plot(k[m], v[m], color=col, lw=2.4)
    ax.plot(kh, [100 * (x - 1) for x in OURS["cdmref"][tag]], "o--",
            color=col, lw=1.7, ms=8, mec="w", mew=1.6)
    k, v = MUR[f"nonlinear_{tag}"]
    m = (k >= 0.30) & (k <= 22)
    bx.plot(k[m], v[m], color=col, lw=2.2)
    bx.plot(q, mat(P, q), "--", color=col, lw=1.9)

ax.legend(handles=[
    Line2D([], [], color="k", lw=2.4, label="Murgia+2019, published (GADGET-III)"),
    Line2D([], [], color="k", lw=1.7, ls="--", marker="o", ms=8, mec="w",
           label="this work (SWIFT + QLA)"),
    Line2D([], [], color=C2, lw=0, marker="o", ms=9, label=r"$M_{\rm PBH}=10^2\,M_\odot$"),
    Line2D([], [], color=C3, lw=0, marker="o", ms=9,
           label=r"$10^3\,M_\odot$, $f_{\rm PBH}=1$")],
    frameon=True, framealpha=0.95, loc="upper left", fontsize=10)
ax.set_ylim(-16, 132)
ax.set_ylabel(r"1D flux power excess over $\Lambda$CDM  [%]")
bx.set_ylim(-0.4, 7.6)
bx.set_ylabel("non-linear matter power excess  [%]")
bx.set_xlabel(r"$k$  [ $h\,$Mpc$^{-1}$, comoving ]")
bx.text(0.02, 0.92, "solid: Murgia+2019   dashed: this work (Pylians, z=5 snapshots)",
        transform=bx.transAxes, fontsize=9.5, color="0.35")
ax.set_title("SWIFT reproduces Murgia+2019 at z = 5 - independent code, same models",
             fontsize=12.5, loc="left")

os.makedirs(FIGS, exist_ok=True)
fig.savefig(os.path.join(FIGS, "murgia_vs_swift_z5.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(FIGS, "murgia_vs_swift_z5.png"), dpi=200, bbox_inches="tight")
plt.show()

# %% [markdown]
# ## 6. To play with
#
# The only thing worth varying is the target the runs are normalised to.  On
# the cluster:
#
# ```bash
# for T in 1.5 2.0 2.5 2.84299 3.2; do
#   python stages/04_p1d.py cache/cache_murgia_cdm_z5.0.npz cache/cache_murgia_M3_z5.0.npz \
#       --labels cdm M3 --norm taueff --tau-eff $T --out figures/sweep_$T
# done
# ```
#
# then point `FILES` at those `.txt` and re-run cells 2 and 3.  That measures
# directly how much of `R(k)` is normalisation and how much is physics, which
# is the question still open.

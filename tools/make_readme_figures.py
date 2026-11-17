#!/usr/bin/env python3
"""Build the explanatory figures used in README.md from the bundled examples.

usage (from the repository root):  python tools/make_readme_figures.py
writes docs/figures/{concept,model_selection,results_overview}.png and copies
each example's published subevent figure to docs/figures/<name>_subevents.png
"""
import os
import re
import shutil
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

REPO = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
EX = os.path.join(REPO, "examples")
OUT = os.path.join(REPO, "docs", "figures")
os.makedirs(OUT, exist_ok=True)

EVENTS = [  # name, label
    ("chile", "2024 Chile–Argentina  Mw 7.4  (123 km deep)"),
    ("california", "2024 offshore N. California  Mw 7.0"),
    ("venezuela", "2026 coastal Venezuela  Mw 7.2"),
    ("kamchatka", "2025 Kamchatka megathrust  Mw 8.8"),
]
SIG = 2 * np.sqrt(2 * np.log(10))   # duration = full width at 10% of the peak

plt.rcParams.update({"font.size": 11, "axes.spines.top": False, "axes.spines.right": False})


def load(name):
    r = os.path.join(EX, name, "reference")
    model = np.loadtxt(os.path.join(r, "Input.model"), ndmin=2)   # t x y T vr theta z
    fm = np.loadtxt(os.path.join(r, "fm.dat"), ndmin=2)[:, 1:7]
    m0 = np.array([np.sqrt((m[0]**2 + 2*m[1]**2 + 2*m[2]**2 + m[3]**2 + 2*m[4]**2 + m[5]**2) / 2)
                   for m in fm])                                     # x 1e27 dyne-cm
    return model, m0


def pulse(t, t0, dur, m0):
    s = dur / SIG
    return m0 / (np.sqrt(2 * np.pi) * s) * np.exp(-(t - t0) ** 2 / (2 * s * s))


def lcurve(name):
    vals = {}
    for line in open(os.path.join(EX, name, "reference", "lcurve.txt")):
        m = re.match(r"\s*(\d+)\s+([0-9.]+)\s", line)
        if m:
            vals[int(m.group(1))] = float(m.group(2))
    return vals


# ------------------------------------------------------------------ concept
def concept():
    model, m0 = load("kamchatka")
    t, x, y, T = model[:, 0], model[:, 1], model[:, 2], model[:, 3]
    n = len(t)
    cols = plt.cm.viridis(np.linspace(0.05, 0.9, n))
    # rupture direction: from the hypocenter (subevent 1 at 0,0) toward the
    # moment-weighted centroid; station A lies in that direction
    cen = np.average(np.c_[x, y], axis=0, weights=m0)
    u = cen / np.linalg.norm(cen)
    azA = np.degrees(np.arctan2(u[0], u[1])) % 360            # azimuth from north
    azB = (azA + 180) % 360
    c = 15.0          # apparent horizontal speed of teleseismic P across the source (km/s)

    fig = plt.figure(figsize=(14, 4.8))
    ax1 = fig.add_axes([0.065, 0.13, 0.23, 0.72])
    ax2 = fig.add_axes([0.37, 0.13, 0.27, 0.72])
    ax3 = fig.add_axes([0.70, 0.13, 0.28, 0.72])

    for i in range(n):
        ax1.scatter(x[i], y[i], s=40 + 900 * m0[i] / m0.max(), color=cols[i], ec="k", lw=0.6, zorder=3)
    # numbers beside the circles; co-located subevents get opposite sides
    for i in range(n):
        near = [j for j in range(n) if j != i and np.hypot(x[j] - x[i], y[j] - y[i]) < 25]
        side = 1 if not near or i < min(near) else -1
        ax1.annotate(f"{i+1}", (x[i], y[i]), xytext=(side * 14, 0), textcoords="offset points",
                     ha="center", va="center", fontsize=10, fontweight="bold", zorder=4)
    ax1.plot(0, 0, "k*", ms=13, zorder=5)
    span = max(np.ptp(x), np.ptp(y)) * 0.75
    mid = np.array([x.mean(), y.mean()])
    reach = np.max(np.abs((np.c_[x, y] - mid) @ u)) + 0.42 * span
    for sgn, lab in ((1, "to station A"), (-1, "to station B")):
        a0 = mid + sgn * reach * u
        a1 = mid + sgn * (reach + 0.35 * span) * u
        ax1.annotate("", a1, a0, arrowprops=dict(arrowstyle="-|>", lw=2, color="0.25"))
        ax1.annotate(lab, a1, xytext=(sgn * u[0] * 30, sgn * u[1] * 30), textcoords="offset points",
                     ha="center", va="center", fontsize=10, color="0.25")
    lim = reach + 0.75 * span
    ax1.set_xlim(mid[0] - lim, mid[0] + lim); ax1.set_ylim(mid[1] - lim, mid[1] + lim)
    ax1.set_aspect("equal"); ax1.set_xlabel("east (km)"); ax1.set_ylabel("north (km)")
    ax1.set_title("1. A big earthquake = a few\n    subevents (★ = where it started)", loc="left", fontsize=12)

    tt = np.linspace(-10, t.max() + T.max() + 20, 2000)
    tot = np.zeros_like(tt)
    for i in range(n):
        p = pulse(tt, t[i], T[i], m0[i]); tot += p
        ax2.fill_between(tt, 0, p, color=cols[i], alpha=0.75, lw=0)
        ax2.text(t[i], p.max() * 1.03, str(i + 1), ha="center", fontsize=9)
    ax2.plot(tt, tot, "k", lw=1.6, label="total")
    ax2.set_yticks([]); ax2.set_xlabel("time after the earthquake started (s)")
    ax2.set_ylabel("moment release rate")
    ax2.set_title("2. Each subevent releases its moment\n    as a smooth pulse (time, duration, size)",
                  loc="left", fontsize=12)
    ax2.legend(frameon=False, loc="upper right")

    for az, lab, ls, off in ((azA, "station A: rupture runs toward it", "-", 0.0),
                             (azB, "station B: rupture runs away from it", "--", 0.0)):
        dx, dy = np.sin(np.radians(az)), np.cos(np.radians(az))
        shift = -(x * dx + y * dy) / c            # earlier arrival when closer
        app = sum(pulse(tt, t[i] + shift[i], T[i], m0[i]) for i in range(n))
        ax3.plot(tt, app, ls, color="k", lw=1.6, label=lab)
    ax3.set_yticks([]); ax3.set_xlabel("time at the station after the first arrival (s)")
    ax3.set_title("3. Stations in different directions see the\n    pulses compressed or stretched in time",
                  loc="left", fontsize=12)
    ax3.set_ylim(0, None)
    ax3.set_ylim(0, ax3.get_ylim()[1] * 1.32)
    ax3.legend(frameon=False, loc="upper center", fontsize=9.5, ncol=1)
    ax3.text(0.02, -0.25, "Fitting the waveforms at ~100-200 stations at once pins down where, when and "
             "how each subevent broke.", transform=ax3.transAxes, fontsize=10, color="0.3")
    fig.savefig(os.path.join(OUT, "concept.png"), dpi=150)
    plt.close(fig)


# ------------------------------------------------------------------ model selection
def model_selection():
    """One panel per example: best misfit against the number of subevents, the selection (star) and
    the dashed line of a misfit 5% lower per added subevent, the cost of a subevent (README section
    4.6). The curve touches the line at the star and stays above it elsewhere, except that a number
    tied with the star can lie up to 0.5% below it. The 5% is 2 k sigma^2 of
    programs/prepare_SubMIT/plot_lcurve.py; the selection is read from examples/NAME/selected_n,
    which that script decided."""
    log_misfit_drop_per_subevent = 0.05
    line_colour = "#2a78d6"
    figure, axes = plt.subplots(1, len(EVENTS), figsize=(2.6 * len(EVENTS), 2.9))
    for axis, (example_name, example_label) in zip(axes, EVENTS):
        best_misfit_of = lcurve(example_name)
        subevent_counts = np.array(sorted(best_misfit_of))
        best_misfits = np.array([best_misfit_of[count] for count in subevent_counts])
        selected_count = int(open(os.path.join(EX, example_name, "selected_n")).read())
        selected_misfit = best_misfit_of[selected_count]
        line_counts = np.linspace(subevent_counts[0] - 0.3, subevent_counts[-1] + 0.3, 50)
        line_misfits = selected_misfit * np.exp(-log_misfit_drop_per_subevent * (line_counts - selected_count))
        axis.plot(line_counts, line_misfits, "--", color="0.55", lw=0.8, zorder=1)
        axis.plot(subevent_counts, best_misfits, "-o", color=line_colour, lw=1.6, ms=4.5)
        axis.plot(selected_count, selected_misfit, "*", ms=15, color=line_colour, mec="white", mew=0.8, zorder=5)
        axis.set_xticks(subevent_counts)
        axis.set_xlim(line_counts[0], line_counts[-1])
        axis.set_ylim(min(best_misfits.min(), line_misfits.min()) * 0.95, best_misfits.max() * 1.04)
        axis.tick_params(labelsize=9)
        axis.grid(axis="y", ls="--", lw=0.4, color="0.85")
        axis.set_axisbelow(True)
        axis.set_xlabel("subevents", fontsize=10)
        # the event name without magnitude and depth: the README text gives those
        axis.text(0.97, 0.95, example_name.capitalize(), transform=axis.transAxes, ha="right", va="top", fontsize=10)
    axes[0].set_ylabel("penalized misfit", fontsize=10)
    figure.tight_layout()
    figure.savefig(os.path.join(OUT, "model_selection.png"), dpi=300)
    plt.close(figure)


# ------------------------------------------------------------------ results overview
def results_overview():
    fig, axs = plt.subplots(1, len(EVENTS), figsize=(3.3 * len(EVENTS), 3.4))
    for ax, (name, lab) in zip(axs, EVENTS):
        model, m0 = load(name)
        t, T = model[:, 0], model[:, 3]
        n = len(t)
        cols = plt.cm.viridis(np.linspace(0.05, 0.9, n))
        tt = np.linspace(-0.05 * (t.max() + T.max()), 1.15 * (t.max() + T.max()), 1500)
        tot = np.zeros_like(tt)
        for i in range(n):
            p = pulse(tt, t[i], T[i], m0[i]) * 1e27; tot += p
            ax.fill_between(tt, 0, p, color=cols[i], alpha=0.75, lw=0)
        ax.plot(tt, tot, "k", lw=1.4)
        mw = 2 / 3 * np.log10(m0.sum() * 1e27) - 10.7
        ax.set_title(lab.replace("  (", "\n(").replace("  Mw", "\nMw"), fontsize=10, loc="left")
        ax.text(0.97, 0.93, f"{n} subevents", transform=ax.transAxes, ha="right", fontsize=10,
                fontweight="bold")
        ax.set_yticks([]); ax.set_xlabel("time (s)")
        ax.set_ylim(0, tot.max() * 1.12)
    axs[0].set_ylabel("moment rate")
    fig.suptitle("Moment-rate functions of the example earthquakes (colored = individual subevents)",
                 x=0.01, ha="left", fontsize=12)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT, "results_overview.png"), dpi=150)
    plt.close(fig)


if __name__ == "__main__":
    concept()
    model_selection()
    results_overview()
    for name, _ in EVENTS:
        ref = os.path.join(EX, name, "reference")
        for src, dst in (("subevents", "subevents"), ("stations", "stations"), ("histoplot_py", "hist"),
                         ("fits_P", "fits_P"), ("fits_SH", "fits_SH"), ("fits_rayl", "fits_rayl")):
            shutil.copy(os.path.join(ref, f"{src}.png"), os.path.join(OUT, f"{name}_{dst}.png"))
    print("wrote", sorted(os.listdir(OUT)))

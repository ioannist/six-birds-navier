#!/usr/bin/env python3
"""Render the main-text paper figures from the manifest NPZ artifacts.

This script performs no simulation. Several figures need slow producers, so it requires a
prior full run of `ns_paper_figures.py --include-slow` (the tracked artifacts satisfy this). It reads the companion .npz files written by the
manifest producers (python/scripts/ns_paper_figures.py) and draws publication versions
into paper/figures/. Every plotted number comes from those arrays.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402

ROOT = Path(__file__).resolve().parents[2]
ART = ROOT / "python" / "artifacts" / "paper"
OUT = ROOT / "paper" / "figures"

BLUE = "#2a78d6"
ORANGE = "#eb6834"
AQUA = "#1baf7a"
INK = "#0b0b0b"
MUTED = "#52514e"
GRID = "#d9d8d4"

plt.rcParams.update(
    {
        "font.size": 9.5,
        "axes.edgecolor": MUTED,
        "axes.labelcolor": INK,
        "axes.titlesize": 9.5,
        "axes.titlecolor": INK,
        "xtick.color": MUTED,
        "ytick.color": MUTED,
        "legend.frameon": False,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "pdf.fonttype": 42,
    }
)


def _load(name: str) -> np.lib.npyio.NpzFile:
    path = ART / f"{name}_data.npz"
    if not path.exists():
        raise SystemExit(
            f"missing artifact {path}; run python/scripts/ns_paper_figures.py "
            "--include-slow --include-optional first (several main-text figures use slow producers)"
        )
    return np.load(path)


def _grid(ax) -> None:
    ax.grid(True, color=GRID, lw=0.6, alpha=0.8)
    ax.set_axisbelow(True)


def toy_zeno() -> Path:
    d = _load("nsbu12_toy_zeno_gallery")
    J = d["j_plot"]
    fig, ax = plt.subplots(figsize=(6.4, 3.0))
    series = [
        (d["z1"], BLUE, r"$\mathrm{Cap}(j)=2^j$, $w=1$: DIV fails"),
        (d["z2"], ORANGE, r"$\mathrm{Cap}(j)=1$, $w(j)=2^{-j}$: WORK fails"),
        (d["nz"], AQUA, r"$\mathrm{Cap}(j)=j+1$, $w=1$: No-Zeno"),
    ]
    for (y, col, lab), ls in zip(series, ["-", (0, (4, 3)), "-"]):
        ax.plot(J, y, color=col, lw=2, ls=ls, label=lab)
    ax.set_xscale("log")
    ax.axhline(2.0, color=MUTED, lw=0.8, ls=":")
    ax.text(J[-1], 2.0, "  limit 2", va="center", ha="left", color=MUTED, fontsize=8.5)
    ax.set_xlabel(r"number of levels $J$")
    ax.set_ylabel(r"crossing time $t_J$")
    ax.set_xlim(J[0], J[-1])
    ax.set_ylim(0, None)
    _grid(ax)
    ax.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    out = OUT / "toy_zeno_gallery.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out


def legal_vs_illegal() -> Path:
    d = _load("ns_spt_illegal_vs_legal")
    j = d["j_legal"]
    if not np.array_equal(j, d["j_illegal"]):
        raise SystemExit("legal and illegal shell indices differ")
    kl, ki = d["kappa_legal"], d["kappa_illegal"]
    thr = float(d["legal_threshold"])
    fig, ax = plt.subplots(figsize=(6.0, 3.0))
    width = 0.36
    floor = 0.3
    for off, vals, col, lab in [(-width / 2, kl, BLUE, "legal baseline"),
                                (width / 2, ki, ORANGE, "localized wavepacket")]:
        shown = np.where(vals > 0, vals, np.nan)
        ax.bar(j + off, shown - floor, bottom=floor, width=width * 0.92, color=col, label=lab)
    ax.axhline(thr, color=INK, lw=1, ls="--")
    ax.text(j[0] - 0.45, thr * 1.12, rf"threshold $\kappa_\star={thr:g}$", color=INK, fontsize=8.5)
    imax = int(np.argmax(ki))
    ax.text(j[imax] + width / 2, ki[imax] * 1.25, rf"$\kappa_{{{int(j[imax])}}}\approx{ki[imax]:.0f}$",
            ha="center", va="bottom", color=INK, fontsize=8.5)
    ax.set_yscale("log")
    ax.set_ylim(floor, 2e3)
    ax.set_xticks(j)
    ax.set_xlabel(r"dyadic shell $j$")
    ax.set_ylabel(r"concentration factor $\kappa_j$")
    _grid(ax)
    ax.legend(loc="upper left", fontsize=8.5)
    fig.tight_layout()
    out = OUT / "spt_illegal_vs_legal.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out


def channel_diagnostics() -> Path:
    c = _load("ns_channel_count")
    d = _load("ns_channel_delocalization")
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 3.0))
    ax = axes[0]
    jc = c["j_values"]
    ax.plot(jc, c["m95_pooled"], "o-", color=BLUE, lw=2, ms=5, label="pooled over time and seeds")
    ax.plot(jc, c["m95_window_median"], "s-", color=ORANGE, lw=2, ms=5, label="median over windows")
    ax.fill_between(jc, c["m95_window_p25"], c["m95_window_p75"], color=ORANGE, alpha=0.18, lw=0)
    ax.set_xlabel(r"dyadic shell $j$")
    ax.set_ylabel(r"rank $m_{95}$")
    ax.set_xticks(jc)
    top = float(max(np.max(c["m95_pooled"]), np.max(c["m95_window_p75"]))) + 1.5
    ax.set_ylim(0, top)
    ax.set_yticks(np.arange(0, int(top) + 1))
    _grid(ax)
    ax.legend(loc="lower right", fontsize=8)
    ax.set_title("(a) effective channel count", loc="left")

    ax = axes[1]
    jd = d["j_values"]
    ax.semilogy(jd, d["eta_hat"], "o-", color=BLUE, lw=2, ms=5, label=r"channel span, $\hat\eta_j$")
    ax.semilogy(jd, d["eta_snapshot_bound"], "s--", color=ORANGE, lw=2, ms=5,
                label="full snapshots, residual charged")
    ax.set_xlabel(r"dyadic shell $j$")
    ax.set_ylabel("local energy fraction bound")
    ax.set_xticks(jd)
    _grid(ax)
    ax.legend(loc="lower left", fontsize=8)
    ax.set_title("(b) anti-localization bounds", loc="left")
    fig.tight_layout()
    out = OUT / "channel_diagnostics.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out


def legal_sweep() -> Path:
    d = _load("ns_spt_legal_sweep")
    j = d["j_sorted"]
    kap, burst = d["worst_kappa"], d["worst_burst"]
    occupied = kap > 0
    meta = __import__("json").loads(bytes(d["meta_json"]).decode())
    fig, axes = plt.subplots(1, 2, figsize=(7.2, 2.9))
    panels = [
        (axes[0], kap, r"worst $\kappa_j$ over sweep", "(a) concentration factor", float(meta["kappa_threshold"]), True),
        (axes[1], burst, "worst burst ratio over sweep", "(b) burst ratio", float(meta["burst_threshold"]), False),
    ]
    for ax, vals, ylab, title, thr, show_thr in panels:
        ax.plot(j[occupied], vals[occupied], "o-", color=BLUE, lw=2, ms=5, label="occupied shells")
        if np.any(~occupied):
            ax.plot(j[~occupied], vals[~occupied], "o", mfc="white", mec=MUTED, ms=6, label="empty shell")
        if show_thr:
            ax.axhline(thr, color=INK, lw=1, ls="--")
            ax.text(j[0], thr * 0.93, rf"threshold {thr:g}", va="top", fontsize=8.5, color=INK)
            ax.set_ylim(0, thr * 1.1)
        else:
            ax.set_ylim(0, max(1.5, float(np.max(vals)) * 1.25))
            ax.text(j[0], ax.get_ylim()[1] * 0.95, rf"threshold {thr:g} (off scale)", va="top",
                    fontsize=8.5, color=MUTED)
        ax.set_xlabel(r"dyadic shell $j$")
        ax.set_ylabel(ylab)
        ax.set_xticks(j)
        _grid(ax)
        ax.set_title(title, loc="left")
    axes[0].legend(loc="center left", fontsize=8)
    fig.tight_layout()
    out = OUT / "spt_legal_sweep.pdf"
    fig.savefig(out)
    plt.close(fig)
    return out


def mu_convergence() -> Path:
    d = _load("ns16_mu_convergence_3d")
    mu, err = d["mu_vals"], d["err_vals"]
    slope = float(np.polyfit(np.log(mu), np.log(err), 1)[0])
    fig, ax = plt.subplots(figsize=(3.6, 2.9))
    ax.loglog(mu, err, "o-", color=BLUE, lw=2, ms=5, label="measured")
    ref = err[-1] * mu / mu[-1]
    ax.loglog(mu, ref, ":", color=MUTED, lw=1.4, label="slope 1")
    ax.set_xlabel(r"hyperviscosity $\mu$")
    ax.set_ylabel("relative low-mode difference")
    ax.text(0.04, 0.9, f"fitted slope {slope:.2f}", transform=ax.transAxes, color=MUTED, fontsize=8.5)
    _grid(ax)
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    out = OUT / "mu_convergence.pdf"
    fig.savefig(out)
    plt.close(fig)
    print(f"mu_convergence fitted slope = {slope:.4f}")
    return out


FIGURES = {
    "toy_zeno": toy_zeno,
    "legal_vs_illegal": legal_vs_illegal,
    "channel_diagnostics": channel_diagnostics,
    "legal_sweep": legal_sweep,
    "mu_convergence": mu_convergence,
}


def main() -> int:
    global ART
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", nargs="*", choices=sorted(FIGURES), default=None)
    parser.add_argument("--artifact-dir", default=str(ART), help="directory holding the *_data.npz files")
    args = parser.parse_args()
    ART = Path(args.artifact_dir).resolve()
    OUT.mkdir(parents=True, exist_ok=True)
    for name in args.only or FIGURES:
        print(f"wrote {FIGURES[name]().relative_to(ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

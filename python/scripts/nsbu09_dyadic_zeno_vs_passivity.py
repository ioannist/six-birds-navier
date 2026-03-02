import sys
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt

ROOT = Path(__file__).resolve().parents[2]
SRC = ROOT / "python" / "src"
sys.path.insert(0, str(SRC))

from nswave.dyadic_models import simulate_kp


ART = ROOT / "python" / "artifacts"


def main() -> int:
    J = 20
    dt = 1e-4
    t_max = 0.3
    nu = 1e-3
    alpha = 1.25
    lam = 2.0
    mu_list = [-1e-5, -5e-6, 0.0, 5e-6, 1e-5]

    u0 = np.zeros(J)
    u0[0] = 1.0
    u0[1] = 0.5

    forcing = np.zeros(J)
    forcing[0] = 0.2

    ART.mkdir(parents=True, exist_ok=True)

    print("NS-BU-09 dyadic Zeno vs passivity")
    print(f"J={J} dt={dt} t_max={t_max} nu={nu} alpha={alpha} lam={lam}")
    print("forcing: f0=0.2 (low-shell)")
    print("mu | zeno_gamma | n_crossings | mean_P_sgs | passive_flag | notes")

    zeno_vals = []
    mu_vals = []

    for mu in mu_list:
        res = simulate_kp(J=J, dt=dt, t_max=t_max, nu=nu, mu=mu, alpha=alpha, u0=u0, lam=lam, forcing=forcing)
        t = res["t"]
        u_hist = res["u_hist"]
        E_hist = res["E_hist"]
        P_hist = res["sgs_power_hist"]

        thresh = 1e-8 * E_hist[0]
        E_shell = 0.5 * u_hist * u_hist
        j_star = np.array([np.max(np.where(E_shell[i] >= thresh)[0]) if np.any(E_shell[i] >= thresh) else -1 for i in range(E_shell.shape[0])])
        j_mon = np.maximum.accumulate(j_star)

        j_vals = np.arange(np.min(j_mon), np.max(j_mon) + 1)
        t_reach = []
        for j in j_vals:
            hit = np.where(j_mon >= j)[0]
            if hit.size:
                t_reach.append(t[hit[0]])
        t_reach = np.array(t_reach)

        delta_t = np.diff(t_reach)
        notes = []
        if t_reach.size < len(j_vals):
            notes.append("plateau")

        if delta_t.size >= 4 and np.all(delta_t > 0):
            j_mid = j_vals[1:]
            y = np.log2(delta_t)
            x = j_mid[: len(y)]
            slope, intercept = np.polyfit(x, y, 1)
            gamma = max(-slope, 0.0)
        else:
            gamma = float("nan")
            notes.append("few-crossings")

        mean_P = float(np.mean(P_hist))
        passive_flag = mean_P >= 0

        zeno_vals.append(gamma)
        mu_vals.append(mu)

        print(
            f"{mu: .1e} | {gamma:9.4f} | {delta_t.size:11d} | {mean_P:10.4e} | {passive_flag!s:11} | {';'.join(notes) if notes else '-'}"
        )

    plt.figure(figsize=(6, 4))
    plt.plot(mu_vals, zeno_vals, marker="o")
    plt.axhline(0.0, color="k", linestyle="--", linewidth=0.8)
    plt.xlabel("mu")
    plt.ylabel("zeno_gamma")
    plt.title("NS-BU-09: Zeno diagnostic vs SGS sign")
    plt.tight_layout()

    out_path = ART / "nsbu09_dyadic_zeno.png"
    plt.savefig(out_path, dpi=150)
    print(f"saved {out_path}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

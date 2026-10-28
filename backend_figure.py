"""Figure 4 for Section 5.5, drawn from results.json."""
import json, numpy as np, matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from erlang import mmc
R = json.load(open("results.json"))
plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 7.6, "axes.spines.top": False, "axes.spines.right": False})
DARK, MID, ACC, GREY = "#1f3b57", "#5b7a99", "#b5532e", "#9aa9b8"
fig, ax = plt.subplots(2, 1, figsize=(3.5, 6.2), dpi=300)
cvs = [0.0, 0.5, 1.0, 2.0]
for lam, ls, mk in ((90, "-", "o"), (110, "--", "s")):
    reg = [100 * next(s for s in R["S"] if s["lam"] == lam and s["cv"] == c)["reg"]["attain"] for c in cvs]
    pool = [100 * next(s for s in R["S"] if s["lam"] == lam and s["cv"] == c)["pool"]["attain"] for c in cvs]
    ax[0].plot(cvs, reg, ls, marker=mk, color=GREY, ms=3.5, lw=1.3, label=f"regional 2×4, λ = {lam}")
    ax[0].plot(cvs, pool, ls, marker=mk, color=DARK, ms=3.5, lw=1.3, label=f"pooled 1×8, λ = {lam}")
ax[0].axvline(1.0, color=ACC, lw=0.7, ls=":")
ax[0].set_xticks(cvs); ax[0].set_xlabel("Service-time coefficient of variation"); ax[0].set_ylabel("Deadline attainment ≤ 250 ms (%)")
ax[0].set_ylim(45, 101); ax[0].legend(frameon=False, fontsize=6.3, loc="lower left")
ax[0].set_title("(a) Service-time variability", fontsize=8.4, loc="left")
labels, vals, cols = [], [], []
def add(l, v, c): labels.append(l); vals.append(v); cols.append(c)
I = {i["split"]: i for i in R["I"]}; N = {n["prefer_local"]: n for n in R["N"]}
add("Regional, 50/50", I[0.5]["reg"]["mean_ms"], GREY)
add("Regional, 60/40", I[0.6]["reg"]["mean_ms"], GREY)
add("Pooled, ideal", I[0.5]["pool"]["mean_ms"], DARK)
add("Pooled, 30 ms on all\n(Section 5.2)", mmc(90, 15, 8, .25)["R_ms"] + 30, ACC)
add(f"Pooled, 30 ms cross-region\nonly ({100*N[False]['cross_frac']:.0f}% cross)", N[False]["mean_ms"], MID)
add(f"Pooled, locality-preferring\n({100*N[True]['cross_frac']:.0f}% cross)", N[True]["mean_ms"], MID)
y = np.arange(len(vals))[::-1]
ax[1].barh(y, vals, color=cols, height=0.62)
for yi, v in zip(y, vals): ax[1].text(v + 2, yi, f"{v:.1f}", va="center", fontsize=6.6)
ax[1].set_yticks(y, labels, fontsize=6.4); ax[1].set_xlabel("Mean inference response at λ = 90 (ms)"); ax[1].set_xlim(0, 175)
ax[1].set_title("(b) Imbalance and locality", fontsize=8.4, loc="left")
fig.tight_layout(h_pad=1.4); fig.savefig("fig4.png"); print("fig4.png written")

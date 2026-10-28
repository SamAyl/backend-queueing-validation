"""Discrete-event validation of the Erlang C analysis (Section 5.5 additions).

FCFS multi-server queues are simulated request by request. Regional partitioning = two independent
4-replica queues; pooled dispatch = one 8-replica queue. Experiments:
 V  exponential service (validates Table 3 against Erlang C)
 S  service-time variability (coefficient of variation 0, 0.5, 1, 2; lognormal / deterministic)
 I  regional load imbalance (50/50, 60/40, 70/30)
 B  bursty arrivals (two-state Markov-modulated Poisson process with the same mean rate)
 N  locality: 30 ms penalty only for requests served by the other region, random vs locality-preferring dispatch
Run: python backend_sim.py -> results.json
"""
import heapq, json, math
import numpy as np
from erlang import mmc

MU, T_DEAD, N_REQ, REPS = 15.0, 0.250, 200_000, 5
rng = np.random.default_rng(31)

def service(n, cv):
    m = 1 / MU
    if cv == 0: return np.full(n, m)
    if cv == 1: return rng.exponential(m, n)
    s2 = math.log(1 + cv * cv)
    return rng.lognormal(math.log(m) - s2 / 2, math.sqrt(s2), n)

def poisson_arrivals(lam, n):
    return np.cumsum(rng.exponential(1 / lam, n))

def mmpp_arrivals(lam, n, hi=1.5, lo=0.5, sojourn=10.0):
    t, out, state = 0.0, [], 0
    end_state = rng.exponential(sojourn)
    while len(out) < n:
        rate = lam * (hi if state == 0 else lo)
        dt = rng.exponential(1 / rate)
        if t + dt > end_state:
            t = end_state; state ^= 1; end_state = t + rng.exponential(sojourn); continue
        t += dt; out.append(t)
    return np.array(out)

def fcfs(arr, svc, c, penalty=None):
    """Returns response times; penalty(i, server) adds transfer time for cross-region service."""
    free = [0.0] * c; heapq.heapify(free)
    resp = np.empty(len(arr))
    for i, (a, s) in enumerate(zip(arr, svc)):
        f = heapq.heappop(free); start = max(a, f)
        heapq.heappush(free, start + s); resp[i] = start + s - a
    return resp

def fcfs_regions(arr, origin, svc, c_per=4, penalty=0.030, prefer_local=False):
    """Pooled 8-replica FCFS dispatch across two regions. A cross-region request reaches its server
    after `penalty` seconds. Without locality preference the earliest-free server is chosen regardless of
    region; with it, the server giving the earliest service start including the penalty is chosen."""
    free = [[0.0] * c_per for _ in range(2)]
    resp = np.empty(len(arr)); cross = 0
    for i in range(len(arr)):
        a, o, s = arr[i], origin[i], svc[i]
        best = None
        for r in (0, 1):
            fr = free[r]; j = min(range(c_per), key=fr.__getitem__); f = fr[j]
            p = penalty if r != o else 0.0
            start = max(a + p, f)
            key = start if prefer_local else max(a, f) + (0.0 if r == o else 1e-9)
            if best is None or key < best[0]:
                best = (key, r, j, start)
        _, r, j, start = best
        free[r][j] = start + s; resp[i] = start + s - a; cross += r != o
    return resp, cross / len(arr)

def stats(resp, T=T_DEAD):
    return dict(mean_ms=1000 * float(resp.mean()), attain=float((resp <= T).mean()), p99_ms=1000 * float(np.quantile(resp, 0.99)))

def regional(lam, cv=1, split=0.5, arrivals="poisson"):
    out = []
    for share in (split, 1 - split):
        l = lam * share; n = int(N_REQ * share)
        a = poisson_arrivals(l, n) if arrivals == "poisson" else mmpp_arrivals(l, n)
        out.append(fcfs(a, service(n, cv), 4))
    return np.concatenate(out)
def pooled(lam, cv=1, arrivals="poisson"):
    a = poisson_arrivals(lam, N_REQ) if arrivals == "poisson" else mmpp_arrivals(lam, N_REQ)
    return fcfs(a, service(N_REQ, cv), 8)
def agg(fn, *args, **kw):
    runs = [stats(fn(*args, **kw)) for _ in range(REPS)]
    m = {k: float(np.mean([r[k] for r in runs])) for k in runs[0]}
    m.update({k + "_ci": float(1.96 * np.std([r[k] for r in runs], ddof=1) / math.sqrt(REPS)) for k in runs[0]})
    return m

OUT = {}
# V: validation against Erlang C
OUT["V"] = []
for lam in (40, 60, 80, 90, 110):
    er, ep = mmc(lam / 2, MU, 4, T_DEAD), mmc(lam, MU, 8, T_DEAD)
    sr, sp = agg(regional, lam), agg(pooled, lam)
    OUT["V"].append(dict(lam=lam, erl_reg_ms=er["R_ms"], sim_reg_ms=sr["mean_ms"], sim_reg_ms_ci=sr["mean_ms_ci"],
                         erl_pool_ms=ep["R_ms"], sim_pool_ms=sp["mean_ms"], sim_pool_ms_ci=sp["mean_ms_ci"],
                         erl_reg_att=er["attain"], sim_reg_att=sr["attain"], erl_pool_att=ep["attain"], sim_pool_att=sp["attain"]))
    print("V", OUT["V"][-1])
# S: service variability
OUT["S"] = []
for lam in (90, 110):
    for cv in (0.0, 0.5, 1.0, 2.0):
        sr, sp = agg(regional, lam, cv=cv), agg(pooled, lam, cv=cv)
        # Lee-Longton style approximation: Wq(M/G/c) ~ (1+cv^2)/2 * Wq(M/M/c)
        er, ep = mmc(lam / 2, MU, 4, T_DEAD), mmc(lam, MU, 8, T_DEAD)
        f = (1 + cv * cv) / 2
        OUT["S"].append(dict(lam=lam, cv=cv, reg=sr, pool=sp, approx_reg_ms=er["wq_ms"] * f + 1000 / MU, approx_pool_ms=ep["wq_ms"] * f + 1000 / MU))
        print("S", lam, cv, round(sr["mean_ms"], 1), round(sp["mean_ms"], 1), round(100 * sr["attain"], 1), round(100 * sp["attain"], 1))
# I: imbalance
OUT["I"] = []
for split in (0.5, 0.6, 0.7):
    sr, sp = agg(regional, 90, split=split), agg(pooled, 90)
    OUT["I"].append(dict(split=split, reg=sr, pool=sp)); print("I", split, round(sr["mean_ms"], 1), round(100 * sr["attain"], 1), round(sp["mean_ms"], 1))
# B: bursty arrivals
OUT["B"] = []
for lam in (80, 90):
    sr, sp = agg(regional, lam, arrivals="mmpp"), agg(pooled, lam, arrivals="mmpp")
    pr, pp = agg(regional, lam), agg(pooled, lam)
    OUT["B"].append(dict(lam=lam, reg=sr, pool=sp, reg_poisson=pr, pool_poisson=pp)); print("B", lam, round(sr["mean_ms"], 1), round(sp["mean_ms"], 1), round(100 * sr["attain"], 1), round(100 * sp["attain"], 1))
# N: locality penalty only on cross-region service
OUT["N"] = []
for prefer in (False, True):
    runs, crosses = [], []
    for _ in range(REPS):
        a = poisson_arrivals(90, N_REQ); o = (rng.random(N_REQ) < 0.5).astype(int)
        r, cr = fcfs_regions(a, o, service(N_REQ, 1), prefer_local=prefer); runs.append(stats(r)); crosses.append(cr)
    m = {k: float(np.mean([x[k] for x in runs])) for k in runs[0]}; m["cross_frac"] = float(np.mean(crosses))
    OUT["N"].append(dict(prefer_local=prefer, **m)); print("N", prefer, m)
json.dump(OUT, open("results.json", "w"), indent=1)
print("done")

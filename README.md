# Scalable Backend Architectures for AI-Driven Enterprise Software — validation code

Code and outputs for the Erlang C calculations and the discrete-event simulation (Section 5.5) of:

> S. Ayloo, *Scalable Backend Architectures for AI-Driven Enterprise Software Applications* (manuscript under submission, 2026).

| File | Paper | What it does |
|---|---|---|
| `erlang.py` | Eqs. 1–3, Tables 3–4 | Exact Erlang C mean response and deadline attainment; prints every analytical value in the paper |
| `backend_sim.py` | Section 5.5, Tables 5–7 | FCFS multi-server discrete-event simulation: validation (V), service-time variability (S), regional imbalance (I), bursty MMPP arrivals (B), cross-region locality penalty (N) |
| `backend_figure.py` | Figure 4 | Draws the figure from `results.json` |

## Run

```bash
pip install -r requirements.txt
python erlang.py         # analytical values (instant)
python backend_sim.py    # simulation, ~5 min -> results.json
python backend_figure.py # figures/fig4.png
```

Seeds are fixed; `results.json` holds the outputs reported in the paper.
These are models with stated parameters, not measurements of an enterprise deployment.

## License

MIT

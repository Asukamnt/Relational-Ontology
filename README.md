# Environmental correction to the Hubble tension (public minimal package)

This repository is a **minimal public release** focusing only on the *Hubble-tension environment-bias* pipeline:

- a runnable dataset + fitting code (`math/hubble_env_scan/`)
- a tiny in-repo toy chain for sign / linearity (`math/spectral_graph/run_relational_hubble_beta_poisson_3d_v0.py`)
- an arXiv-ready LaTeX draft + build script (`papers/`)

It is intentionally narrow in scope: **a reproducible, falsifiable research note** (not a claim of a solved tension).

## TL;DR (model)

We parameterize an observer-environment effect:

\[
H_0(\delta)=H_{0,\mathrm{ref}}\,(1-\beta\,\delta).
\]

The anchor-free, “hard content” is the ratio law:

\[
\frac{H_{0,A}}{H_{0,B}}=\frac{1-\beta\,\delta_A}{1-\beta\,\delta_B}.
\]

Current reproducible fit (v1 direct convention; see JSON for details): **β ≈ 0.171** (MC 16/50/84: 0.139 / 0.169 / 0.201).

## Quick reproduce (3 commands)

From the repository root:

```bash
pip install -r requirements.txt
python math/hubble_env_scan/run_hubble_env_scan_v1.py
python math/spectral_graph/run_relational_hubble_beta_poisson_3d_v0.py
```

Expected outputs:

- `math/hubble_env_scan/results_v1.json`
- `math/hubble_env_scan/plot_v1.png`
- `math/spectral_graph/relational_hubble_beta_poisson_3d_results_v0.json`
- `math/spectral_graph/relational_hubble_beta_poisson_3d_plot_v0.png`

## Where to read

- One-page summary: `ONEPAGER_HUBBLE_TENSION.md`
- English draft: `papers/hubble_tension_environment_v0_en.md`
- arXiv LaTeX draft: `papers/hubble_tension_environment_arxiv_v0.tex`
- Repro steps: `REPRODUCE_HUBBLE_TENSION.md`

## License / citation

- License: MIT (`LICENSE`)
- Citation metadata: `CITATION.cff`


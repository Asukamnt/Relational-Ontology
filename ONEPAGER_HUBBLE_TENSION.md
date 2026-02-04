# One-page summary (v0): Environmental correction to the Hubble tension

**Scope.** This is a *reproducible, falsifiable research note* (not a claim of a solved tension). The goal is to freeze a minimal hypothesis and a runnable pipeline so that new data can decisively support, refine, or refute it.

---

## Hypothesis (minimal parameterization)

We model the locally inferred Hubble constant as an observer-environment effect:

\[
H_0(\delta)=H_{0,\mathrm{ref}}\,(1-\beta\,\delta),
\]

where:

- \(H_{0,\mathrm{ref}}\): reference expansion rate for a statistically typical observer (\(\delta=0\)); in practice anchored to CMB / Planck 2018 (this is a convention / unit scale, not a “God’s-eye” quantity)
- \(\delta\): large-scale density contrast around the observer (definition/scale depends on source)
- \(\beta>0\): density–expansion coupling to be estimated from multi-environment points

**Relational form (recommended).** The hard, anchor-free content is the relative bias law:
\[
\frac{H_{0,A}}{H_{0,B}}=\frac{1-\beta\,\delta_A}{1-\beta\,\delta_B}.
\]

**Interpretation.**

- Underdense environment (\(\delta<0\)) \(\Rightarrow\) higher inferred local \(H_0\)
- Overdense environment (\(\delta>0\)) \(\Rightarrow\) lower inferred local \(H_0\)

---

## Current reproducible numbers (v1 “direct” convention)

Pipeline output (fit points n=3; details in JSON):

- Anchor: \(H_{0,\mathrm{ref}} = 67.36\) km/s/Mpc (Planck 2018)
- Fit: \(\beta \approx 0.171 \pm 0.024\)
- MC uncertainty (16/50/84): **0.139 / 0.169 / 0.201**

Key figure:

- `math/hubble_env_scan/plot_v1.png`

---

## Frozen forecast points (example; stored in `results_v1.json`)

Note on absolute numbers. The km/s/Mpc values below assume a fixed reference anchor
\(H_{0,\mathrm{ref}}=67.36\) km/s/Mpc (Planck 2018). This is a *convention* used to
express the relational law in physical units. The anchor-free, relational content
is the ratio law above; for comparisons across different anchoring conventions,
use \(H_0/H_{0,\mathrm{ref}}\) or pairwise ratios.

Using \(\beta \approx 0.171\) (direct convention), the pipeline freezes example predictions:

- Cold Spot-like underdensity:
  - \(\delta=-0.70\) \(\Rightarrow\) \(H_0\approx 75.4\)
  - \(\delta=-0.60\) \(\Rightarrow\) \(H_0\approx 74.3\)
- Typical cluster:
  - \(\delta=+0.30\) \(\Rightarrow\) \(H_0\approx 63.9\)
  - \(\delta=+0.50\) \(\Rightarrow\) \(H_0\approx 61.6\)
- Coma-like strong overdensity:
  - \(\delta=+1.00\) \(\Rightarrow\) \(H_0\approx 55.8\)
  - \(\delta=+2.00\) \(\Rightarrow\) \(H_0\approx 44.3\)

MC intervals for each forecast point are recorded in:

- `math/hubble_env_scan/results_v1.json` (`predictions_mc_direct`)

---

## External consistency (direction / order-of-magnitude)

The repository includes an external N-body reference check (kept as a separate \(\delta\) convention):

- Gavas et al. (2024, \(\delta_{\mathrm{den}}\) within 10 Mpc/h) gives \(\beta \sim 0.10\) (trendline fit),
  which is directionally consistent while not directly comparable due to definition/scale differences.

---

## Reproducibility (3 commands)

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

---

## Falsification criteria (v0 convention)

Any of the following outcomes requires revising or abandoning the model:

- Multi-environment fits yield \(\beta \le 0\) with significance \(\ge 2\sigma\).
- Using a single explicitly stated \(\delta\) definition/scale, the monotonic ordering implied by the ratio law must hold:
  if \(\delta_A<\delta_B\), then \(H_{0,A}/H_{0,B} > 1\) (equivalently \(H_{0,A} > H_{0,B}\)).
  Systematic violations at \(\ge 2\sigma\) falsify the model.
- In overdense environments (\(\delta>0\)), measured local \(H_0\) is systematically **above** \(H_{0,\mathrm{ref}}\)
  at \(\ge 2\sigma\) after controlling for systematics.
- Residuals are strongly nonlinear and cannot be absorbed by a first-order correction (requiring a scale-dependent or nonlinear upgrade; if still failing, the path is falsified).

---

## Submission-ready bundle

An arXiv-ready TeX bundle (main TeX + required figures) can be built as:

```bash
python papers/build_arxiv_hubble_tension_environment.py
```

Output:

- `papers/_arxiv_hubble_tension_environment/hubble_tension_environment_arxiv_v0.zip`

English draft (for endorsers / arXiv):

- `papers/hubble_tension_environment_arxiv_v0.tex`
- `papers/hubble_tension_environment_v0_en.md`


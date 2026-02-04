# Environmental Correction to the Hubble Tension: A Reproducible, Falsifiable v0 Relational-Network Model

Author: Baiyi Wang  
Date: 2026-01-31  
Version: v0 (research note / preprint draft; not peer reviewed)

---

## Abstract

The Hubble tension refers to the persistent discrepancy between early-universe (CMB) and late-universe (distance-ladder / SNe) measurements of the Hubble constant \(H_0\). This note proposes a minimal, testable “observer-environment effect” model: **the locally inferred \(H_0\) depends systematically on the large-scale density contrast \(\delta\) around the observer**. In a linear parameterization,
\[
H_0(\delta)=H_{0,\mathrm{ref}}\,(1-\beta\,\delta),
\]
where \(H_{0,\mathrm{ref}}\) is a reference expansion rate (conventionally anchored to the CMB/Planck value; i.e., a normalization choice) and \(\beta>0\) is a density–expansion coupling coefficient. The model predicts higher local \(H_0\) in underdense regions and lower local \(H_0\) in overdense regions. This repository provides a reproducible v0/v1 pipeline (CSV dataset, fitting scripts, JSON outputs, and figures) and explicit falsification criteria. At present, external N-body simulations and selected observational results are **not inconsistent** with the predicted sign and order of magnitude, while the available point set and the definition/scale of \(\delta\) remain the dominant limitations. The goal of this note is therefore narrow: to freeze a simple hypothesis in a reproducible, falsifiable form, enabling future data to decisively support, refine, or refute it.

---

## 1. Problem statement

Typical representative values (for context):

- Planck 2018 (CMB): \(H_0 \approx 67.4\ \mathrm{km\,s^{-1}\,Mpc^{-1}}\)
- SH0ES (distance ladder): \(H_0 \approx 73.2\ \mathrm{km\,s^{-1}\,Mpc^{-1}}\)

The relative difference is \(\sim 9\%\), at multi-\(\sigma\) significance.

This note does **not** claim a definitive solution. Its objective is:

- to state a **minimal testable parameterization**;
- to provide a **reproducible pipeline**;
- to define **clear falsification criteria**.

---

## 2. Model

### 2.1 Definitions

- \(H_{0,\mathrm{ref}}\): reference expansion rate (conventionally anchored to CMB/Planck; normalization choice)
- \(\delta\): (observer) density contrast relative to the mean, \(\delta=(\rho-\bar\rho)/\bar\rho\) (definition depends on the source; see §3)
- \(\beta\): coupling coefficient to be estimated from multi-environment datasets

### 2.2 Main equation

\[
H_0(\delta)=H_{0,\mathrm{ref}}\,(1-\beta\,\delta).
\]

Interpretation:

- \(\delta<0\) (underdense) \(\Rightarrow H_0(\delta) > H_{0,\mathrm{ref}}\)
- \(\delta>0\) (overdense) \(\Rightarrow H_0(\delta) < H_{0,\mathrm{ref}}\)

### 2.3 Linear-theory baseline (comparison)

In standard linear theory one often obtains \(\delta H/H \sim -(f/3)\delta\), with \(f\simeq \Omega_m^{0.55}\). For Planck-like parameters \(\Omega_m\approx 0.315\),
\[
\beta_{\mathrm{lin}}\approx \frac{f}{3}\approx 0.177.
\]
This serves as an external order-of-magnitude baseline for \(\beta\).

### 2.4 A relational-operator toy derivation (operational)

This repository also contains a minimal toy chain (not a full cosmological normalization) to demonstrate the sign and near-linearity:

1. Use a 3D torus graph where the graph Laplacian \(L\approx -\nabla^2\) in the continuum limit.
2. Define a compensated density-contrast field \(\delta_{\mathrm{rel}}\) and solve a Poisson response:
   \[
   L\phi=-\delta_{\mathrm{rel}}.
   \]
3. Define an entropic/potential-gradient flow proxy:
   \[
   \mathbf v=-f\nabla\phi.
   \]
4. Estimate the regression bias \(\delta H/H\) using the same shell/window estimator as in the pipeline, yielding \(\beta \approx f/3\).

Script: `math/spectral_graph/run_relational_hubble_beta_poisson_3d_v0.py`  
Outputs: `relational_hubble_beta_poisson_3d_results_v0.json` / `..._plot_v0.png`

---

## 3. Data and “\(\delta\)” conventions (v0/v1)

The dominant practical limitation in v0 is that \(\delta\) values are taken from heterogeneous sources (different smoothing scales, definitions, bias corrections). The pipeline therefore:

- records the data source and notes;
- keeps “fit points” and “forecast points” separated;
- treats \(\delta\)-definition unification as a v1/v2 priority.

See:

- `math/hubble_env_scan/hubble_env_data_v1.csv`
- `math/hubble_env_scan/README.md`

---

## 4. Results (current reproducible outputs)

Running:

```bash
python math/hubble_env_scan/run_hubble_env_scan_v1.py
```

produces:

- `math/hubble_env_scan/results_v1.json`
- `math/hubble_env_scan/plot_v1.png`

The current v1 direct-fit (n=3) summary (see JSON for details):

- \(H_{0,\mathrm{ref}}=67.36\) (Planck/CMB anchor; convention)
- \(\beta \approx 0.171 \pm 0.024\)
- MC (16/50/84): 0.139 / 0.169 / 0.201

Forecast points (frozen in the pipeline) include, for example:

- Cold Spot-like underdensity (\(\delta\approx-0.6\sim-0.7\)): \(H_0 \approx 74\sim 75\)
- Typical cluster (\(\delta\approx+0.3\sim+0.5\)): \(H_0 \approx 62\sim 64\)

---

## 5. Falsification criteria (v0 convention)

Any of the following outcomes would require revising or abandoning the model:

- Multi-environment fits yield \(\beta \le 0\) with significance \(\ge 2\sigma\).
- In overdense environments (\(\delta>0\)), measured local \(H_0\) is systematically **above** \(H_{0,\mathrm{ref}}\) at \(\ge 2\sigma\) after controlling for systematics.
- Residuals are strongly nonlinear and cannot be absorbed by a first-order correction (necessitating a scale-dependent or nonlinear upgrade; if still failing, the path is falsified).

See: `docs/verification/predictions.md` and `docs/physics/hubble_tension.md`.

---

## 6. Reproducibility (minimal)

```bash
pip install -r requirements.txt
python math/hubble_env_scan/run_hubble_env_scan_v1.py
python math/spectral_graph/run_relational_hubble_beta_poisson_3d_v0.py
```

See also: `REPRODUCE_HUBBLE_TENSION.md`.


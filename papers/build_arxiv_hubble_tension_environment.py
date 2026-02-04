"""
Build a minimal arXiv upload bundle for the Hubble-tension environment-correction note.

Why this exists
---------------
The repository contains many large result files. Uploading the whole repository to arXiv is
unnecessary and can exceed size limits or complicate compilation.

This script assembles a small folder with:
  - LaTeX main file
  - Only the required figures (copied into fig/)
and then creates a ZIP file suitable for arXiv upload.

Usage (from repo root)
---------------------
python papers/build_arxiv_hubble_tension_environment.py

Outputs
-------
papers/_arxiv_hubble_tension_environment/
  - hubble_tension_environment_arxiv_v0.tex
  - fig/*.png
  - hubble_tension_environment_arxiv_v0.zip
"""

from __future__ import annotations

import shutil
import zipfile
from dataclasses import dataclass
from pathlib import Path


@dataclass(frozen=True)
class FigureSpec:
    src: Path
    dst_name: str


def _repo_root() -> Path:
    # This file lives at <repo>/papers/build_arxiv_hubble_tension_environment.py
    return Path(__file__).resolve().parents[1]


def main() -> None:
    root = _repo_root()

    tex_src = root / "papers" / "hubble_tension_environment_arxiv_v0.tex"

    out_dir = root / "papers" / "_arxiv_hubble_tension_environment"
    fig_dir = out_dir / "fig"
    out_zip = out_dir / "hubble_tension_environment_arxiv_v0.zip"

    figures = [
        FigureSpec(
            src=root / "math" / "hubble_env_scan" / "plot_v1.png",
            dst_name="hubble_env_scan_plot_v1.png",
        ),
        FigureSpec(
            src=root / "math" / "spectral_graph" / "relational_hubble_beta_poisson_3d_plot_v0.png",
            dst_name="relational_hubble_beta_poisson_3d_plot_v0.png",
        ),
    ]

    missing = [p for p in [tex_src] if not p.exists()]
    missing += [f.src for f in figures if not f.src.exists()]
    if missing:
        msg = "Missing required files:\n" + "\n".join(f"  - {p}" for p in missing)
        raise FileNotFoundError(msg)

    # Ensure the bundle is minimal and not polluted by stale files from prior runs.
    if out_dir.exists():
        shutil.rmtree(out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    fig_dir.mkdir(parents=True, exist_ok=True)

    shutil.copy2(tex_src, out_dir / "hubble_tension_environment_arxiv_v0.tex")

    for fig in figures:
        shutil.copy2(fig.src, fig_dir / fig.dst_name)

    # Build zip (flat folder structure, preserving fig/ subdir)
    if out_zip.exists():
        out_zip.unlink()

    with zipfile.ZipFile(out_zip, mode="w", compression=zipfile.ZIP_DEFLATED) as z:
        z.write(out_dir / "hubble_tension_environment_arxiv_v0.tex", arcname="hubble_tension_environment_arxiv_v0.tex")
        for fig in figures:
            z.write(fig_dir / fig.dst_name, arcname=f"fig/{fig.dst_name}")

    print(f"Wrote folder: {out_dir}")
    print(f"Wrote zip:    {out_zip}")


if __name__ == "__main__":
    main()


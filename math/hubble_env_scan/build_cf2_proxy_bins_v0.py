"""
Build proxy (delta_g*) environment points from Cosmicflows-2 (CF2) + 2M++ (v0)
=============================================================================

Why this exists
--------------
Our main-line "direct delta" dataset is still tiny (CMB + local void). To quickly
expand the (delta, H0) coverage without hand-collecting many papers, we can use:

- Cosmicflows-2 (Tully+ 2013) distances + CMB-frame velocities  -> H0 ~ V/D
- 2M++ reconstructed density field (Carrick+ 2015)              -> delta_g*

This produces a *proxy track* (delta_kind=delta_g_star_2mpp). It is not a
replacement for literature "direct delta" points, but it is an efficient way
to check:

- correlation direction (overdense -> lower inferred H0?)
- order-of-magnitude beta in H0(delta)=H0_global*(1-beta*delta)

Outputs
-------
1) A "bins-only" CSV in the hubble_env_data_v1 schema (appendable).
2) (Optional) A merged dataset file (base v1 + proxy bins).

Notes / caveats
---------------
- H0 is computed as Vcmba / Dist for individual objects; peculiar velocities
  can dominate at small distances. Use --dist-min / --vmin to mitigate this.
- H0_err is propagated from distance fractional error and an optional
  peculiar-velocity term (default 250 km/s).
- delta_g* is a luminosity-weighted galaxy overdensity proxy, not delta_m.
"""

from __future__ import annotations

import argparse
import csv
import math
import urllib.request
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Sequence, Tuple

import numpy as np

from twompp_density_lookup_v0 import lookup_delta_g_star


VIZIER_ASU_TSV = "https://vizier.cfa.harvard.edu/viz-bin/asu-tsv"
CF2_TABLE1 = "J/AJ/146/86/table1"  # galaxy distances
CF2_TABLE2 = "J/AJ/146/86/table2"  # group distances


@dataclass
class Cf2Row:
    glon_deg: float
    glat_deg: float
    dist_mpc: float
    vcmba_kms: float
    dist_frac_err: Optional[float] = None
    n1: Optional[int] = None
    n2: Optional[int] = None
    vel_disp_sigma_kms: Optional[float] = None
    pva_kms: Optional[float] = None


def _cache_dir() -> Path:
    d = Path(__file__).parent / "_cache"
    d.mkdir(parents=True, exist_ok=True)
    return d


def _download_if_missing(url: str, path: Path) -> None:
    if path.exists():
        return
    print(f"Downloading CF2 from VizieR -> {path}")
    print(f"Source: {url}")
    urllib.request.urlretrieve(url, path)  # noqa: S310 (intentional external download)


def _build_cf2_request_url(*, table: str, out_max: int) -> str:
    # Only request columns we actually need.
    # Column names follow VizieR table metadata.
    table = str(table).strip()
    if table.endswith("/table2"):
        cols = "N1,<Dist>,N2,GLON,GLAT,<Vcmba>,<PVa>,sigma,Group,Ng"
    else:
        cols = "Dist,Err,Vcmba,GLON,GLAT"
    return f"{VIZIER_ASU_TSV}?-source={table}&-out.max={int(out_max)}&-out={cols}"


def _read_cf2_asu_tsv(path: Path, *, table: str, group_base_dist_frac_err: float) -> List[Cf2Row]:
    """
    Parse a VizieR asu-tsv response.

    The format includes many '#' metadata lines, then a tab-separated table with:
      header row (column names),
      units row,
      separator row,
      data rows.
    """
    rows: List[Cf2Row] = []
    with path.open("r", encoding="utf-8", errors="replace") as f:
        # Skip metadata lines starting with '#'
        lines = [ln.rstrip("\n") for ln in f if not ln.startswith("#")]

    # Find the header line that starts the TSV table
    header_idx = None
    for i, ln in enumerate(lines):
        if "\t" in ln and "GLON" in ln and "GLAT" in ln:
            header_idx = i
            break
    if header_idx is None:
        raise ValueError("Could not find TSV header in VizieR response.")

    # The table starts at header_idx; data begins after 3 lines (header, units, dashed)
    table_lines = lines[header_idx:]
    if len(table_lines) < 4:
        return rows

    tsv = "\n".join(table_lines)
    reader = csv.DictReader(tsv.splitlines(), delimiter="\t")
    # The DictReader will treat the dashed separator as a data row; we skip non-numeric Dist.
    table = str(table).strip()
    is_group = table.endswith("/table2")

    for r in reader:
        try:
            glon = float((r.get("GLON") or "").strip())
            glat = float((r.get("GLAT") or "").strip())
        except Exception:
            continue
        if not (math.isfinite(glon) and math.isfinite(glat)):
            continue

        if is_group:
            # Group table uses bracketed column names.
            try:
                dist = float((r.get("<Dist>") or "").strip())
                v = float((r.get("<Vcmba>") or "").strip())
            except Exception:
                continue
            if not (math.isfinite(dist) and math.isfinite(v)) or dist <= 0:
                continue

            # Approximate fractional distance error: typical CF2 per-galaxy err / sqrt(N1).
            n1 = None
            try:
                n1 = int(float((r.get("N1") or "").strip()))
            except Exception:
                n1 = None
            dist_frac_err = None
            if n1 is not None and n1 > 0:
                dist_frac_err = float(group_base_dist_frac_err) / math.sqrt(float(n1))

            n2 = None
            try:
                n2 = int(float((r.get("N2") or "").strip()))
            except Exception:
                n2 = None

            sigma = None
            try:
                sigma = float((r.get("sigma") or "").strip())
            except Exception:
                sigma = None

            pva = None
            try:
                pva = float((r.get("<PVa>") or "").strip())
            except Exception:
                pva = None

            rows.append(
                Cf2Row(
                    glon_deg=glon,
                    glat_deg=glat,
                    dist_mpc=dist,
                    vcmba_kms=v,
                    dist_frac_err=dist_frac_err,
                    n1=n1,
                    n2=n2,
                    vel_disp_sigma_kms=sigma,
                    pva_kms=pva,
                )
            )
        else:
            try:
                dist = float((r.get("Dist") or "").strip())
                err = float((r.get("Err") or "").strip())
                v = float((r.get("Vcmba") or "").strip())
            except Exception:
                continue
            if not (math.isfinite(dist) and math.isfinite(err) and math.isfinite(v)) or dist <= 0:
                continue
            rows.append(Cf2Row(glon_deg=glon, glat_deg=glat, dist_mpc=dist, vcmba_kms=v, dist_frac_err=err))
    return rows


def _weighted_mean(x: np.ndarray, w: np.ndarray) -> float:
    s = float(np.sum(w))
    if not math.isfinite(s) or s <= 0:
        return float("nan")
    return float(np.sum(w * x) / s)


def _weighted_se_mean(w: np.ndarray) -> float:
    # Standard error of weighted mean under independent errors: sqrt(1/sum(w))
    s = float(np.sum(w))
    if not math.isfinite(s) or s <= 0:
        return float("nan")
    return float(math.sqrt(1.0 / s))


def _assign_category(delta: float) -> str:
    if delta <= -0.10:
        return "void"
    if delta >= 0.10:
        return "cluster"
    return "supercluster"


def _write_v1_rows(path: Path, rows: List[Dict[str, object]]) -> None:
    fieldnames = [
        "env_id",
        "env_name",
        "category",
        "data_kind",
        "use_for_fit",
        "is_global_anchor",
        "delta",
        "delta_err",
        "H0",
        "H0_err",
        "redshift_range",
        "scale_mpc",
        "source",
        "notes",
        "delta_kind",
        "delta_scale_mpc_h",
        "delta_source",
    ]
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=fieldnames)
        w.writeheader()
        for r in rows:
            w.writerow({k: r.get(k) for k in fieldnames})


def _merge_v1_dataset(base_path: Path, extra_rows: List[Dict[str, object]], out_path: Path) -> None:
    with base_path.open("r", newline="", encoding="utf-8") as f:
        base = list(csv.DictReader(f))
    if not base:
        raise ValueError(f"Base dataset is empty: {base_path}")

    # Preserve column order from base file, but ensure v1 fields exist.
    base_fieldnames = list(base[0].keys())
    for k in [
        "env_id",
        "env_name",
        "category",
        "data_kind",
        "use_for_fit",
        "is_global_anchor",
        "delta",
        "delta_err",
        "H0",
        "H0_err",
        "redshift_range",
        "scale_mpc",
        "source",
        "notes",
        "delta_kind",
        "delta_scale_mpc_h",
        "delta_source",
    ]:
        if k not in base_fieldnames:
            base_fieldnames.append(k)

    merged = base + [{k: r.get(k, "") for k in base_fieldnames} for r in extra_rows]
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as f:
        w = csv.DictWriter(f, fieldnames=base_fieldnames)
        w.writeheader()
        for r in merged:
            w.writerow({k: r.get(k, "") for k in base_fieldnames})


def main() -> None:
    ap = argparse.ArgumentParser(description="Build CF2+2M++ proxy bins (v0).")
    ap.add_argument(
        "--table",
        type=str,
        default="table2",
        choices=["table1", "table2"],
        help="Which CF2 table to use: table1 (galaxies) or table2 (groups). Default: table2.",
    )
    ap.add_argument("--out-max", type=int, default=12000, help="VizieR row cap (CF2 has 8315 rows)")
    ap.add_argument(
        "--force-download",
        action="store_true",
        help="Redownload CF2 even if a cached file exists.",
    )
    ap.add_argument("--smooth", type=float, default=50.0, help="Total Gaussian smoothing sigma for 2M++ [Mpc/h] (>=4)")
    ap.add_argument("--dist-min", type=float, default=40.0, help="Min CF2 distance [Mpc] (reduce peculiar-velocity dominance)")
    ap.add_argument("--dist-max", type=float, default=250.0, help="Max CF2 distance [Mpc]")
    ap.add_argument("--err-max", type=float, default=0.20, help="Max CF2 fractional distance error")
    ap.add_argument("--vmin", type=float, default=3000.0, help="Min Vcmba [km/s]")
    ap.add_argument("--vmax", type=float, default=20000.0, help="Max Vcmba [km/s] (2M++ cube limit)")
    ap.add_argument("--peculiar-vel-err", type=float, default=250.0, help="Additive velocity error [km/s] for H0_err propagation")
    ap.add_argument(
        "--velocity-mode",
        type=str,
        default="raw",
        choices=["raw", "hubble_flow"],
        help="Velocity used for H0=V/Dist. raw: use Vcmba. hubble_flow: use (Vcmba - PVa) (table2 only).",
    )
    ap.add_argument(
        "--h-ref",
        type=float,
        default=0.6736,
        help="Reference h=H0/100 used to convert CF2 Dist[Mpc] -> r[Mpc/h] for 2M++ lookup (default: Planck 67.36/100).",
    )
    ap.add_argument(
        "--r-from",
        type=str,
        default="dist",
        choices=["dist", "vcmba"],
        help="How to set r (Mpc/h) for the 2M++ lookup: dist (Dist*h_ref) or vcmba (Vcmba/100).",
    )
    ap.add_argument(
        "--group-base-dist-frac-err",
        type=float,
        default=0.20,
        help="(table2 only) Base fractional distance error for a single galaxy; group error ~ base/sqrt(N1).",
    )
    ap.add_argument(
        "--group-min-n1",
        type=int,
        default=2,
        help="(table2 only) Minimum N1 (count of distance galaxies) to keep a group.",
    )
    ap.add_argument("--n-bins", type=int, default=7, help="Number of delta bins (quantile bins)")
    ap.add_argument("--min-per-bin", type=int, default=200, help="Minimum count per bin")
    ap.add_argument(
        "--base-data-file",
        type=str,
        default=str(Path(__file__).parent / "hubble_env_data_v1.csv"),
        help="Base v1 dataset to merge with (optional)",
    )
    ap.add_argument(
        "--out-bins-file",
        type=str,
        default=str(Path(__file__).parent / "hubble_env_cf2_proxy_bins_v0.csv"),
        help="Output CSV (bins only) in hubble_env_data_v1 schema",
    )
    ap.add_argument(
        "--out-merged-file",
        type=str,
        default="",
        help="Optional output merged dataset (base + bins). If empty, do not write merged file.",
    )
    args = ap.parse_args()

    # Cache filename depends on out_max to avoid reusing a smaller preview download.
    table = CF2_TABLE2 if str(args.table).strip() == "table2" else CF2_TABLE1
    cache_path = _cache_dir() / f"cf2_{Path(table).name}_asu_tsv_outmax_{int(args.out_max)}.tsv"
    url = _build_cf2_request_url(table=table, out_max=int(args.out_max))
    if bool(args.force_download) and cache_path.exists():
        cache_path.unlink()
    _download_if_missing(url, cache_path)

    raw = _read_cf2_asu_tsv(cache_path, table=table, group_base_dist_frac_err=float(args.group_base_dist_frac_err))
    print(f"Loaded CF2 rows: {len(raw)}")

    # Compute per-object proxy values
    deltas: List[float] = []
    h0s: List[float] = []
    h0_errs: List[float] = []

    for r in raw:
        if r.dist_mpc < float(args.dist_min) or r.dist_mpc > float(args.dist_max):
            continue
        v_for_selection = float(r.vcmba_kms)
        if str(args.velocity_mode) == "hubble_flow" and table.endswith("/table2") and r.pva_kms is not None:
            v_for_selection = float(r.vcmba_kms) - float(r.pva_kms)
        if v_for_selection < float(args.vmin) or v_for_selection > float(args.vmax):
            continue
        if table.endswith("/table2"):
            if r.n1 is None or int(r.n1) < int(args.group_min_n1):
                continue
        if r.dist_frac_err is not None:
            if r.dist_frac_err <= 0 or r.dist_frac_err > float(args.err_max):
                continue

        if str(args.r_from) == "vcmba":
            # Map a velocity scale to a comoving distance scale in Mpc/h: r = V/100
            r_mpc_h = float(r.vcmba_kms) / 100.0
        else:
            # Decouple the 2M++ lookup position from the derived H0=V/D.
            r_mpc_h = float(r.dist_mpc) * float(args.h_ref)
        if r_mpc_h <= 0 or r_mpc_h > 200.0:
            continue

        try:
            lookup = lookup_delta_g_star(
                l_deg=float(r.glon_deg),
                b_deg=float(r.glat_deg),
                r_mpc_h=float(r_mpc_h),
                smoothing_gaussian_mpc_h=float(args.smooth),
            )
            delta_g = float(lookup["delta_g_star"])
        except Exception:
            continue
        if not math.isfinite(delta_g):
            continue

        v_for_h0 = float(r.vcmba_kms)
        if str(args.velocity_mode) == "hubble_flow" and table.endswith("/table2") and r.pva_kms is not None:
            v_for_h0 = float(r.vcmba_kms) - float(r.pva_kms)
        if not math.isfinite(v_for_h0) or v_for_h0 <= 0:
            continue

        H0 = float(v_for_h0) / float(r.dist_mpc)
        # Propagate distance fractional error, plus an optional velocity term.
        sigma_v = float(args.peculiar_vel_err)
        if r.vel_disp_sigma_kms is not None and math.isfinite(float(r.vel_disp_sigma_kms)):
            sigma_v = math.sqrt(sigma_v * sigma_v + float(r.vel_disp_sigma_kms) * float(r.vel_disp_sigma_kms))
        dist_frac_err = float(r.dist_frac_err) if r.dist_frac_err is not None else float(args.err_max)
        H0_err = math.sqrt((H0 * dist_frac_err) ** 2 + (sigma_v / float(r.dist_mpc)) ** 2)

        if not (math.isfinite(H0) and math.isfinite(H0_err) and H0_err > 0):
            continue

        deltas.append(delta_g)
        h0s.append(H0)
        h0_errs.append(H0_err)

    if len(deltas) < int(args.n_bins) * int(args.min_per_bin):
        print(
            f"Warning: filtered rows={len(deltas)} may be too small for stable binning "
            f"(n_bins={int(args.n_bins)}, min_per_bin={int(args.min_per_bin)})."
        )
    else:
        print(f"Filtered rows: {len(deltas)}")

    deltas_np = np.asarray(deltas, dtype=float)
    h0_np = np.asarray(h0s, dtype=float)
    h0e_np = np.asarray(h0_errs, dtype=float)

    # Quantile bins to ensure enough points in each bin.
    qs = np.linspace(0.0, 1.0, int(args.n_bins) + 1)
    edges = np.quantile(deltas_np, qs)
    # Ensure strictly increasing edges (handle ties)
    edges = np.unique(edges)
    if edges.size < 4:
        raise RuntimeError("Not enough delta variation in filtered data to build bins.")

    bin_rows: List[Dict[str, object]] = []
    total_used = 0
    for bi in range(edges.size - 1):
        lo = float(edges[bi])
        hi = float(edges[bi + 1])
        if bi == edges.size - 2:
            m = (deltas_np >= lo) & (deltas_np <= hi)
        else:
            m = (deltas_np >= lo) & (deltas_np < hi)
        idx = np.where(m)[0]
        n = int(idx.size)
        if n < int(args.min_per_bin):
            continue
        total_used += n

        d = deltas_np[idx]
        y = h0_np[idx]
        ye = h0e_np[idx]
        w = 1.0 / (ye * ye)

        d_mean = _weighted_mean(d, w)
        y_mean = _weighted_mean(y, w)
        y_se = _weighted_se_mean(w)

        if not (math.isfinite(d_mean) and math.isfinite(y_mean) and math.isfinite(y_se)):
            continue

        env_id = f"cf2_proxy_bin_{bi+1:02d}_sigma{float(args.smooth):.0f}"
        env_name = f"CF2 proxy bin {bi+1} (δ_g*∈[{lo:.3f},{hi:.3f}], n={n})"
        bin_rows.append(
            {
                "env_id": env_id,
                "env_name": env_name,
                "category": _assign_category(float(d_mean)),
                "data_kind": "observed",
                "use_for_fit": True,
                "is_global_anchor": False,
                "delta": float(d_mean),
                "delta_err": None,
                "H0": float(y_mean),
                "H0_err": float(y_se),
                "redshift_range": "local",
                "scale_mpc": f"bin from CF2; smooth={float(args.smooth):.0f} Mpc/h",
                "source": "Cosmicflows-2 (Tully et al. 2013; VizieR J/AJ/146/86) + 2M++ (Carrick et al. 2015)",
                "notes": f"Quantile-binned mean; H0=Vcmba/Dist with peculiar_vel_err={float(args.peculiar_vel_err):.0f} km/s",
                "delta_kind": "delta_g_star_2mpp",
                "delta_scale_mpc_h": float(args.smooth),
                "delta_source": str(cache_path),
            }
        )

    print(f"Built bins: {len(bin_rows)} (using {total_used} objects)")

    out_bins = Path(args.out_bins_file)
    _write_v1_rows(out_bins, bin_rows)
    print(f"Wrote bins-only v1 rows: {out_bins}")

    if str(args.out_merged_file).strip():
        base_path = Path(args.base_data_file)
        out_merged = Path(args.out_merged_file)
        _merge_v1_dataset(base_path, bin_rows, out_merged)
        print(f"Wrote merged dataset: {out_merged}")


if __name__ == "__main__":
    main()


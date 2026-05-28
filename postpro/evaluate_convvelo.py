import argparse
from pathlib import Path

import dask
import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar

try:
    from build_dataset import build_convvelo_dataset
    from derive import fourier_dy_fast, get_banded_derivative_operators
except ImportError:
    from .build_dataset import build_convvelo_dataset
    from .derive import fourier_dy_fast, get_banded_derivative_operators

def _uses_full_reconstruction(cvv: xr.Dataset) -> bool:
    return "v_cross_u" in cvv


def _with_mean_derivatives(cvv: xr.Dataset) -> xr.Dataset:
    result = cvv.copy()
    d0_banded, d1_banded, _ = get_banded_derivative_operators(result.attrs, result.y)
    use_full_reconstruction = _uses_full_reconstruction(result)

    result["u_mean"] = result["mean_u"]
    result["dudy_mean"] = fourier_dy_fast(result["u_mean"].chunk({"y": -1}), d1_banded, d0_banded)

    if use_full_reconstruction:
        if "mean_v" not in result or "mean_w" not in result:
            raise ValueError("Full convvelo reconstruction requires mean_v and mean_w profiles.")
        result["v_mean"] = result["mean_v"]
        result["w_mean"] = result["mean_w"]
        result["dvdy_mean"] = fourier_dy_fast(result["v_mean"].chunk({"y": -1}), d1_banded, d0_banded)
        result["dwdy_mean"] = fourier_dy_fast(result["w_mean"].chunk({"y": -1}), d1_banded, d0_banded)

    if "mean_t" in result:
        result["t_mean"] = result["mean_t"]
        result["dtdy_mean"] = fourier_dy_fast(result["t_mean"].chunk({"y": -1}), d1_banded, d0_banded)

    return result


def reconstruct_parts(cvv: xr.Dataset) -> xr.Dataset:
    use_full_reconstruction = _uses_full_reconstruction(cvv)
    result = _with_mean_derivatives(cvv)
    result.attrs["convvelo_reconstruction_mode"] = "full" if use_full_reconstruction else "minimal"

    result["u_p1"] = result["u_cross_p"]
    result["u_p2"] = result["u_cross_uu"] - 2.0 * result["u_mean"] * result["u_cross_u"]
    result["u_p3"] = result["u_cross_uw"] - result["u_mean"] * result["u_cross_w"] - (
        result["w_mean"] * result["u_cross_u"] if use_full_reconstruction else 0.0
    )
    result["u_p4"] = result["u_cross_v"]
    result["u_p5"] = result["u_cross_dyyu"]
    result["u_p6"] = (
        result["u_cross_dyuv"]
        - result["dudy_mean"] * result["u_cross_v"]
        - result["u_mean"] * result["u_cross_dyv"]
        - (result["dvdy_mean"] * result["u_cross_u"] if use_full_reconstruction else 0.0)
        - (result["v_mean"] * result["u_cross_dyu"] if use_full_reconstruction else 0.0)
    )
    result["u_p7"] = result["u_cross_u"]

    result["v_p1"] = result["v_cross_uv"] - result["u_mean"] * result["v_cross_v"] - (
        result["v_mean"] * result["v_cross_u"] if use_full_reconstruction else 0.0
    )
    result["v_p2"] = result["v_cross_vw"] - (
        result["v_mean"] * result["v_cross_w"] + result["w_mean"] * result["v_cross_v"]
        if use_full_reconstruction
        else 0.0
    )
    result["v_p3"] = result["v_cross_dpdy"]
    result["v_p4"] = result["v_cross_dyyv"]
    result["v_p5"] = result["v_cross_dyvv"] - (
        2.0 * result["dvdy_mean"] * result["v_cross_v"] + 2.0 * result["v_mean"] * result["v_cross_dyv"]
        if use_full_reconstruction
        else 0.0
    )
    result["v_p6"] = result["v_cross_v"]

    result["w_p1"] = result["w_cross_uw"] - result["u_mean"] * result["w_cross_w"] - (
        result["w_mean"] * result["w_cross_u"] if use_full_reconstruction else 0.0
    )
    result["w_p2"] = result["w_cross_p"]
    result["w_p3"] = result["w_cross_ww"] - (
        2.0 * result["w_mean"] * result["w_cross_w"] if use_full_reconstruction else 0.0
    )
    result["w_p4"] = result["w_cross_dyyw"]
    result["w_p5"] = result["w_cross_dyvw"] - (
        result["dvdy_mean"] * result["w_cross_w"]
        + result["v_mean"] * result["w_cross_dyw"]
        + result["dwdy_mean"] * result["w_cross_v"]
        + result["w_mean"] * result["w_cross_dyv"]
        if use_full_reconstruction
        else 0.0
    )
    result["w_p6"] = result["w_cross_w"]

    if "mean_t" in result:
        result["t_p1"] = result["t_cross_tu"] - result["u_mean"] * result["t_cross_t"] - result["t_mean"] * result["t_cross_u"]
        if "t_cross_tw" in result:
            result["t_p2"] = result["t_cross_tw"] - result["t_mean"] * result["t_cross_w"] - (
                result["w_mean"] * result["t_cross_t"] if use_full_reconstruction else 0.0
            )
        result["t_p3"] = result["t_cross_v"]
        result["t_p4"] = result["t_cross_dyyt"]
        result["t_p5"] = (
            result["t_cross_dytv"]
            - result["dtdy_mean"] * result["t_cross_v"]
            - result["t_mean"] * result["t_cross_dyv"]
            - (result["dvdy_mean"] * result["t_cross_t"] if use_full_reconstruction and "t_cross_dyt" in result else 0.0)
            - (result["v_mean"] * result["t_cross_dyt"] if use_full_reconstruction and "t_cross_dyt" in result else 0.0)
        )
        result["t_p6"] = result["t_cross_t"]
        result["t_forcing"] = xr.zeros_like(result["t_p6"])

    scale = 1.0 / (result.attrs["alfa0"] * result.attrs["beta0"])
    result["energy_u"] = (scale * result["u_cross_u"]).real
    result["energy_v"] = (scale * result["v_cross_v"]).real
    result["energy_w"] = (scale * result["w_cross_w"]).real
    if "t_cross_t" in result:
        result["energy_t"] = (scale * result["t_cross_t"]).real

    return result


def contrib_u(cvv: xr.Dataset, contrib: xr.Dataset | None = None) -> xr.Dataset:
    if contrib is None:
        contrib = cvv.copy()

    contrib["u_c0"] = cvv["u_mean"].real
    contrib["u_c1"] = cvv.kx_folded * cvv["u_p1"].real
    contrib["u_c2"] = cvv.kx_folded * cvv["u_p2"].real
    contrib["u_c3"] = cvv.kz * cvv["u_p3"].real
    contrib["u_c4"] = cvv["dudy_mean"].real * cvv["u_p4"].imag
    contrib["u_c5"] = -cvv.attrs["ni"] * cvv["u_p5"].imag
    contrib["u_c6"] = cvv["u_p6"].imag
    contrib["u_c7"] = cvv.kx_folded * cvv["u_p7"].real
    contrib["u_uc"] = (
        (contrib["u_c1"] + contrib["u_c2"] + contrib["u_c3"] + contrib["u_c4"] + contrib["u_c5"] + contrib["u_c6"])
        / contrib["u_c7"]
        + contrib["u_c0"]
    )
    return contrib


def contrib_v(cvv: xr.Dataset, contrib: xr.Dataset | None = None) -> xr.Dataset:
    if contrib is None:
        contrib = cvv.copy()

    contrib["v_c0"] = cvv["u_mean"].real
    contrib["v_c1"] = cvv.kx_folded * cvv["v_p1"].real
    contrib["v_c2"] = cvv.kz * cvv["v_p2"].real
    contrib["v_c3"] = cvv["v_p3"].imag
    contrib["v_c4"] = -cvv.attrs["ni"] * cvv["v_p4"].imag
    contrib["v_c5"] = cvv["v_p5"].imag
    contrib["v_c6"] = cvv.kx_folded * cvv["v_p6"].real
    contrib["v_uc"] = (
        (contrib["v_c1"] + contrib["v_c2"] + contrib["v_c3"] + contrib["v_c4"] + contrib["v_c5"])
        / contrib["v_c6"]
        + contrib["v_c0"]
    )
    return contrib


def contrib_w(cvv: xr.Dataset, contrib: xr.Dataset | None = None) -> xr.Dataset:
    if contrib is None:
        contrib = cvv.copy()

    contrib["w_c0"] = cvv["u_mean"].real
    contrib["w_c1"] = cvv.kx_folded * cvv["w_p1"].real
    contrib["w_c2"] = cvv.kz * cvv["w_p2"].real
    contrib["w_c3"] = cvv.kz * cvv["w_p3"].real
    contrib["w_c4"] = -cvv.attrs["ni"] * cvv["w_p4"].imag
    contrib["w_c5"] = cvv["w_p5"].imag
    contrib["w_c6"] = cvv.kx_folded * cvv["w_p6"].real
    contrib["w_uc"] = (
        (contrib["w_c1"] + contrib["w_c2"] + contrib["w_c3"] + contrib["w_c4"] + contrib["w_c5"])
        / contrib["w_c6"]
        + contrib["w_c0"]
    )
    return contrib


def contrib_t(cvv: xr.Dataset, contrib: xr.Dataset | None = None) -> xr.Dataset:
    if contrib is None:
        contrib = cvv.copy()

    contrib["t_c0"] = cvv["u_mean"].real
    contrib["t_c1"] = cvv.kx_folded * cvv["t_p1"].real
    contrib["t_c2"] = cvv.kz * cvv["t_p2"].real
    contrib["t_c3"] = cvv["dtdy_mean"].real * cvv["t_p3"].imag
    contrib["t_c4"] = -cvv.attrs["ni"] / cvv["Pr"] * cvv["t_p4"].imag
    contrib["t_c5"] = cvv["t_p5"].imag
    contrib["t_c6"] = cvv.kx_folded * cvv["t_p6"].real
    contrib["forcing"] = -cvv["t_forcing"].imag
    contrib["t_uc"] = (
        (contrib["t_c1"] + contrib["t_c2"] + contrib["t_c3"] + contrib["t_c4"] + contrib["t_c5"] + contrib["forcing"])
        / contrib["t_c6"]
        + contrib["t_c0"]
    )
    return contrib


def evaluate_convvelo(cvv: xr.Dataset) -> xr.Dataset:
    result = reconstruct_parts(cvv)
    result = contrib_u(result, result)
    result = contrib_v(result, result)
    result = contrib_w(result, result)
    if "mean_t" in result:
        result = contrib_t(result, result)
    return result


def _resolve_zarr_path(path: Path) -> Path:
    return path if path.suffix == ".zarr" else path / "convvelo.zarr"


def _ensure_convvelo_store(path: Path, *, chunks_x: int = -1) -> xr.Dataset:
    if path.exists():
        return xr.open_zarr(path, consolidated=False, chunks={})

    data = build_convvelo_dataset(path.parent, chunks_x=chunks_x)
    data.to_zarr(path, mode="w", consolidated=False)
    return xr.open_zarr(path, consolidated=False)


def _evaluated_uc_fields(cvv: xr.Dataset) -> xr.Dataset:
    result = evaluate_convvelo(cvv)
    field_names = ["u_uc", "v_uc", "w_uc"]
    if "t_uc" in result:
        field_names.append("t_uc")
    return result[field_names]


def _compute_with_progress(delayed_obj, *, description: str = "Computing"):
    try:
        from dask.distributed import get_client, progress

        client = get_client()
    except (ImportError, ValueError):
        print(f"{description} with local progress bar...")
        with ProgressBar():
            return dask.compute(delayed_obj)[0]

    print(f"{description} with distributed progress...")
    future = client.compute(delayed_obj)
    progress(future)
    return future.result()



if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Evaluate convvelo and append only the reconstructed uc fields to convvelo.zarr."
    )
    parser.add_argument("directory", type=Path, help="Case directory containing convvelo.zarr or the zarr store itself")
    parser.add_argument("--chunks-x", type=int, default=-1, help="Chunk size for kx_folded if convvelo.zarr must be built first")
    args = parser.parse_args()

    from dask.distributed import Client, LocalCluster
    cluster = LocalCluster(n_workers=8, threads_per_worker=2)
    client = Client(cluster)

    zarr_path = _resolve_zarr_path(args.directory)
    convvelo = _ensure_convvelo_store(zarr_path, chunks_x=args.chunks_x)
    result = _evaluated_uc_fields(convvelo)
    print(result)
    write = convvelo.assign(result).to_zarr(zarr_path, mode="a", consolidated=False, compute=False)
    _compute_with_progress(write, description=f"Appending evaluated convvelo fields to {zarr_path}")
    print(f"Appended evaluated convvelo fields to {zarr_path}")

    client.close()
    cluster.close()

"""Export time-averaged raw statistics as headerless Fourier fields.

The output is arranged for later MPI-IO reading from Fortran: velocity
statistics precede scalar statistics, all statistics for scalar 1 precede all
statistics for scalar 2, and each individual complex field has wall-normal
index ``y`` contiguous in the file.
"""

import argparse
from pathlib import Path
import sys
from typing import BinaryIO, Dict, Any

REPO_ROOT = Path(__file__).resolve().parents[1]
POSTPRO_DIR = REPO_ROOT / "postpro"
for path in (REPO_ROOT / "tools", POSTPRO_DIR):
    path_text = str(path)
    if path_text not in sys.path:
        sys.path.insert(0, path_text)

import numpy as np
import xarray as xr
from dask.diagnostics import ProgressBar


from ffts import physical, fourier
from derive import fourier_dy_fast, get_banded_derivative_operators, setup_derivatives
from build_dataset import build_convvelo_dataset
import logging


OUTPUT_DTYPE = np.dtype("<c16")

logging.getLogger("distributed.shuffle._scheduler_plugin").setLevel(logging.ERROR)


SCALAR_PART_KEYS = [f"t_p{i}" for i in range(1, 7)]
SCALAR_RAW_STATISTIC_KEYS = [
    "t_cross_t",
    "t_cross_u",
    "t_cross_v",
    "t_cross_w",
    "t_cross_tu",
    "t_cross_tw",
    "t_cross_dyyt",
    "t_cross_dytv",
    "t_cross_dyt",
    "t_cross_dyv",
]
VELOCITY_RAW_STATISTIC_KEYS = [
    "u_cross_u",
    "u_cross_dyu",
    "u_cross_v",
    "u_cross_dyv",
    "u_cross_w",
    "u_cross_dyw",
    "u_cross_dyyu",
    "v_cross_u",
    "v_cross_dyu",
    "v_cross_v",
    "v_cross_dyv",
    "v_cross_w",
    "v_cross_dyw",
    "v_cross_dyyv",
    "w_cross_u",
    "w_cross_dyu",
    "w_cross_v",
    "w_cross_dyv",
    "w_cross_w",
    "w_cross_dyw",
    "w_cross_dyyw",
    "u_cross_p",
    "v_cross_dpdy",
    "w_cross_p",
    "u_cross_uu",
    "u_cross_uw",
    "v_cross_uv",
    "v_cross_vw",
    "w_cross_uw",
    "w_cross_ww",
    "u_cross_dyuv",
    "v_cross_dyvv",
    "w_cross_dyvw",
]

def compute_raw_statistics(
    raw_data: xr.Dataset,
    config: Dict[str, Any],
    D0_banded: xr.DataArray,
    D1_banded: xr.DataArray,
    D2_banded: xr.DataArray,
    include_velocity: bool = False,
) -> xr.Dataset:
    """Build raw co-spectra needed for scalar and, optionally, velocity parts."""

    if "t" not in raw_data:
        raise ValueError("Scalar raw statistics require a `t` field.")

    chunk_config = {"y": -1, "kz": 10}
    phys = {key: physical(raw_data[key], config) for key in "tuvw"}
    theta_star = raw_data["t"].conj()
    theta_star_yslab = theta_star.chunk(chunk_config)

    theta_flux = {key: fourier(phys["t"] * phys[key], config) for key in "uvw"}
    dtheta_dyy = fourier_dy_fast(raw_data["t"].chunk(chunk_config), D2_banded, D0_banded)
    dtheta_dy = fourier_dy_fast(raw_data["t"].chunk(chunk_config), D1_banded, D0_banded)
    dvdy = fourier_dy_fast(raw_data["v"].chunk(chunk_config), D1_banded, D0_banded)
    dthetav_dy = fourier_dy_fast(theta_flux["v"].chunk(chunk_config), D1_banded, D0_banded)

    stats = xr.Dataset(coords=raw_data.coords)
    stats["t_cross_t"] = theta_star * raw_data["t"]
    stats["t_cross_u"] = theta_star * raw_data["u"]
    stats["t_cross_v"] = theta_star * raw_data["v"]
    stats["t_cross_w"] = theta_star * raw_data["w"]
    stats["t_cross_tu"] = theta_star * theta_flux["u"]
    stats["t_cross_tw"] = theta_star * theta_flux["w"]
    stats["t_cross_dyyt"] = theta_star_yslab * dtheta_dyy
    stats["t_cross_dytv"] = theta_star_yslab * dthetav_dy
    stats["t_cross_dyt"] = theta_star_yslab * dtheta_dy
    stats["t_cross_dyv"] = theta_star_yslab * dvdy

    if include_velocity:
        first_dy = {
            "u": fourier_dy_fast(raw_data["u"].chunk(chunk_config), D1_banded, D0_banded),
            "v": dvdy,
            "w": fourier_dy_fast(raw_data["w"].chunk(chunk_config), D1_banded, D0_banded),
        }
        second_dy = {
            key: fourier_dy_fast(raw_data[key].chunk(chunk_config), D2_banded, D0_banded)
            for key in "uvw"
        }
        rst = {
            key: fourier(phys[key[0]] * phys[key[1]], config)
            for key in ["uu", "uv", "uw", "vv", "vw", "ww"]
        }
        drstdy = {
            key: fourier_dy_fast(rst[key].chunk(chunk_config), D1_banded, D0_banded)
            for key in ["uv", "vv", "vw"]
        }

        for outer in "uvw":
            outer_star = raw_data[outer].conj()
            outer_star_yslab = outer_star.chunk(chunk_config)
            for inner in "uvw":
                stats[f"{outer}_cross_{inner}"] = outer_star * raw_data[inner]
                stats[f"{outer}_cross_dy{inner}"] = outer_star_yslab * first_dy[inner]
            stats[f"{outer}_cross_dyy{outer}"] = outer_star_yslab * second_dy[outer]

        stats["u_cross_p"] = raw_data["u"].conj() * raw_data["p"]
        stats["v_cross_dpdy"] = raw_data["v"].conj() * raw_data["dpdy"]
        stats["w_cross_p"] = raw_data["w"].conj() * raw_data["p"]
        stats["u_cross_uu"] = raw_data["u"].conj() * rst["uu"]
        stats["u_cross_uw"] = raw_data["u"].conj() * rst["uw"]
        stats["v_cross_uv"] = raw_data["v"].conj() * rst["uv"]
        stats["v_cross_vw"] = raw_data["v"].conj() * rst["vw"]
        stats["w_cross_uw"] = raw_data["w"].conj() * rst["uw"]
        stats["w_cross_ww"] = raw_data["w"].conj() * rst["ww"]
        stats["u_cross_dyuv"] = raw_data["u"].conj().chunk(chunk_config) * drstdy["uv"]
        stats["v_cross_dyvv"] = raw_data["v"].conj().chunk(chunk_config) * drstdy["vv"]
        stats["w_cross_dyvw"] = raw_data["w"].conj().chunk(chunk_config) * drstdy["vw"]

    return stats


def _select_indices(indices: list[int] | None, size: int, label: str) -> list[int]:
    """Normalize optional repeated CLI indices and retain the requested order."""

    selected = list(range(size)) if indices is None else indices
    normalized = []
    for index in selected:
        value = index + size if index < 0 else index
        if value < 0 or value >= size:
            raise ValueError(f"{label} index {index} is outside [0, {size - 1}]")
        normalized.append(value)
    return normalized


def _write_fortran_field(
    stream: BinaryIO,
    field: xr.DataArray,
) -> None:
    required_dims = {"kx_folded", "kz", "y"}
    if set(field.dims) != required_dims:
        raise ValueError(f"Expected dimensions {required_dims}, received {field.dims}")

    nz = (field.sizes["kz"] - 1) // 2

    # Dataset kz order is [0, 1, ..., nz, -nz, ..., -1].
    # Reorder to solver/Fortran order [-nz, ..., -1, 0, ..., nz].
    field_fortran_kz = xr.concat(
        [
            field.isel(kz=slice(nz + 1, None)),
            field.isel(kz=slice(0, nz + 1)),
        ],
        dim="kz",
    )

    # C-order (kx, kz, y) bytes equal Fortran-order (y, kz, kx) bytes.
    values = np.asarray(
        field_fortran_kz.transpose("kx_folded", "kz", "y").values,
        dtype=OUTPUT_DTYPE,
    )
    values.tofile(stream)

def _write_statistics(
    stream: BinaryIO,
    raw_stats: xr.Dataset,
    statistics: list[str],
    label: str,
) -> None:
    """Time-average and append selected statistics in declared order."""

    for statistic in statistics:
        print(f"Computing {label}, {statistic}")
        with ProgressBar():
            field = raw_stats[statistic].mean("time").compute()
        _write_fortran_field(stream, field)


def export_raw_statistics(
    dataset: xr.Dataset,
    output_path: Path,
    scalar_indices: list[int] | None,
    time_indices: list[int] | None,
) -> None:
    """Compute and write velocity-first, scalar-major raw statistics."""

    if "isc" not in dataset["t"].dims:
        raise ValueError("Expected scalar field `t` to contain an `isc` dimension.")

    scalars = _select_indices(scalar_indices, dataset.sizes["isc"], "scalar")
    times = _select_indices(time_indices, dataset.sizes["time"], "time")
    D0_banded, D1_banded, D2_banded = get_banded_derivative_operators(
        dataset.attrs, dataset.y
    )
    field_bytes = int(
        dataset.sizes["y"]
        * dataset.sizes["kz"]
        * dataset.sizes["kx_folded"]
        * OUTPUT_DTYPE.itemsize
    )
    field_count = len(VELOCITY_RAW_STATISTIC_KEYS) + len(scalars) * len(SCALAR_RAW_STATISTIC_KEYS)
    expected_total_bytes = field_count * field_bytes

    print(
        f"Time-averaging {len(times)} snapshots and writing {field_count} fields "
        f"({expected_total_bytes / 1024**3:.3f} GiB expected)"
    )
    selected_data = dataset.isel(time=times, isc=scalars)
    raw_stats = compute_raw_statistics(
        selected_data,
        dataset.attrs,
        D0_banded,
        D1_banded,
        D2_banded,
        include_velocity=True,
    )

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("wb") as stream:
        # Velocity terms have no `isc` dimension; scalar terms retain it.
        _write_statistics(stream, raw_stats, VELOCITY_RAW_STATISTIC_KEYS, "velocity")

        for scalar_position, scalar_index in enumerate(scalars):
            _write_statistics(
                stream,
                raw_stats.isel(isc=scalar_position, drop=True),
                SCALAR_RAW_STATISTIC_KEYS,
                f"scalar {scalar_index + 1}",
            )

    print(f"Wrote raw-statistics fields to {output_path}")


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Write time-averaged raw co-spectra as headerless Fortran-readable fields."
    )
    parser.add_argument(
        "directory",
        type=Path,
        help="Directory containing convvelo.zarr or DNS output snapshots from which it can be assembled.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("raw_statistics.bin"),
        help="Binary output path, relative to DIRECTORY unless absolute.",
    )
    parser.add_argument(
        "--scalar-index",
        type=int,
        action="append",
        dest="scalar_indices",
        help="Zero-based scalar index to export; repeat to select several. Default: all.",
    )
    parser.add_argument(
        "--time-index",
        type=int,
        action="append",
        dest="time_indices",
        help="Zero-based time index contributing to the average; repeat to select several. Default: all.",
    )
    args = parser.parse_args()

    output_path = args.output if args.output.is_absolute() else args.directory / args.output

    zarr_path = args.directory / "convvelo.zarr"
    if not zarr_path.exists():
        print(f"{zarr_path} does not exist; assembling an in-memory dataset from DNS snapshots")
        dataset = build_convvelo_dataset(args.directory, pressure=True)
    else:
        dataset = xr.open_zarr(zarr_path, consolidated=False)
    dataset = dataset.chunk({"time": 1, "kz": -1})
    missing_pressure_fields = {"p", "dpdy"} - set(dataset.data_vars)
    if missing_pressure_fields:
        raise ValueError(
            f"{zarr_path} is missing {sorted(missing_pressure_fields)}; "
            "rebuild it with pressure fields enabled."
        )
    export_raw_statistics(
        dataset,
        output_path,
        args.scalar_indices,
        args.time_indices,
    )


if __name__ == "__main__":
    main()

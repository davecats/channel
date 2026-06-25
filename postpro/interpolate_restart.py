#!/usr/bin/env python3
"""Interpolate a restart field onto the mesh described by ./dns.in."""

from __future__ import annotations

import argparse
import configparser
from pathlib import Path

import dask
import numpy as np
from dask.diagnostics import ProgressBar
from scipy.interpolate import interp1d


NOISE_AMPLITUDE = 5.0e-6
OUTPUT_NAME = "Dati.cart.interp"
HEADER_BYTES = 3 * np.dtype(np.int32).itemsize + 7 * np.dtype(np.float64).itemsize
COMPLEX_BYTES = np.dtype(np.complex128).itemsize


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Interpolate a restart file onto the mesh described by ./dns.in "
            "and write ./Dati.cart.interp in the current directory."
        )
    )
    parser.add_argument("restart_file", type=Path, help="Path to the original Dati.cart*.out file")
    parser.add_argument(
        "--input-scalar",
        type=int,
        default=None,
        help="1-based scalar index from the input file to duplicate into all output scalar fields.",
    )
    parser.add_argument(
        "--chunks-x",
        type=int,
        default=16,
        help="Number of kx planes per dask task. Use <=0 for one full-width chunk.",
    )
    return parser.parse_args()


def read_restart_header(path: Path) -> tuple[dict[str, float | int], int]:
    with path.open("rb") as handle:
        header = {
            "nx": np.fromfile(handle, dtype=np.int32, count=1)[0],
            "ny": np.fromfile(handle, dtype=np.int32, count=1)[0],
            "nz": np.fromfile(handle, dtype=np.int32, count=1)[0],
            "alfa0": np.fromfile(handle, dtype=np.float64, count=1)[0],
            "beta0": np.fromfile(handle, dtype=np.float64, count=1)[0],
            "ni": np.fromfile(handle, dtype=np.float64, count=1)[0],
            "a": np.fromfile(handle, dtype=np.float64, count=1)[0],
            "ymin": np.fromfile(handle, dtype=np.float64, count=1)[0],
            "ymax": np.fromfile(handle, dtype=np.float64, count=1)[0],
            "time": np.fromfile(handle, dtype=np.float64, count=1)[0],
        }

    base_points = (int(header["nx"]) + 1) * (2 * int(header["nz"]) + 1) * (int(header["ny"]) + 3)
    payload_bytes = path.stat().st_size - HEADER_BYTES
    if base_points == 0 or payload_bytes <= 0 or payload_bytes % (base_points * COMPLEX_BYTES) != 0:
        raise ValueError(f"{path} does not contain a valid restart payload")
    return header, payload_bytes // (base_points * COMPLEX_BYTES)


def read_target_config() -> dict[str, float | int]:
    cfg = configparser.ConfigParser()
    with Path("dns.in").open("r", encoding="utf-8") as handle:
        cfg.read_file(handle)

    return {
        "nx": cfg.getint("mesh", "nx"),
        "ny": cfg.getint("mesh", "ny"),
        "nz": cfg.getint("mesh", "nz"),
        "alfa0": cfg.getfloat("mesh", "alfa0"),
        "beta0": cfg.getfloat("mesh", "beta0"),
        "ni": cfg.getfloat("velocity", "ni"),
        "a": cfg.getfloat("mesh", "stretching"),
        "ymin": cfg.getfloat("mesh", "ymin"),
        "ymax": cfg.getfloat("mesh", "ymax"),
        "nphi": cfg.getint("scalars", "nphi"),
    }


def build_y_grid(ny: int, ymin: float, ymax: float, a: float) -> np.ndarray:
    iy = np.arange(-1, ny + 2, dtype=np.float64)
    return ymin + 0.5 * (ymax - ymin) * (np.tanh(a * (2.0 * iy / ny - 1.0)) / np.tanh(a) + 1.0)


def make_noise(shape: tuple[int, ...]) -> np.ndarray:
    phase = np.random.random(shape) * 2.0 * np.pi
    return NOISE_AMPLITUDE * np.exp(1.0j * phase)


def resolve_input_scalar(nphi_in: int, nphi_out: int, input_scalar: int | None) -> int | None:
    if nphi_in > 0 and nphi_out > 0 and input_scalar is None:
        input_scalar = int(input(f"Input file contains {nphi_in} scalar(s). Select one to use [1-{nphi_in}]: "))

    if nphi_in == 0:
        return None
    if input_scalar is not None and not (1 <= input_scalar <= nphi_in):
        raise ValueError(f"--input-scalar must be between 1 and {nphi_in}")
    return input_scalar


def interpolate_restart_dask(
    source_path: Path,
    source_header: dict[str, float | int],
    n_components_in: int,
    target: dict[str, float | int],
    output_path: Path,
    input_scalar: int | None,
    chunks_x: int,
) -> None:
    nphi_in = n_components_in - 3
    nphi_out = int(target["nphi"])
    input_scalar = resolve_input_scalar(nphi_in, nphi_out, input_scalar)

    y_old = build_y_grid(
        int(source_header["ny"]),
        float(source_header["ymin"]),
        float(source_header["ymax"]),
        float(source_header["a"]),
    )
    y_new = build_y_grid(
        int(target["ny"]),
        float(target["ymin"]),
        float(target["ymax"]),
        float(target["a"]),
    )
    nxf_out = int(target["nx"]) + 1
    if chunks_x <= 0:
        chunks_x = nxf_out

    @dask.delayed
    def _write_header() -> None:
        with output_path.open("wb") as handle:
            np.array(target["nx"], dtype=np.int32).tofile(handle)
            np.array(target["ny"], dtype=np.int32).tofile(handle)
            np.array(target["nz"], dtype=np.int32).tofile(handle)
            np.array(target["alfa0"], dtype=np.float64).tofile(handle)
            np.array(target["beta0"], dtype=np.float64).tofile(handle)
            np.array(1 / target["ni"], dtype=np.float64).tofile(handle)
            np.array(target["a"], dtype=np.float64).tofile(handle)
            np.array(target["ymin"], dtype=np.float64).tofile(handle)
            np.array(target["ymax"], dtype=np.float64).tofile(handle)
            np.array(0, dtype=np.float64).tofile(handle)

    @dask.delayed
    def _interpolate_and_write_component(component_out: int, x_start: int, x_stop: int, _header_done: object) -> None:
        src_nxf = int(source_header["nx"]) + 1
        src_nzf = 2 * int(source_header["nz"]) + 1
        src_nyf = int(source_header["ny"]) + 3
        dst_nxf = int(target["nx"]) + 1
        dst_nzf = 2 * int(target["nz"]) + 1
        dst_nyf = int(target["ny"]) + 3

        slab = make_noise((x_stop - x_start, dst_nzf, dst_nyf))
        nx_common = min(int(source_header["nx"]), int(target["nx"]))

        if x_start <= nx_common and (component_out < 3 or nphi_in > 0):
            x_common_stop = min(x_stop, nx_common + 1)
            src_nz = int(source_header["nz"])
            dst_nz = int(target["nz"])
            nz_common = min(src_nz, dst_nz)
            component_in = component_out if component_out < 3 else 2 + int(input_scalar)
            source_plane_bytes = src_nzf * src_nyf * COMPLEX_BYTES
            source_offset = HEADER_BYTES + (component_in * src_nxf + x_start) * source_plane_bytes
            source_count = (x_common_stop - x_start) * src_nzf * src_nyf

            with source_path.open("rb") as handle:
                handle.seek(source_offset)
                source_slab = np.fromfile(handle, dtype=np.complex128, count=source_count)

            source_slab = source_slab.reshape((x_common_stop - x_start, src_nzf, src_nyf))
            interp_values = interp1d(y_old, source_slab[:, src_nz - nz_common:src_nz + nz_common + 1, :], axis=-1)(y_new)
            slab[: x_common_stop - x_start, dst_nz - nz_common:dst_nz + nz_common + 1, :] = interp_values

        dest_plane_bytes = dst_nzf * dst_nyf * COMPLEX_BYTES
        dest_offset = HEADER_BYTES + (component_out * dst_nxf + x_start) * dest_plane_bytes
        with output_path.open("r+b") as handle:
            handle.seek(dest_offset)
            np.ascontiguousarray(slab).tofile(handle)

    header_task = _write_header()
    write_tasks = []
    for x_start in range(0, nxf_out, chunks_x):
        x_stop = min(x_start + chunks_x, nxf_out)
        for component_out in range(3 + nphi_out):
            write_tasks.append(_interpolate_and_write_component(component_out, x_start, x_stop, header_task))

    with ProgressBar():
        dask.compute(*write_tasks)


def main() -> None:
    args = parse_args()
    source_header, n_components_in = read_restart_header(args.restart_file)
    target = read_target_config()
    interpolate_restart_dask(
        args.restart_file,
        source_header,
        n_components_in,
        target,
        Path.cwd() / OUTPUT_NAME,
        args.input_scalar,
        args.chunks_x,
    )


if __name__ == "__main__":
    main()

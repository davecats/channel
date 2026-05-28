#!/usr/bin/env python3
"""Interpolate a restart field onto the mesh described by ./dns.in."""

from __future__ import annotations

import argparse
import configparser
from pathlib import Path

import numpy as np
from scipy.interpolate import interp1d


NOISE_AMPLITUDE = 5.0e-6
OUTPUT_NAME = "Dati.cart.interp"


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
    return parser.parse_args()


def read_restart(path: Path) -> tuple[dict[str, float | int], np.ndarray]:
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

        remaining = np.fromfile(handle, dtype=np.complex128)

    base_points = (header["nx"] + 1) * (2 * header["nz"] + 1) * (header["ny"] + 3)
    if base_points == 0 or remaining.size % base_points != 0:
        raise ValueError(f"{path} does not contain a valid restart payload")

    n_components = remaining.size // base_points
    field = remaining.reshape(
        (n_components, header["nx"] + 1, 2 * header["nz"] + 1, header["ny"] + 3)
    )
    return header, field


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


def interpolate_field(
    source_header: dict[str, float | int],
    source_field: np.ndarray,
    target: dict[str, float | int],
    input_scalar: int | None,
) -> np.ndarray:
    nphi_out = int(target["nphi"])
    out = make_noise((3 + nphi_out, int(target["nx"]) + 1, 2 * int(target["nz"]) + 1, int(target["ny"]) + 3))

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

    nx_common = min(int(source_header["nx"]), int(target["nx"]))
    nz_common = min(int(source_header["nz"]), int(target["nz"]))
    nphi_in = source_field.shape[0] - 3
    src_nz = int(source_header["nz"])
    dst_nz = int(target["nz"])

    if nphi_in > 0 and nphi_out > 0 and input_scalar is None:
        input_scalar = int(input(f"Input file contains {nphi_in} scalar(s). Select one to use [1-{nphi_in}]: "))

    if nphi_in == 0:
        input_scalar = None
    elif input_scalar is not None and not (1 <= input_scalar <= nphi_in):
        raise ValueError(f"--input-scalar must be between 1 and {nphi_in}")

    src_slice = np.s_[:, : nx_common + 1, src_nz - nz_common:src_nz + nz_common + 1, :]
    dst_slice = np.s_[:, : nx_common + 1, dst_nz - nz_common:dst_nz + nz_common + 1, :]
    interp = interp1d(y_old, source_field[src_slice], axis=-1)
    interp_values = interp(y_new)

    out[dst_slice][0:3] = interp_values[0:3]

    if nphi_in > 0 and nphi_out > 0:
        comp_in = 2 + input_scalar
        scalar_block = interp_values[comp_in]
        out[3:3 + nphi_out, : nx_common + 1, dst_nz - nz_common:dst_nz + nz_common + 1, :] = scalar_block[None, ...]

    return out


def write_restart(path: Path, target: dict[str, float | int], time: float, field: np.ndarray) -> None:
    with path.open("wb") as handle:
        np.array(target["nx"], dtype=np.int32).tofile(handle)
        np.array(target["ny"], dtype=np.int32).tofile(handle)
        np.array(target["nz"], dtype=np.int32).tofile(handle)
        np.array(target["alfa0"], dtype=np.float64).tofile(handle)
        np.array(target["beta0"], dtype=np.float64).tofile(handle)
        np.array(target["ni"], dtype=np.float64).tofile(handle)
        np.array(target["a"], dtype=np.float64).tofile(handle)
        np.array(target["ymin"], dtype=np.float64).tofile(handle)
        np.array(target["ymax"], dtype=np.float64).tofile(handle)
        np.array(time, dtype=np.float64).tofile(handle)
        field.tofile(handle)


def main() -> None:
    args = parse_args()
    source_header, source_field = read_restart(args.restart_file)
    target = read_target_config()
    out_field = interpolate_field(source_header, source_field, target, args.input_scalar)
    write_restart(Path.cwd() / OUTPUT_NAME, target, 0, out_field)


if __name__ == "__main__":
    main()

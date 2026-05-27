import argparse
import configparser
from pathlib import Path
import warnings

import dask
import dask.array as da
from dask.diagnostics import ProgressBar
from dask.distributed import Client, LocalCluster
import numpy as np
import xarray as xr

from utils import fftfit, mesh2k

CONVVELO_HEADER_BYTES = 2 * np.dtype(np.float64).itemsize + np.dtype(np.int64).itemsize


def _get_ini_config(path: Path) -> configparser.ConfigParser:
    cfg = configparser.ConfigParser(inline_comment_prefixes=("!", "#"))
    cfg.read(path)
    return cfg


def get_sim_metadata(directory):
    path = Path(directory) / "dns.in"
    cfg = _get_ini_config(path)
    nx = cfg.getint("mesh", "nx")
    ny = cfg.getint("mesh", "ny")
    nz = cfg.getint("mesh", "nz")
    alfa0 = cfg.getfloat("mesh", "alfa0")
    beta0 = cfg.getfloat("mesh", "beta0")
    ni = 1.0 / cfg.getfloat("velocity", "ni")
    a = cfg.getfloat("mesh", "stretching")
    ymin = cfg.getfloat("mesh", "ymin")
    ymax = cfg.getfloat("mesh", "ymax")
    nphi = cfg.getint("scalars", "nphi", fallback=0)
    pr_array = np.fromstring(cfg.get("scalars", "pr", fallback=""), sep=" ")
    if nphi == 0:
        pr_array = np.array([], dtype=np.float64)

    def next_fft_size(n):
        while not fftfit(n):
            n += 1
        return n

    nxd = next_fft_size(3 * (nx + 1) // 2)
    nzd = next_fft_size(3 * nz)

    y = ymin + 0.5 * (ymax - ymin) * (
        np.tanh(a * (2 * np.arange(-1, ny + 2) / ny - 1)) / np.tanh(a) + 1
    )
    z = np.arange(nzd + 1) / (nzd + 1) * (2 * np.pi / beta0)
    x = np.arange(2 * nxd + 1) / (2 * nxd + 1) * (2 * np.pi / alfa0)
    kx_folded = mesh2k(alfa0, nx)[:nx + 1]
    kz = mesh2k(beta0, nz)

    return dict(
        nx=nx, ny=ny, nz=nz,
        alfa0=alfa0, beta0=beta0,
        ni=ni, a=a, ymin=ymin, ymax=ymax,
        nxd=nxd, nzd=nzd, pr=pr_array,
        x=x, y=y, z=z,
        kx_folded=kx_folded,
        kz=kz,
    )


def _get_convvelo_runtime_config(directory):
    cfg = _get_ini_config(Path(directory) / "dns.in")
    if not cfg.has_section("convvelo"):
        return None
    return {
        "mode": cfg.get("convvelo", "output_mode", fallback="full").lower(),
        "output_file": "convvelo.bin",
        "dt_write": cfg.getfloat("convvelo", "dt_write", fallback=-1.0),
    }


def _read_convvelo_field_sidecar(layout_base_path):
    sidecar = layout_base_path.with_suffix(layout_base_path.suffix + ".fields")
    if not sidecar.exists():
        raise FileNotFoundError(f"Missing convvelo field layout file for {layout_base_path}")

    parsed = {}
    with open(sidecar) as stream:
        for raw_line in stream:
            line = raw_line.split("#", 1)[0].strip()
            if not line:
                continue
            if ":" not in line:
                raise ValueError(f"Malformed convvelo field sidecar line in {sidecar}: {raw_line.rstrip()}")
            key, values = line.split(":", 1)
            parsed[key.strip().lower()] = values.split()

    required = {"profile_fields", "velocity_fields", "scalar_fields"}
    missing = required - set(parsed)
    if missing:
        raise ValueError(f"Convvelo field sidecar {sidecar} is missing sections: {sorted(missing)}")

    return {
        "sidecar": str(sidecar),
        "profile_fields": parsed["profile_fields"],
        "velocity_fields": parsed["velocity_fields"],
        "scalar_fields": parsed["scalar_fields"],
    }


def _convvelo_profile_header_slots(nphi):
    return 3 + nphi


def _convvelo_expected_size(metadata, velocity_fields, scalar_fields):
    nyf = metadata["ny"] + 3
    nxf = metadata["nx"] + 1
    nzf = 2 * metadata["nz"] + 1
    nphi = len(metadata["pr"])
    profile_bytes = _convvelo_profile_header_slots(nphi) * nyf * np.dtype(np.complex128).itemsize
    field_bytes = (len(velocity_fields) + nphi * len(scalar_fields)) * nxf * nzf * nyf * np.dtype(np.complex128).itemsize
    return CONVVELO_HEADER_BYTES + profile_bytes + field_bytes


def _read_convvelo_header(path):
    with open(path, "rb") as stream:
        average_start_time, average_end_time = np.fromfile(stream, dtype=np.float64, count=2)
        average_count = np.fromfile(stream, dtype=np.int64, count=1)[0]
    return {
        "average_start_time": float(average_start_time),
        "average_end_time": float(average_end_time),
        "average_count": int(average_count),
    }


def _find_convvelo_runtime_files(directory):
    config = _get_convvelo_runtime_config(directory)
    if config is None:
        return [], None

    base_path = Path(directory) / config["output_file"]
    numbered = sorted(base_path.parent.glob(f"{base_path.stem}.*{base_path.suffix}"))
    return numbered, config


def _read_convvelo_info(runtime_path, metadata, velocity_fields, scalar_fields):
    _sanity_check_convvelo_runtime(runtime_path, metadata, velocity_fields, scalar_fields)
    header = _read_convvelo_header(runtime_path)
    return {
        "path": runtime_path,
        "start_time": header["average_start_time"],
        "end_time": header["average_end_time"],
        "nfields": header["average_count"],
    }


def _average_convvelo_steps(data):
    if "step" not in data.dims:
        return data

    print("Convvelo averaging windows:")
    for step, start_time, end_time, count in zip(
        np.asarray(data["step"].values),
        np.asarray(data["start_time"].values),
        np.asarray(data["end_time"].values),
        np.asarray(data["nfields"].values),
        strict=True,
    ):
        print(f"  step={step}: start_time={start_time:.14f}, end_time={end_time:.14f}, nfields={int(count)}")

    weights = data["nfields"].astype(np.float64)
    averaged = data.drop_vars(["start_time", "end_time", "nfields"]).weighted(weights).mean("step")
    averaged.attrs = data.attrs
    averaged["start_time"] = data["start_time"].isel(step=0, drop=True)
    averaged["end_time"] = data["end_time"].isel(step=-1, drop=True)
    averaged["nfields"] = data["nfields"].sum("step")
    return averaged


def _sanity_check_convvelo_runtime(runtime_path, metadata, velocity_fields, scalar_fields):
    actual_size = runtime_path.stat().st_size
    expected_size = _convvelo_expected_size(metadata, velocity_fields, scalar_fields)
    if actual_size != expected_size:
        raise ValueError(
            f"{runtime_path} has size {actual_size}, expected {expected_size} for convvelo runtime output"
        )

    header = _read_convvelo_header(runtime_path)
    if header["average_count"] <= 0:
        raise ValueError(f"{runtime_path} stores a non-positive average_count")

    nyf = metadata["ny"] + 3
    nzf = 2 * metadata["nz"] + 1
    nxf = metadata["nx"] + 1
    profile_bytes = nyf * np.dtype(np.complex128).itemsize
    field_bytes = nxf * nzf * nyf * np.dtype(np.complex128).itemsize
    nphi = len(metadata["pr"])

    with open(runtime_path, "rb") as stream:
        for profile_slot, name in enumerate(("mean_u", "mean_v", "mean_w")):
            stream.seek(CONVVELO_HEADER_BYTES + profile_slot * profile_bytes)
            profile = np.fromfile(stream, dtype=np.complex128, count=nyf)
            if profile.size != nyf or not np.isfinite(profile).all():
                raise ValueError(f"{runtime_path} failed sanity check for {name}")

        for i_phi in range(nphi):
            stream.seek(CONVVELO_HEADER_BYTES + (3 + i_phi) * profile_bytes)
            profile = np.fromfile(stream, dtype=np.complex128, count=nyf)
            if profile.size != nyf or not np.isfinite(profile).all():
                raise ValueError(f"{runtime_path} failed sanity check for mean_t")

        stream.seek(CONVVELO_HEADER_BYTES + _convvelo_profile_header_slots(nphi) * profile_bytes)
        first_velocity = np.fromfile(stream, dtype=np.complex128, count=min(nxf * nzf * nyf, 64))
        if first_velocity.size == 0 or not np.isfinite(first_velocity).all():
            raise ValueError(f"{runtime_path} failed sanity check for first velocity field")

        if nphi > 0:
            first_scalar_offset = (
                CONVVELO_HEADER_BYTES
                + _convvelo_profile_header_slots(nphi) * profile_bytes
                + len(velocity_fields) * field_bytes
            )
            stream.seek(first_scalar_offset)
            first_scalar = np.fromfile(stream, dtype=np.complex128, count=min(nxf * nzf * nyf, 64))
            if first_scalar.size == 0 or not np.isfinite(first_scalar).all():
                raise ValueError(f"{runtime_path} failed sanity check for first scalar field")


def load_convvelo_runtime(directory, metadata):
    runtime_files, config = _find_convvelo_runtime_files(directory)
    if not runtime_files:
        return None

    layout_base_path = Path(directory) / config["output_file"]
    field_info = _read_convvelo_field_sidecar(layout_base_path)
    velocity_fields = field_info["velocity_fields"]
    scalar_fields = field_info["scalar_fields"]
    nphi = len(metadata["pr"])
    nyf = metadata["ny"] + 3
    nzf = 2 * metadata["nz"] + 1
    nxf = metadata["nx"] + 1
    profile_bytes = nyf * np.dtype(np.complex128).itemsize
    field_bytes = nxf * nzf * nyf * np.dtype(np.complex128).itemsize
    field_offset = CONVVELO_HEADER_BYTES + _convvelo_profile_header_slots(nphi) * profile_bytes

    file_info = [_read_convvelo_info(path, metadata, velocity_fields, scalar_fields) for path in runtime_files]
    paths = [item["path"] for item in file_info]
    step_values = np.arange(len(paths), dtype=np.int64)

    def _load_complex_slice(path, offset, shape):
        with open(path, "rb") as stream:
            stream.seek(offset)
            buf = stream.read(int(np.prod(shape)) * np.dtype(np.complex128).itemsize)
        return np.frombuffer(buf, dtype=np.complex128).reshape(shape)

    def _stack_delayed(offsets, shape):
        return da.stack([
            da.from_delayed(
                dask.delayed(_load_complex_slice)(path, offset, shape),
                shape=shape,
                dtype=np.complex128,
            )
            for path, offset in zip(paths, offsets, strict=True)
        ])

    coords = {
        "step": step_values,
        "isc": np.arange(nphi),
        "Pr": ("isc", metadata["pr"]),
        "kx_folded": metadata["kx_folded"],
        "kz": metadata["kz"],
        "y": metadata["y"],
    }

    profile_offsets = [CONVVELO_HEADER_BYTES + i * profile_bytes for i in range(3 + nphi)]
    profile_data = {
        "mean_u": (("step", "y"), _stack_delayed([profile_offsets[0]] * len(paths), (nyf,))),
        "mean_v": (("step", "y"), _stack_delayed([profile_offsets[1]] * len(paths), (nyf,))),
        "mean_w": (("step", "y"), _stack_delayed([profile_offsets[2]] * len(paths), (nyf,))),
    }

    if nphi > 0:
        profile_data["mean_t"] = (
            ("step", "isc", "y"),
            da.stack(
                [
                    _stack_delayed([profile_offsets[3 + i_phi]] * len(paths), (nyf,))
                    for i_phi in range(nphi)
                ],
                axis=1,
            ),
        )

    spectral_fields = {
        name: (
            ("step", "kx_folded", "kz", "y"),
            _stack_delayed([field_offset + i_field * field_bytes] * len(paths), (nxf, nzf, nyf)),
        )
        for i_field, name in enumerate(velocity_fields)
    }
    spectral_fields.update({
        name: (
            ("step", "isc", "kx_folded", "kz", "y"),
            da.stack(
                [
                    _stack_delayed(
                        [field_offset + (len(velocity_fields) + i_phi * len(scalar_fields) + i_field) * field_bytes] * len(paths),
                        (nxf, nzf, nyf),
                    )
                    for i_phi in range(nphi)
                ],
                axis=1,
            ),
        )
        for i_field, name in enumerate(scalar_fields)
    })

    ds = xr.Dataset(
        {
            **profile_data,
            **spectral_fields,
            "start_time": ("step", np.array([item["start_time"] for item in file_info], dtype=np.float64)),
            "end_time": ("step", np.array([item["end_time"] for item in file_info], dtype=np.float64)),
            "nfields": ("step", np.array([item["nfields"] for item in file_info], dtype=np.int64)),
        },
        coords=coords,
    ).roll(kz=-metadata["nz"])
    ds = ds.sortby("start_time")

    ds.attrs = {
        "convvelo_runtime_mode": config["mode"],
        "convvelo_runtime_sidecar": field_info["sidecar"] or "",
        "convvelo_runtime_files": tuple(str(path) for path in paths),
        "convvelo_velocity_fields": tuple(velocity_fields),
        "convvelo_scalar_fields": tuple(scalar_fields),
        **{k: v for k, v in metadata.items() if not isinstance(v, np.ndarray)},
    }
    return _average_convvelo_steps(ds)


def build_convvelo_dataset(directory, *, pressure=False, command_semaphore=None):
    if pressure:
        warnings.warn("pressure=True is ignored: convvelo datasets now come only from online runtime files.")
    if command_semaphore is not None:
        warnings.warn("command_semaphore is ignored: offline field reconstruction has been removed.")

    directory = Path(directory)
    sim_meta = get_sim_metadata(directory)
    data = load_convvelo_runtime(directory, sim_meta)
    if data is None:
        raise FileNotFoundError(f"No convvelo runtime files found in {directory}")
    data["Pr"] = ("isc", sim_meta["pr"])
    return data


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Convert convvelo runtime files to a step-indexed xarray Dataset and save as Zarr.")
    parser.add_argument("directory", type=Path, help="Directory containing convvelo runtime output and dns.in")
    parser.add_argument("--chunks-x", type=int, default=-1, help="Chunk size for kx_folded on write")
    parser.add_argument("--chunks-y", type=int, default=-1, help="Chunk size for y on write")
    args = parser.parse_args()

    cluster = LocalCluster(n_workers=8, threads_per_worker=2)
    client = Client(cluster)

    data = build_convvelo_dataset(args.directory)
    print(data)

    with ProgressBar():
        chunk_map = {"y": args.chunks_y, "kx_folded": args.chunks_x}
        data = data.chunk({name: size for name, size in chunk_map.items() if name in data.dims})
        data.to_zarr(args.directory / "convvelo.zarr", mode="w", consolidated=False)

    client.close()
    cluster.close()

#!/usr/bin/env python3
import argparse
import math
import struct
import sys
from pathlib import Path


def fail(message: str) -> int:
    print(f"ERROR: {message}", file=sys.stderr)
    return 1


def profile_value(physical_y: int, component: int) -> complex:
    return complex(100.0 * component + physical_y,
                   -10.0 * component + 0.125 * physical_y)


def field_value(physical_y: int, iz: int, ix: int, field_index: int) -> complex:
    return complex(1000.0 * field_index + 10.0 * ix + 0.5 * iz + 0.03125 * physical_y,
                   -200.0 * field_index + 0.25 * ix - 0.0625 * iz + 0.0078125 * physical_y)


def check_complex(got: complex, expected: complex, atol: float, label: str) -> int:
    diff = abs(got - expected)
    if math.isnan(got.real) or math.isnan(got.imag) or diff > atol:
        return fail(f"{label}: got ({got.real:.17e},{got.imag:.17e}), "
                    f"expected ({expected.real:.17e},{expected.imag:.17e}), diff={diff:.6e}")
    return 0


def main() -> int:
    parser = argparse.ArgumentParser(description="Verify a synthetic convvelo raw binary file.")
    parser.add_argument("filename", type=Path)
    parser.add_argument("--nx", type=int, required=True)
    parser.add_argument("--ny", type=int, required=True)
    parser.add_argument("--nz", type=int, required=True)
    parser.add_argument("--nphi", type=int, required=True)
    parser.add_argument("--atol", type=float, default=0.0)
    args = parser.parse_args()

    data = args.filename.read_bytes()
    ny_count = args.ny + 3
    z_count = 2 * args.nz + 1
    x_count = args.nx + 1
    n_profiles = 3 + args.nphi
    n_fields = 33 + args.nphi * 10
    header_bytes = 24
    profile_bytes = 16 * ny_count
    field_bytes = 16 * ny_count * z_count * x_count
    expected_size = header_bytes + n_profiles * profile_bytes + n_fields * field_bytes

    if len(data) != expected_size:
        return fail(f"size mismatch: got {len(data)} bytes, expected {expected_size}")

    start_time, end_time = struct.unpack_from("<2d", data, 0)
    sample_count = struct.unpack_from("<q", data, 16)[0]
    if abs(start_time - 1.25) > args.atol or abs(end_time - 2.50) > args.atol:
        return fail(f"header time mismatch: got start={start_time}, end={end_time}")
    if sample_count != 7:
        return fail(f"sample count mismatch: got {sample_count}, expected 7")

    offset = header_bytes
    for profile_index in range(n_profiles):
        component = profile_index + 1
        for y_file in range(ny_count):
            physical_y = y_file - 1
            real, imag = struct.unpack_from("<dd", data, offset)
            status = check_complex(complex(real, imag), profile_value(physical_y, component),
                                   args.atol, f"profile={profile_index} y={physical_y}")
            if status:
                return status
            offset += 16

    for field_zero_based in range(n_fields):
        field_index = field_zero_based + 1
        for ix in range(x_count):
            for z_file in range(z_count):
                iz = z_file - args.nz
                for y_file in range(ny_count):
                    physical_y = y_file - 1
                    real, imag = struct.unpack_from("<dd", data, offset)
                    status = check_complex(complex(real, imag), field_value(physical_y, iz, ix, field_index),
                                           args.atol,
                                           f"field={field_zero_based} y={physical_y} iz={iz} ix={ix}")
                    if status:
                        return status
                    offset += 16

    print(f"synthetic convvelo raw file verified: {args.filename}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

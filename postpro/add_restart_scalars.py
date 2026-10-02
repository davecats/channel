#!/usr/bin/env python3
"""Rewrite a restart file with a different set of passive scalars.

The restart header records the mesh but not nPhi, and the solver builds its
MPI-IO view for `3 + nPhi` components, so raising nPhi against an older file
reads past the end of it.  The reader now refuses that; this is how you get a
file it will accept.

On disk the field is one array, `(ny+3, 2*nz+1, nx+1, 3+nPhi)` complex doubles
in Fortran order, so the component index is the slowest -- each component is a
contiguous block and this is a copy, not a transform.  Nothing is interpolated;
the mesh is unchanged.  Use `interpolate_restart.py` when the mesh changes.

    # three scalars -> six, the new three starting from copies of the old three
    ./add_restart_scalars.py Dati.cart.out Dati.cart.6phi.out --scalars 1 2 3 1 2 3

    # one scalar out of three, for a cheap test case
    ./add_restart_scalars.py Dati.cart.out Dati.cart.1phi.out --scalars 2

    # six scalars, the new three starting from zero instead
    ./add_restart_scalars.py Dati.cart.out Dati.cart.6phi.out --scalars 1 2 3 0 0 0

`--scalars` lists the 1-based scalar indices of the *input* in the order the
output should hold them; `0` means a zero-filled scalar.  Starting the isoflux
scalars from copies of the matching isothermal ones is the point of the tool:
it puts them near equilibrium, instead of relaxing a whole mean profile, which
at Pr = 0.025 is the difference between a short transient and a very long one.
"""

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np

HEADER_BYTES = 3 * np.dtype(np.int32).itemsize + 7 * np.dtype(np.float64).itemsize
COMPLEX_BYTES = np.dtype(np.complex128).itemsize
CHUNK_ELEMENTS = 1 << 24  # 256 MiB of complex128 per read


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("input", type=Path, help="restart file to read (Dati.cart*.out)")
    parser.add_argument("output", type=Path, help="restart file to write; must not exist")
    parser.add_argument("--scalars", type=int, nargs="+", required=True, metavar="N",
                        help="1-based input scalar indices, in output order; 0 for a zero-filled scalar")
    parser.add_argument("--force", action="store_true", help="overwrite the output if it exists")
    return parser.parse_args()


def read_header(path: Path) -> tuple[dict[str, int | float], int]:
    """The header fields, and how many components the file actually holds."""
    with path.open("rb") as handle:
        ints = np.fromfile(handle, dtype=np.int32, count=3)
        reals = np.fromfile(handle, dtype=np.float64, count=7)
    if ints.size != 3 or reals.size != 7:
        raise SystemExit(f"{path}: too short to hold a restart header")
    header = dict(zip(("nx", "ny", "nz"), (int(v) for v in ints)))
    header.update(zip(("alfa0", "beta0", "ni", "a", "ymin", "ymax", "time"),
                      (float(v) for v in reals)))

    elements = component_elements(header)
    payload = path.stat().st_size - HEADER_BYTES
    if payload <= 0 or payload % (elements * COMPLEX_BYTES) != 0:
        raise SystemExit(f"{path}: {payload} payload bytes is not a whole number of "
                         f"{elements * COMPLEX_BYTES}-byte components")
    return header, payload // (elements * COMPLEX_BYTES)


def component_elements(header: dict[str, int | float]) -> int:
    return (int(header["ny"]) + 3) * (2 * int(header["nz"]) + 1) * (int(header["nx"]) + 1)


def copy_component(src, dst, source_index: int | None, elements: int) -> None:
    """Append one component to dst: a copy of source_index, or zeros."""
    if source_index is not None:
        src.seek(HEADER_BYTES + source_index * elements * COMPLEX_BYTES)
    remaining = elements
    while remaining > 0:
        count = min(remaining, CHUNK_ELEMENTS)
        if source_index is None:
            block = np.zeros(count, dtype=np.complex128)
        else:
            block = np.fromfile(src, dtype=np.complex128, count=count)
            if block.size != count:
                raise SystemExit("input ended early; is another process writing it?")
        block.tofile(dst)
        remaining -= count


def main() -> None:
    args = parse_args()
    if args.output.exists() and not args.force:
        raise SystemExit(f"{args.output} exists; pass --force to overwrite")

    header, n_components = read_header(args.input)
    nphi_in = n_components - 3
    if nphi_in < 0:
        raise SystemExit(f"{args.input} holds {n_components} components, fewer than the three velocities")
    for index in args.scalars:
        if not 0 <= index <= nphi_in:
            raise SystemExit(f"--scalars {index}: the input has {nphi_in} scalar(s), so 0..{nphi_in}")

    elements = component_elements(header)
    print(f"{args.input}: nx = {header['nx']}, ny = {header['ny']}, nz = {header['nz']}, "
          f"time = {header['time']}, nPhi = {nphi_in}")
    print(f"{args.output}: nPhi = {len(args.scalars)}, scalars from "
          + ", ".join("zero" if i == 0 else f"input {i}" for i in args.scalars))

    with args.input.open("rb") as src, args.output.open("wb") as dst:
        # The header is copied verbatim: it says nothing about nPhi, which is
        # exactly why the size has to be right.
        src.seek(0)
        dst.write(src.read(HEADER_BYTES))
        for component in range(3):
            copy_component(src, dst, component, elements)
        for index in args.scalars:
            copy_component(src, dst, None if index == 0 else 2 + index, elements)

    written = args.output.stat().st_size
    expected = HEADER_BYTES + (3 + len(args.scalars)) * elements * COMPLEX_BYTES
    if written != expected:
        raise SystemExit(f"wrote {written} bytes, expected {expected}")
    print(f"wrote {written} bytes ({3 + len(args.scalars)} components)")


if __name__ == "__main__":
    main()

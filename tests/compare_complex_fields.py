#!/usr/bin/env python3
import sys
import numpy as np


def main() -> int:
    if len(sys.argv) != 4:
        print("Usage: compare_complex_fields.py <actual> <reference> <tol>")
        return 2

    actual_path, reference_path, tol_text = sys.argv[1:]
    tol = float(tol_text)

    actual = np.fromfile(actual_path, dtype=np.complex128)
    reference = np.fromfile(reference_path, dtype=np.complex128)

    if actual.shape != reference.shape:
        print(
            f"Shape mismatch: {actual_path} has {actual.size} values, "
            f"{reference_path} has {reference.size}"
        )
        return 1

    diff = np.abs(actual - reference)
    max_err = np.max(diff) if diff.size else 0.0
    print(f"{actual_path}: max abs error = {max_err}")

    if not np.isfinite(max_err) or max_err >= tol:
        return 1

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

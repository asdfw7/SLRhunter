#!/usr/bin/env python3
"""Wrapper for the bundled matrixer executable."""

from __future__ import annotations

import subprocess
import sys
from pathlib import Path


def resolve_matrixer() -> str:
    """Find the bundled matrixer binary."""
    pkg_bin = Path(__file__).resolve().parent / "bin" / "matrixer"
    if pkg_bin.exists():
        return str(pkg_bin)

    raise FileNotFoundError(
        "matrixer executable was not found. Please install the package or add matrixer to PATH."
    )


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else list(argv)
    exe = resolve_matrixer()
    proc = subprocess.run([exe, *argv], check=False)
    return int(proc.returncode)


if __name__ == "__main__":
    raise SystemExit(main())

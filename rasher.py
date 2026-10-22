"""Backward-compatible entry point: `python rasher.py ...` is the same as `rasher ...`."""

import sys

from rasher.cli import main

if __name__ == "__main__":
    sys.exit(main())

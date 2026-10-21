#!/usr/bin/env python3
"""従来の ``python3 generate.py {ohp,gaia-etch}`` 用ラッパー。``gingaforge`` コマンドと同じです。"""

from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent / "src"))

from gingaforge.cli import main

if __name__ == "__main__":
    raise SystemExit(main())

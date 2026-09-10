#!/usr/bin/env python3
# Copyright (C) 2026 東京大学地文研究会天文部; GPL-2.0-or-later
"""GingaForge: OHPフィルム・Gaiaエッチング用の天の川原盤。"""

import argparse
import sys


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog="方式別の設定・オプション: generate.py {ohp,gaia-etch} --help",
    )
    parser.add_argument("mode", choices=("ohp", "gaia-etch"),
                        help="ohp: Mellinger全天画像 / gaia-etch: Gaia DR3光量")
    parser.add_argument("options", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if args.mode == "gaia-etch":
        from galaxy_transformer import main as generate
    else:
        from build_mellinger_ohp import main as generate
    return generate(args.options)


if __name__ == "__main__":
    raise SystemExit(main())

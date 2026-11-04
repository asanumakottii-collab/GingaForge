# Copyright (C) 2026 東京大学地文研究会天文部; GPL-2.0-or-later
"""GingaForge: OHPフィルム・Gaiaエッチング用の天の川原盤。"""

import argparse
import sys


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog="方式別の設定・オプション: gingaforge {gaia-etch,gaia-ohp} --help",
    )
    parser.add_argument("mode", choices=("gaia-etch", "gaia-ohp"),
                        help="gaia-etch: Gaia DR3光量のSVG・PDF / "
                             "gaia-ohp: Gaia DR3光量のOHPフィルム印刷用SVG・PDF（天の川強調あり）")
    parser.add_argument("options", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    from .gaia_etch import main as generate
    if args.mode == "gaia-ohp":
        return generate(args.options, "gaia-ohp.properties", prog="gingaforge gaia-ohp")
    return generate(args.options)

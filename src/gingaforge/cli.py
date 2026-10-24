# Copyright (C) 2026 東京大学地文研究会天文部; GPL-2.0-or-later
"""GingaForge: OHPフィルム・Gaiaエッチング用の天の川原盤。"""

import argparse
import sys


def main(argv=None):
    argv = sys.argv[1:] if argv is None else argv
    parser = argparse.ArgumentParser(
        description=__doc__,
        epilog="方式別の設定・オプション: gingaforge {mellinger-ohp,gaia-etch,gaia-ohp,install-starnet} --help",
    )
    parser.add_argument("mode", choices=("mellinger-ohp", "gaia-etch", "gaia-ohp", "install-starnet"),
                        help="mellinger-ohp: Mellinger全天画像のPDF / gaia-etch: Gaia DR3光量のSVG・PDF / "
                             "gaia-ohp: Gaia DR3光量のOHPフィルム印刷用SVG・PDF（天の川強調あり） / "
                             "install-starnet: mellinger-ohpの星の除去に使うStarNet2を公式配布元から導入")
    parser.add_argument("options", nargs=argparse.REMAINDER)
    args = parser.parse_args(argv)
    if args.mode == "gaia-etch":
        from .gaia_etch import main as generate
    elif args.mode == "gaia-ohp":
        from .gaia_etch import main as gaia_etch
        return gaia_etch(args.options, "gaia-ohp.properties", prog="gingaforge gaia-ohp")
    elif args.mode == "install-starnet":
        from .starnet_setup import main as generate
    else:
        from .ohp import main as generate
    return generate(args.options)

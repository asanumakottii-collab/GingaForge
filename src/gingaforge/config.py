# -*- coding: utf-8 -*-
"""Javaのプロパティファイル形式(key = value, '#'または'!'で始まる行はコメント)を読み込むユーティリティ。"""

from pathlib import Path

# 同梱の設定ファイル・source・outputを置くリポジトリのルート（src/gingaforge から2階層上）。
PROJECT_ROOT = Path(__file__).resolve().parents[2]


def load_properties(path):
    """プロパティファイルを読み込み、文字列の辞書として返します。"""
    props = {}
    with open(path, encoding="utf-8") as f:
        for raw_line in f:
            line = raw_line.strip()
            if not line or line.startswith("#") or line.startswith("!"):
                continue
            for sep in ("=", ":"):
                if sep in line:
                    key, value = line.split(sep, 1)
                    props[key.strip()] = value.strip()
                    break
    return props

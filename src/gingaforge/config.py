# -*- coding: utf-8 -*-
"""Javaのプロパティファイル形式(key = value, '#'または'!'で始まる行はコメント)を読み込むユーティリティ。"""

import math
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


def output_directory(config_path, *parts):
    """生成物の保存先。設定ファイルと同じフォルダの ``output/<設定ファイル名>/`` に、設定ごとに分ける。"""
    config_path = Path(config_path)
    return config_path.parent.joinpath("output", config_path.stem, *parts)


def positive_number(props, key, default, *, allow_zero=False):
    try:
        value = float(props.get(key, default))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key}: 数値を指定してください。") from exc
    if not math.isfinite(value) or (value < 0 if allow_zero else value <= 0):
        raise ValueError(f"{key}: 有限の{'0以上' if allow_zero else '正'}の数値が必要です。")
    return value

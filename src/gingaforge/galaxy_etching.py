# Copyright (C) 2026 東京大学地文研究会天文部; GPL-2.0-or-later
"""Gaiaの局所光量を、最小径・金属幅を満たす円形穴へ量子化する。

穴は個々の恒星ではなく局所的な合計光量の標本。全穴を同じ直径にし、
固定ピッチの候補格子の一部を選ぶことで、原盤全体の最小金属幅を保証する。
Gaiaの集計ビンの光量は、原盤上でビンとセルが重なる面積に比例して分配する。
天の川強調を有効にすると、背景に相当する一様な明るさを差し引いてから穴を配置する。
"""

from collections import defaultdict
from array import array
from dataclasses import dataclass
import math
import random

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from .config import positive_number
from .models import SpherePosition


def _overlaps(start, end, cell):
    """区間と幅cellのセルの重なりを ``(セル番号, 区間に対する割合)`` で返す。幅0の区間は点として扱う。"""
    low, high = sorted((start, end))
    # セル境界にごく近い端は丸め誤差とみなし、隣のセルにわずかな光量を作らない。
    first = math.floor(low / cell + 1e-9)
    last = max(first, math.ceil(high / cell - 1e-9) - 1)
    if first == last:
        return ((first, 1.),)
    pieces = [(i, min(high, (i + 1) * cell) - max(low, i * cell)) for i in range(first, last + 1)]
    total = sum(length for _, length in pieces)
    return [(i, length / total) for i, length in pieces]


class _Holes:
    """大量の穴を座標配列で保持し、``(dir, index, x, y)`` として返す。"""
    def __init__(self):
        self.coordinates = array("d")

    def append(self, x, y, region):
        self.coordinates.extend((x, y, region))

    def __len__(self):
        return len(self.coordinates) // 3

    def __iter__(self):
        for i in range(0, len(self.coordinates), 3):
            x, y, region = self.coordinates[i:i + 3]
            yield (*divmod(int(region), 2), x, y)


@dataclass(frozen=True)
class EtchingSettings:
    min_diameter_mm: float = 0.02
    min_web_mm: float = 0.02
    cell_mm: float = 0.2
    flux_gain: float = 1.
    seed: int = 0
    scale: float = 1.
    # 天の川強調。既定は無効で、光量に比例した従来の穴面積になる。
    milky_way: bool = False
    milky_way_black_percentile: float = 40.
    milky_way_anchor_percentile: float = 95.
    milky_way_gamma: float = 1.
    milky_way_smoothing_mm: float = 1.

    @classmethod
    def from_properties(cls, props):
        # 金属エッチング用は galaxy.etch.*、OHP印刷用は galaxy.ohp.* と書ける。中身は同じ設定。
        prefixes = {p for p in ("galaxy.etch", "galaxy.ohp") if any(k.startswith(p + ".") for k in props)}
        if len(prefixes) > 1:
            raise ValueError("galaxy.etch.* と galaxy.ohp.* は同じ設定です。どちらか一方だけで書いてください。")
        prefix = prefixes.pop() if prefixes else "galaxy.etch"
        diameter = positive_number(props, f"{prefix}.min-hole-diameter-mm", .02)
        web = positive_number(props, f"{prefix}.min-web-mm", .02)
        cell = positive_number(props, f"{prefix}.cell-size-mm", .2)
        gain = positive_number(props, f"{prefix}.flux-gain", 1.)
        scale = positive_number(props, "scale", 1.)
        try:
            seed = int(props.get(f"{prefix}.seed", 0))
        except (ValueError, TypeError) as exc:
            raise ValueError(f"{prefix}.seed: 整数を指定してください。") from exc
        if cell < diameter + web:
            raise ValueError(f"{prefix}.cell-size-mm は最小穴径＋最小金属幅以上にしてください。")
        if cell / (diameter + web) > 100:
            raise ValueError("1セルの候補点が多すぎます。cell-size-mm / (穴径+金属幅) <= 100 としてください。")
        milky_way = props.get(f"{prefix}.milky-way", "no").strip().lower() in {"y", "yes", "true"}
        black = positive_number(props, f"{prefix}.milky-way.black-percentile", 40., allow_zero=True)
        anchor = positive_number(props, f"{prefix}.milky-way.anchor-percentile", 95.)
        gamma = positive_number(props, f"{prefix}.milky-way.gamma", 1.)
        smoothing = positive_number(props, f"{prefix}.milky-way.smoothing-mm", 1.)
        if not black < anchor <= 100:
            raise ValueError(f"{prefix}.milky-way: 0 <= black-percentile < anchor-percentile <= 100 としてください。")
        if smoothing / cell > 50:
            raise ValueError(f"{prefix}.milky-way.smoothing-mm は cell-size-mm の50倍以下にしてください。")
        return cls(diameter, web, cell, gain, seed, scale, milky_way, black, anchor, gamma, smoothing)


class EtchingBuilder:
    def __init__(self, transformer, settings, frame_radius):
        """半径 frame_radius の円形枠の内側に穴を置く。"""
        self.transformer = transformer
        self.settings = settings
        self.radius = settings.min_diameter_mm * settings.scale / 2
        # 最終寸法における丸め誤差も安全側へ寄せる。
        self.pitch = (settings.min_diameter_mm + settings.min_web_mm) * settings.scale * (1 + 1e-10)
        self.sites_per_side = max(1, math.floor(settings.cell_mm * settings.scale / self.pitch + 1e-8))
        self.cell = self.sites_per_side * self.pitch
        self.quantum_area = math.pi * self.radius ** 2
        self.reference_area = math.pi * transformer.radius ** 2 * settings.flux_gain
        self.frame_radius = frame_radius
        self.cells = defaultdict(float)
        self.milky_way = None
        self._corner_unit, self._corners = None, {}

    def add_sample(self, position, relative_flux):
        """点 ``position`` の光量を、その点を含むセルに加える。"""
        if not math.isfinite(relative_flux) or relative_flux < 0:
            raise ValueError("光量は有限の非負値が必要です。")
        pp = self.transformer.transform(position)
        if pp is None:
            return
        key = (pp.dir, pp.index, math.floor(pp.xmm / self.cell), math.floor(pp.ymm / self.cell))
        self.cells[key] += self.reference_area * relative_flux

    def add_bin(self, lon_min, lon_max, lat_min, lat_max, relative_flux):
        """銀経・銀緯の範囲（度）に一様に広がる光量を、原盤上で重なる面積に比例してセルへ分配する。

        ビンの中心1点に置くと、ビンとセルの大きさの比が整数でないため、セルごとに入る
        ビンの数が変わり（例えば4個と5個）、一様な空にも格子状の縞が出る。
        原盤上のビンの像は、各辺の両端の平均で決まる長方形とみなす（レンズ中心がドーム中心にあれば厳密）。
        """
        if not math.isfinite(relative_flux) or relative_flux < 0:
            raise ValueError("光量は有限の非負値が必要です。")
        unit = self.transformer.unit_of_galactic((lon_min + lon_max) / 2, (lat_min + lat_max) / 2)
        if unit is None:
            return
        (x00, y00), (x10, y10), (x01, y01), (x11, y11) = (
            self._plate_corner(unit, lon, lat) for lat in (lat_min, lat_max) for lon in (lon_min, lon_max))
        area = self.reference_area * relative_flux
        for ix, x_fraction in _overlaps((x00 + x10) / 2, (x01 + x11) / 2, self.cell):
            for iy, y_fraction in _overlaps((y00 + y01) / 2, (y10 + y11) / 2, self.cell):
                self.cells[(*unit, ix, iy)] += area * x_fraction * y_fraction

    def _plate_corner(self, unit, longitude, latitude):
        """ビンの角の原盤座標。隣り合うビンが共有する角は1回だけ計算する。"""
        if unit != self._corner_unit:
            self._corner_unit, self._corners = unit, {}
        key = (round(longitude * 1e7), round(latitude * 1e7))
        if key not in self._corners:
            pp = self.transformer.transform_unit(SpherePosition.from_galactic(longitude, latitude), *unit)
            self._corners[key] = pp.xmm, pp.ymm
        return self._corners[key]

    def _candidates(self, unit, ix, iy):
        n = self.sites_per_side
        for j in range(n):
            y = (iy * n + j + .5) * self.pitch
            for i in range(n):
                x = (ix * n + i + .5) * self.pitch
                # 枠外で穴が切れたり、周囲の金属幅が不足するのを避ける。
                margin = self.radius + self.settings.min_web_mm * self.settings.scale
                if math.hypot(x, y) + margin > self.frame_radius:
                    continue
                sp = self.transformer.inverse_transform_unit(x, y, *unit)
                if self.transformer.assigned_unit(sp) != unit:
                    continue
                yield x, y

    def _solid_angle_deg2(self, unit, ix, iy):
        """セルが天球上で占める立体角（平方度）。偏心投影機でもドーム中心から見た空の広さ。"""
        x, y = (ix + .5) * self.cell, (iy + .5) * self.cell
        center = self.transformer.inverse_transform_unit(x, y, *unit).to_vector(1.)
        right = self.transformer.inverse_transform_unit(x + self.cell, y, *unit).to_vector(1.)
        below = self.transformer.inverse_transform_unit(x, y + self.cell, *unit).to_vector(1.)
        normal = right.minus(center).cross(below.minus(center))
        return math.sqrt(normal.dot(normal)) * math.degrees(1.) ** 2

    def _enhance_milky_way(self):
        """背景に相当する一様な明るさを差し引き、天の川を残す倍率をセルの光量に掛ける。

        倍率は中央値で平滑化した表面輝度Sから T(S)/S とする。孤立した明るい星は中央値に
        影響せず、元のセル光量に倍率を掛けるので天の川内部の暗黒帯などの細部は残る。
        """
        s = self.settings
        half = max(0, math.ceil((s.milky_way_smoothing_mm * s.scale / self.cell - 1) / 2 - 1e-6))
        window = 2 * half + 1
        keys, smoothed = [], []
        for region in range(len(self.transformer.LONGITUDE_REGIONS)):
            unit = divmod(region, 2)
            unit_keys = [key for key in self.cells if key[:2] == unit]
            if not unit_keys:
                continue
            ix = np.array([key[2] for key in unit_keys])
            iy = np.array([key[3] for key in unit_keys])
            ix, iy = ix - ix.min(), iy - iy.min()
            # データのないセルはNaNとして中央値から除く。
            grid = np.full((ix.max() + 1 + 2 * half, iy.max() + 1 + 2 * half), np.nan)
            grid[ix + half, iy + half] = [self.cells[key] / self._solid_angle_deg2(unit, *key[2:])
                                          for key in unit_keys]
            view = sliding_window_view(grid, (window, window))
            values = np.empty(len(unit_keys))
            for start in range(0, len(unit_keys), 4096):
                part = slice(start, start + 4096)
                values[part] = np.nanmedian(view[ix[part], iy[part]], axis=(1, 2))
            keys += unit_keys
            smoothed.append(values)
        smoothed = np.concatenate(smoothed)
        # 4原盤共通の黒点・基準点。原盤ごとの自動露出はしない。
        black, anchor = np.percentile(smoothed, [s.milky_way_black_percentile,
                                                 s.milky_way_anchor_percentile])
        if not anchor > black:
            raise ValueError("天の川強調: 明るさの分布が一様で黒点と基準点が一致します。"
                             "black-percentileとanchor-percentileを見直してください。")
        target = anchor * np.clip((smoothed - black) / (anchor - black), 0, None) ** s.milky_way_gamma
        gains = np.divide(target, smoothed, out=np.zeros_like(target), where=smoothed > 0)
        cells = {key: self.cells[key] * gain for key, gain in zip(keys, gains.tolist())}
        return cells, {
            "smoothing_window_cells": window,
            # 基準等級の星が1平方度あたり何個分か。
            "black_level": float(black) / self.reference_area,
            "anchor_level": float(anchor) / self.reference_area,
            "dark_cells": int(np.count_nonzero(smoothed <= black)), "cells": len(keys),
        }

    def build(self):
        """穴と、原盤ごとの ``{plate, holes, linear_area_mm2, requested_area_mm2, capacity_loss_area_mm2}`` を返す。

        天の川強調を使った場合は、黒点・基準点などを ``self.milky_way`` に記録する。
        ``linear_area_mm2`` は強調前（G帯光量に比例）の要求面積。
        """
        holes, plates = _Holes(), []
        if self.settings.milky_way:
            cells, self.milky_way = self._enhance_milky_way()
        else:
            cells, self.milky_way = self.cells, None
        for region in range(len(self.transformer.LONGITUDE_REGIONS)):
            unit = divmod(region, 2)
            rng = random.Random(self.settings.seed + region)
            # 蛇行順の誤差繰越。各セルの誤差は最小穴1個未満、原盤全体では半個以内。
            keys = sorted((key for key in cells if key[:2] == unit),
                          key=lambda key: (key[3], key[2] if key[3] % 2 == 0 else -key[2]))
            carry = .5
            linear = requested = achievable = 0.
            hole_count = 0
            for key in keys:
                area = cells[key]
                linear += self.cells[key]
                requested += area
                candidates = list(self._candidates(unit, key[2], key[3]))
                target = min(area / self.quantum_area, len(candidates))
                achievable += target * self.quantum_area
                number = min(len(candidates), math.floor(target + carry))
                carry += target - number
                rng.shuffle(candidates)
                for x, y in candidates[:number]:
                    holes.append(x, y, region)
                hole_count += number
            plates.append({
                "plate": ("N" if unit[0] == 0 else "S") + str(unit[1]),
                "holes": hole_count, "linear_area_mm2": linear, "requested_area_mm2": requested,
                "capacity_loss_area_mm2": max(0., requested - achievable),
            })
        return holes, plates

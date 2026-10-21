"""Gaiaの局所光量を、最小径・金属幅を満たす円形穴へ量子化する。

穴は個々の恒星ではなく局所的な合計光量の標本。全穴を同じ直径にし、
固定ピッチの候補格子の一部を選ぶことで、原盤全体の最小金属幅を保証する。
天の川強調を有効にすると、背景に相当する一様な明るさを差し引いてから穴を配置する。
Copyright (C) 2026 東京大学地文研究会天文部; GPL-2.0-or-later
"""

from collections import defaultdict
from collections.abc import Sequence
from array import array
from dataclasses import asdict, dataclass
import math
import random

import numpy as np
from numpy.lib.stride_tricks import sliding_window_view

from .gaia_catalog import positive_number
from .models import PlatePosition, PlateStar


class _Holes(Sequence):
    """大量の穴は座標配列で保持し、出力時だけPlateStarを作る。"""
    def __init__(self, radius):
        self.radius = radius
        self.coordinates = array("d")

    def append(self, x, y, region):
        self.coordinates.extend((x, y, region))

    def __len__(self):
        return len(self.coordinates) // 3

    def __getitem__(self, index):
        if isinstance(index, slice):
            return [self[i] for i in range(*index.indices(len(self)))]
        if index < 0:
            index += len(self)
        if not 0 <= index < len(self):
            raise IndexError(index)
        x, y, region = self.coordinates[3 * index:3 * index + 3]
        hole = PlateStar()
        hole.p = PlatePosition()
        hole.p.dir, hole.p.index = divmod(int(region), 2)
        hole.p.xmm, hole.p.ymm, hole.rmm = x, y, self.radius
        return hole


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
        diameter = positive_number(props, "galaxy.etch.min-hole-diameter-mm", .02)
        web = positive_number(props, "galaxy.etch.min-web-mm", .02)
        cell = positive_number(props, "galaxy.etch.cell-size-mm", .2)
        gain = positive_number(props, "galaxy.etch.flux-gain", 1.)
        scale = positive_number(props, "scale", 1.)
        try:
            seed = int(props.get("galaxy.etch.seed", 0))
        except (ValueError, TypeError) as exc:
            raise ValueError("galaxy.etch.seed: 整数を指定してください。") from exc
        if cell < diameter + web:
            raise ValueError("galaxy.etch.cell-size-mm は最小穴径＋最小金属幅以上にしてください。")
        if cell / (diameter + web) > 100:
            raise ValueError("1セルの候補点が多すぎます。cell-size-mm / (穴径+金属幅) <= 100 としてください。")
        milky_way = props.get("galaxy.etch.milky-way", "no").strip().lower() in {"y", "yes", "true"}
        black = positive_number(props, "galaxy.etch.milky-way.black-percentile", 40., allow_zero=True)
        anchor = positive_number(props, "galaxy.etch.milky-way.anchor-percentile", 95.)
        gamma = positive_number(props, "galaxy.etch.milky-way.gamma", 1.)
        smoothing = positive_number(props, "galaxy.etch.milky-way.smoothing-mm", 1.)
        if not black < anchor <= 100:
            raise ValueError("galaxy.etch.milky-way: 0 <= black-percentile < anchor-percentile <= 100 としてください。")
        if smoothing / cell > 50:
            raise ValueError("galaxy.etch.milky-way.smoothing-mm は cell-size-mm の50倍以下にしてください。")
        return cls(diameter, web, cell, gain, seed, scale, milky_way, black, anchor, gamma, smoothing)


class EtchingBuilder:
    def __init__(self, transformer, settings, frame_radius):
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
        self.counts = defaultdict(int)
        self.outside_stars = 0

    def add_sample(self, position, relative_flux, count=1):
        if not math.isfinite(relative_flux) or relative_flux < 0 or count < 0:
            raise ValueError("光量と恒星数は有限の非負値が必要です。")
        pp = self.transformer.transform(position)
        if pp is None:
            self.outside_stars += count
            return
        unit = pp.dir, pp.index
        key = (*unit, math.floor(pp.xmm / self.cell), math.floor(pp.ymm / self.cell))
        self.cells[key] += self.reference_area * relative_flux
        self.counts[unit] += count

    def _candidates(self, unit, ix, iy):
        n = self.sites_per_side
        for j in range(n):
            y = (iy * n + j + .5) * self.pitch
            for i in range(n):
                x = (ix * n + i + .5) * self.pitch
                # 枠外で穴が切れたり、周囲の金属幅が不足するのを避ける。
                margin = self.radius + self.settings.min_web_mm * self.settings.scale
                if max(abs(x), abs(y)) + margin > self.frame_radius:
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
        return right.minus(center).cross(below.minus(center)).get_mag() * math.degrees(1.) ** 2

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
            "enabled": True, "smoothing_window_cells": window,
            "smoothing_window_output_mm": window * self.cell,
            # 基準等級の星が1平方度あたり何個分か。
            "black_level": float(black) / self.reference_area,
            "anchor_level": float(anchor) / self.reference_area,
            "level_unit": "reference-magnitude stars per square degree",
            "dark_cells": int(np.count_nonzero(smoothed <= black)), "cells": len(keys),
        }

    def build(self):
        holes, reports = _Holes(self.radius), []
        if self.settings.milky_way:
            cells, milky_way = self._enhance_milky_way()
        else:
            cells, milky_way = self.cells, {"enabled": False}
        for region in range(len(self.transformer.LONGITUDE_REGIONS)):
            unit = divmod(region, 2)
            label = ("N" if unit[0] == 0 else "S") + str(unit[1])
            rng = random.Random(self.settings.seed + region)
            # 蛇行順の誤差繰越。各セルの誤差は最小穴1個未満、原盤全体では半個以内。
            keys = sorted((key for key in cells if key[:2] == unit),
                          key=lambda key: (key[3], key[2] if key[3] % 2 == 0 else -key[2]))
            carry = .5
            linear = requested = achievable = 0.
            saturated_cells = empty_cells = hole_count = 0
            for key in keys:
                area = cells[key]
                linear += self.cells[key]
                requested += area
                candidates = list(self._candidates(unit, key[2], key[3]))
                demand = area / self.quantum_area
                capacity = len(candidates)
                if demand > capacity:
                    saturated_cells += 1
                if not capacity:
                    empty_cells += 1
                target = min(demand, capacity)
                achievable += target * self.quantum_area
                number = min(capacity, math.floor(target + carry))
                carry += target - number
                rng.shuffle(candidates)
                for x, y in candidates[:number]:
                    holes.append(x, y, region)
                hole_count += number
            output_area = hole_count * self.quantum_area
            reports.append({
                "plate": label, "source_stars": self.counts[unit], "cells": len(keys),
                "holes": hole_count, "linear_area_mm2": linear, "requested_area_mm2": requested,
                "achievable_area_mm2": achievable, "output_area_mm2": output_area,
                "capacity_loss_area_mm2": max(0., requested - achievable),
                "quantization_error_area_mm2": output_area - achievable,
                "saturated_cells": saturated_cells, "no_candidate_cells": empty_cells,
            })
        algorithm = "G-band flux sum; equal circular holes on spaced lattice; serpentine error carry"
        limits = [
            "光量→穴面積は均一照明・薄板を仮定。板厚、エッチング誤差、光源、ドームで要実測校正。",
            "G帯の観測光量。V帯/暗所視への変換、減光補正、欠測星・散光星雲の補完はしない。",
            "穴はセル内に再配置した光量標本で、個々の恒星位置ではない。",
            "天球の取得ビンと原盤セルの両方で空間分解能が制限される。",
        ]
        if milky_way["enabled"]:
            algorithm += "; milky-way background subtraction on median-smoothed surface brightness"
            limits.append("天の川強調: 背景を差し引いた表示用の階調。穴面積はG帯光量に比例しない（linear_area_mm2が比例時の値）。")
        report = {
            "algorithm": algorithm,
            "settings": asdict(self.settings), "actual_cell_size_output_mm": self.cell,
            "hole_diameter_output_mm": 2 * self.radius,
            "min_web_output_mm": self.settings.min_web_mm * self.settings.scale,
            "quantum_area_mm2": self.quantum_area, "outside_region_stars": self.outside_stars,
            "milky_way": milky_way,
            "plates": reports,
            "limits": limits,
        }
        return holes, report

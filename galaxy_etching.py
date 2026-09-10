"""Gaiaの局所光量を、最小径・金属幅を満たす円形穴へ量子化する。

穴は個々の恒星ではなく局所的な合計光量の標本。全穴を同じ直径にし、
固定ピッチの候補格子の一部を選ぶことで、原盤全体の最小金属幅を保証する。
Copyright (C) 2026 東京大学地文研究会天文部; GPL-2.0-or-later
"""

from collections import defaultdict
from collections.abc import Sequence
from array import array
from dataclasses import asdict, dataclass
import math
import random

from gaia_catalog import positive_number
from models import PlatePosition, PlateStar


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
        return cls(diameter, web, cell, gain, seed, scale)


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

    def build(self):
        holes, reports = _Holes(self.radius), []
        for region in range(len(self.transformer.LONGITUDE_REGIONS)):
            unit = divmod(region, 2)
            label = ("N" if unit[0] == 0 else "S") + str(unit[1])
            rng = random.Random(self.settings.seed + region)
            # 蛇行順の誤差繰越。各セルの誤差は最小穴1個未満、原盤全体では半個以内。
            keys = sorted((key for key in self.cells if key[:2] == unit),
                          key=lambda key: (key[3], key[2] if key[3] % 2 == 0 else -key[2]))
            carry = .5
            requested = achievable = 0.
            saturated_cells = empty_cells = hole_count = 0
            for key in keys:
                area = self.cells[key]
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
                "holes": hole_count, "requested_area_mm2": requested,
                "achievable_area_mm2": achievable, "output_area_mm2": output_area,
                "capacity_loss_area_mm2": max(0., requested - achievable),
                "quantization_error_area_mm2": output_area - achievable,
                "saturated_cells": saturated_cells, "no_candidate_cells": empty_cells,
            })
        report = {
            "algorithm": "G-band flux sum; equal circular holes on spaced lattice; serpentine error carry",
            "settings": asdict(self.settings), "actual_cell_size_output_mm": self.cell,
            "hole_diameter_output_mm": 2 * self.radius,
            "min_web_output_mm": self.settings.min_web_mm * self.settings.scale,
            "quantum_area_mm2": self.quantum_area, "outside_region_stars": self.outside_stars,
            "plates": reports,
            "limits": [
                "光量→穴面積は均一照明・薄板を仮定。板厚、エッチング誤差、光源、ドームで要実測校正。",
                "G帯の観測光量。V帯/暗所視への変換、減光補正、欠測星・散光星雲の補完はしない。",
                "穴はセル内に再配置した光量標本で、個々の恒星位置ではない。",
                "天球の取得ビンと原盤セルの両方で空間分解能が制限される。",
            ],
        }
        return holes, report

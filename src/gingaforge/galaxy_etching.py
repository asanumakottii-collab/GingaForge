"""Gaiaの局所光量を、最小径・金属幅を満たす円形穴へ量子化する。

穴は個々の恒星ではなく局所的な合計光量の標本。全穴を同じ直径にし、
固定ピッチの候補格子の一部を選ぶことで、原盤全体の最小金属幅を保証する。
Copyright (C) 2026 東京大学地文研究会天文部; GPL-2.0-or-later
"""

from collections import defaultdict
from array import array
from dataclasses import dataclass
import math
import random

from .config import positive_number


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

    def add_sample(self, position, relative_flux):
        if not math.isfinite(relative_flux) or relative_flux < 0:
            raise ValueError("光量は有限の非負値が必要です。")
        pp = self.transformer.transform(position)
        if pp is None:
            return
        key = (pp.dir, pp.index, math.floor(pp.xmm / self.cell), math.floor(pp.ymm / self.cell))
        self.cells[key] += self.reference_area * relative_flux

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
        """穴と、原盤ごとの ``{plate, holes, requested_area_mm2, capacity_loss_area_mm2}`` を返す。"""
        holes, plates = _Holes(), []
        for region in range(len(self.transformer.LONGITUDE_REGIONS)):
            unit = divmod(region, 2)
            rng = random.Random(self.settings.seed + region)
            # 蛇行順の誤差繰越。各セルの誤差は最小穴1個未満、原盤全体では半個以内。
            keys = sorted((key for key in self.cells if key[:2] == unit),
                          key=lambda key: (key[3], key[2] if key[3] % 2 == 0 else -key[2]))
            carry = .5
            requested = achievable = 0.
            hole_count = 0
            for key in keys:
                area = self.cells[key]
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
                "holes": hole_count, "requested_area_mm2": requested,
                "capacity_loss_area_mm2": max(0., requested - achievable),
            })
        return holes, plates

# -*- coding: utf-8 -*-
#
# Orb Transform Library(OTL)
# Copyright (C) 2007,2012,2014 東京大学地文研究会天文部
# Copyright (C) 2026 東京大学地文研究会天文部
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
#
"""OHP・Gaia原盤で共有する担当領域と円筒投影。"""

import math
from .mathvector import MathVector
from .models import PlatePosition, SpherePosition

REGIONS = {"N0": (0., 60.), "N1": (60., 120.),
           "S0": (300., 360.), "S1": (240., 300.)}
LATITUDE_LIMIT = 20.


def lens_distance(image, focus):
    """原板とレンズの距離。

    原板から像（ドーム）までの距離 image と焦点距離 focus から、OTLと同じくレンズの公式
    1/focus = 1/x + 1/(image - x) の小さい方の解 x を返す。
    """
    discriminant = image * image / 4 - image * focus
    if not (focus > 0 and discriminant >= 0):
        raise ValueError("焦点距離が長すぎて、ドームに像を結べません。"
                         "焦点距離を、原板からドームまでの距離の1/4以下にしてください。")
    # 2解の積が image * focus であることを使い、桁落ちを避けて小さい方の解を求める。
    return image * focus / (image / 2 + math.sqrt(discriminant))


class GalaxyTransformer:
    # (N0, N1, S0, S1) の順。N/Sは出力用紙上の位置を表し、銀緯の正負ではない。
    LONGITUDE_REGIONS = tuple(REGIONS.values())
    LATITUDE_LIMIT = LATITUDE_LIMIT
    _ANGLE_EPS = 1e-10  # 度。座標回転の丸め誤差を吸収する。

    def __init__(self, radius, sphere, projector_horizontal, projector_vertical, plate, focus):
        """
        座標変換のパラメータを設定します。
        :param radius: 7.5等星の半径
        :param sphere: ドームの半径
        :param projector_horizontal: ドームの中心から投影機中心へのベクトルの赤道面に平行な成分
        :param projector_vertical: ドームの中心から投影機中心へのベクトルの赤道面に垂直な成分
        :param plate: 投影機中心と原版中心の間の距離
        :param focus: 焦点の距離
        """
        self.radius = radius
        self.sphere = sphere
        self.projector_vertical = projector_vertical
        self.projector_horizontal = projector_horizontal
        self.plate = plate
        self.focus = focus

        # 投影機の位置は従来どおり赤道座標の(0, ±水平距離, ±垂直距離)。
        # 各担当領域の銀経中央・銀緯0度への光線を原盤中心に合わせる。
        # OTLと同じく、原板は投影機中心から plate の位置にあり、レンズはその先に、
        # 原板とレンズの距離がレンズの公式を満たす位置にある。円筒投影の半径は原板とレンズの距離。
        galactic_north = SpherePosition.from_galactic(0., 90.).to_vector(1.)
        self._plate_geometry = {}
        for region, (start, end) in enumerate(self.LONGITUDE_REGIONS):
            dir_, index = divmod(region, 2)
            projector = MathVector(
                0., projector_horizontal * (1 if index == 0 else -1),
                projector_vertical * (1 if dir_ == 0 else -1),
            )
            center = SpherePosition.from_galactic((start + end) / 2., 0.)
            to_center = center.to_vector(sphere).minus(projector)
            forward = to_center.unit_vector()
            east = galactic_north.cross(forward).unit_vector()
            north = forward.cross(east).unit_vector()
            lens = lens_distance(math.sqrt(to_center.dot(to_center)) - plate, focus)
            lens_center = projector.plus(forward.mult_scalar(plate + lens))
            self._plate_geometry[dir_, index] = (lens_center, lens, forward, east, north)

    @classmethod
    def _normalize_longitude(cls, longitude):
        longitude %= 360.
        for boundary in (0., 60., 120., 180., 240., 300., 360.):
            if abs(longitude - boundary) <= cls._ANGLE_EPS:
                return boundary % 360.
        return longitude

    def assigned_unit(self, sp):
        """銀緯±20度の担当 ``(出力位置, 番号)`` を返す。対象外はNone。"""
        return self.unit_of_galactic(*sp.to_galactic())

    def unit_of_galactic(self, longitude, latitude):
        """銀経・銀緯（度）の担当 ``(出力位置, 番号)`` を返す。対象外はNone。"""
        if not abs(latitude) <= self.LATITUDE_LIMIT + self._ANGLE_EPS:
            return None
        longitude = self._normalize_longitude(longitude)
        for region, (start, end) in enumerate(self.LONGITUDE_REGIONS):
            if start <= longitude < end:
                return divmod(region, 2)
        return None

    def transform(self, sp):
        """銀河座標で担当原盤を選び、領域外の位置は除外します。"""
        unit = self.assigned_unit(sp)
        if unit is None:
            return None
        return self.transform_unit(sp, *unit)

    def lens_distance(self, dir_, index):
        """原盤 ``(出力位置, 番号)`` の原板とレンズの距離（円筒投影の半径）。"""
        return self._plate_geometry[dir_, index][1]

    def transform_unit(self, sp, dir_, index):
        """担当領域の中央を基準に、レンズ中心からの円筒投影で原盤座標へ変換します。"""
        lens_center, lens, forward, east, north = self._plate_geometry[dir_, index]
        ray = sp.to_vector(self.sphere).minus(lens_center)
        x, y, z = ray.dot(forward), ray.dot(east), ray.dot(north)
        pp = PlatePosition()
        pp.index = index
        pp.dir = dir_
        pp.xmm = lens * z / math.hypot(x, y)
        pp.ymm = -lens * math.atan2(y, x)
        return pp

    def inverse_transform_unit(self, xmm, ymm, dir_, index):
        """原盤上の点から、同じ円筒投影の光線とドームの交点を求める。"""
        lens_center, lens, forward, east, north = self._plate_geometry[dir_, index]
        angle = -ymm / lens
        ray = (forward.mult_scalar(math.cos(angle))
               .plus(east.mult_scalar(math.sin(angle)))
               .plus(north.mult_scalar(xmm / lens))).unit_vector()
        dot = lens_center.dot(ray)
        discriminant = dot * dot + self.sphere ** 2 - lens_center.dot(lens_center)
        if discriminant < 0:
            raise ValueError("原盤の光線がドームと交差しません。")
        distance = -dot + math.sqrt(discriminant)
        return SpherePosition.from_vector(lens_center.plus(ray.mult_scalar(distance)))

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
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.
#
"""OHP・Gaia原盤で共有する担当領域と円筒投影。"""

import math
from mathvector import MathVector
from models import PlatePosition, SpherePosition

REGIONS = {"N0": (0., 60.), "N1": (60., 120.),
           "S0": (300., 360.), "S1": (240., 300.)}
LATITUDE_LIMIT = 20.


class GalaxyTransformer:
    # (N0, N1, S0, S1) の順。N/Sは出力用紙上の位置を表し、銀緯の正負ではない。
    LONGITUDE_REGIONS = tuple(REGIONS.values())
    LATITUDE_LIMIT = LATITUDE_LIMIT
    _ANGLE_EPS = 1e-10  # 度。座標回転の丸め誤差を吸収する。

    def __init__(self, radius, sphere, projector_horizontal, projector_vertical, plate):
        """
        座標変換のパラメータを設定します。
        :param radius: 7.5等星の半径
        :param sphere: ドームの半径
        :param projector_horizontal: ドームの中心から投影機中心へのベクトルの赤道面に平行な成分
        :param projector_vertical: ドームの中心から投影機中心へのベクトルの赤道面に垂直な成分
        :param plate: 投影機中心と原版中心の間の距離
        """
        self.radius = radius
        self.sphere = sphere
        self.projector_vertical = projector_vertical
        self.projector_horizontal = projector_horizontal
        self.plate = plate

        # 投影機の位置は従来どおり赤道座標の(0, ±水平距離, ±垂直距離)。
        # 各担当領域の銀経中央・銀緯0度への光線を原盤中心に合わせる。
        galactic_north = SpherePosition.from_galactic(0., 90.).to_vector(1.)
        self._plate_geometry = {}
        for region, (start, end) in enumerate(self.LONGITUDE_REGIONS):
            dir_, index = divmod(region, 2)
            projector = MathVector(
                0., projector_horizontal * (1 if index == 0 else -1),
                projector_vertical * (1 if dir_ == 0 else -1),
            )
            center = SpherePosition.from_galactic((start + end) / 2., 0.)
            forward = center.to_vector(sphere).minus(projector).unit_vector()
            east = galactic_north.cross(forward).unit_vector()
            north = forward.cross(east).unit_vector()
            self._plate_geometry[dir_, index] = (projector, forward, east, north)

    @classmethod
    def _normalize_longitude(cls, longitude):
        longitude %= 360.
        for boundary in (0., 60., 120., 180., 240., 300., 360.):
            if abs(longitude - boundary) <= cls._ANGLE_EPS:
                return boundary % 360.
        return longitude

    def assigned_unit(self, sp):
        """銀緯±20度の担当 ``(出力位置, 番号)`` を返す。対象外はNone。"""
        longitude, latitude = sp.to_galactic()
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

    def transform_unit(self, sp, dir_, index):
        """担当領域の中央を基準に、従来の円筒投影で原盤座標へ変換します。"""
        projector, forward, east, north = self._plate_geometry[dir_, index]
        ray = sp.to_vector(self.sphere).minus(projector)
        x, y, z = ray.dot(forward), ray.dot(east), ray.dot(north)
        pp = PlatePosition()
        pp.index = index
        pp.dir = dir_
        pp.xmm = self.plate * z / math.hypot(x, y)
        pp.ymm = -self.plate * math.atan2(y, x)
        return pp

    def inverse_transform_unit(self, xmm, ymm, dir_, index):
        """原盤上の点から、同じ円筒投影の光線とドームの交点を求める。"""
        projector, forward, east, north = self._plate_geometry[dir_, index]
        angle = -ymm / self.plate
        ray = (forward.mult_scalar(math.cos(angle))
               .plus(east.mult_scalar(math.sin(angle)))
               .plus(north.mult_scalar(xmm / self.plate))).unit_vector()
        dot = projector.dot(ray)
        discriminant = dot * dot + self.sphere ** 2 - projector.dot(projector)
        if discriminant < 0:
            raise ValueError("原盤の光線がドームと交差しません。")
        distance = -dot + math.sqrt(discriminant)
        return SpherePosition.from_vector(projector.plus(ray.mult_scalar(distance)))

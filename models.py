# -*- coding: utf-8 -*-
#
# Orb Transform Library(OTL)
# Copyright (C) 2007,2012 東京大学地文研究会天文部
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
"""天球・原盤上の位置とエッチング穴を表すデータクラス群。"""

import math

from mathvector import MathVector


# Hipparcos/ICRSから銀河座標への回転行列（各行は銀河座標軸）。
# 参照: https://github.com/liberfa/erfa/blob/master/src/icrs2g.c
_GALACTIC_AXES = (
    MathVector(-0.05487556041621537, -0.8734370902348850, -0.4838350155487132),
    MathVector(0.4941094278755837, -0.4448296299600112, 0.7469822444972189),
    MathVector(-0.8676661490190047, -0.1980763734312015, 0.4559837761750669),
)


class PlatePosition:
    """原板上での位置を扱うクラスです。"""

    def __init__(self):
        self.xmm = 0.0
        self.ymm = 0.0
        self.dir = 0  # 北天(0)，南天(1)の別
        self.index = 0  # ユニット番号

    def __eq__(self, other):
        if not isinstance(other, PlatePosition):
            return NotImplemented
        # NOTE: 原Javaコードの isSamePlate 判定にあった dir と p.index の
        # 比較をそのまま踏襲している(挙動を変えないための忠実な移植)。
        return (self.xmm == other.xmm and self.ymm == other.ymm
                and self.dir == other.index and self.index == other.index)

    def is_same_plate(self, other):
        """同じ原板上にあるかどうか"""
        return self.dir == other.dir and self.index == other.index

    def distance_from_center(self):
        """原板中心からの距離"""
        return math.sqrt(self.xmm * self.xmm + self.ymm * self.ymm)


class SpherePosition:
    """赤道座標を扱うクラスです。"""

    def __init__(self):
        self.dedeg = 0.0  # 赤緯
        self.radeg = 0.0  # 赤経

    def __eq__(self, other):
        if not isinstance(other, SpherePosition):
            return NotImplemented
        return self.dedeg == other.dedeg and self.radeg == other.radeg

    def to_vector(self, radius):
        """直交座標に変換します。"""
        return MathVector.from_mag_lng_lat(radius, math.radians(self.radeg), math.radians(self.dedeg))

    def to_galactic(self):
        """ICRS赤道座標から ``(銀経[0, 360), 銀緯)`` を度で返します。"""
        vector = self.to_vector(1.).transform(_GALACTIC_AXES)
        longitude = math.degrees(vector.get_lng()) % 360.
        return longitude, math.degrees(vector.get_lat())

    @staticmethod
    def from_galactic(longitude, latitude):
        """度単位の銀経・銀緯からICRS赤道座標を作成します。"""
        vector = MathVector.from_mag_lng_lat(
            1., math.radians(longitude), math.radians(latitude)
        )
        equatorial = (_GALACTIC_AXES[0].mult_scalar(vector.x)
                      .plus(_GALACTIC_AXES[1].mult_scalar(vector.y))
                      .plus(_GALACTIC_AXES[2].mult_scalar(vector.z)))
        return SpherePosition.from_vector(equatorial)

    @staticmethod
    def from_vector(vector):
        """直交座標から新しいSpherePositionインスタンスを作成します。"""
        p = SpherePosition()
        p.radeg = math.degrees(vector.get_lng())
        p.dedeg = math.degrees(vector.get_lat())
        return p


class PlateStar:
    """原板上の恒星を扱うクラスです。"""

    def __init__(self):
        self.rmm = 0.0  # 原板上での半径
        self.p = None  # 原板上での位置 (PlatePosition)
        self.hip_number = None  # Hipparcos HIP番号 (Tycho星などはNone)

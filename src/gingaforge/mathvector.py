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
"""数学ベクトルのデータセットとその計算を扱うモジュールです。(MathVector.javaの移植)"""

import math

import numpy as np


class MathVector:
    """数学ベクトル。"""

    __slots__ = ("_v",)

    def __init__(self, x=0.0, y=0.0, z=0.0):
        self._v = np.array([x, y, z], dtype=float)

    @staticmethod
    def from_mag_lng_lat(mag, lng, lat):
        """極座標から新しいMathVectorオブジェクトを作成します。"""
        return MathVector(mag * math.cos(lng) * math.cos(lat),
                           mag * math.sin(lng) * math.cos(lat),
                           mag * math.sin(lat))

    @staticmethod
    def _of(v):
        result = MathVector()
        result._v = v
        return result

    @property
    def x(self):
        return float(self._v[0])

    @property
    def y(self):
        return float(self._v[1])

    @property
    def z(self):
        return float(self._v[2])

    def get_lng(self):
        """極座標の経度を取得します。"""
        return math.atan2(self.y, self.x)

    def get_lat(self):
        """極座標の緯度を取得します。"""
        return math.atan2(self.z, math.hypot(self.x, self.y))

    def plus(self, operand):
        """このベクトルに指定されたベクトルを加算したベクトルを作成します。"""
        return self._of(self._v + operand._v)

    def minus(self, operand):
        """このベクトルに指定されたベクトルを減算したベクトルを作成します。"""
        return self._of(self._v - operand._v)

    def mult_scalar(self, operand):
        """このベクトルを指定されたスカラー倍したベクトルを作成します。"""
        return self._of(self._v * operand)

    def unit_vector(self):
        """このベクトルの単位ベクトルを作成します。"""
        return self.mult_scalar(1.0 / float(np.linalg.norm(self._v)))

    def dot(self, operand):
        """このベクトルの内積を計算します。"""
        return float(np.dot(self._v, operand._v))

    def cross(self, operand):
        """このベクトルの外積を計算します。"""
        return self._of(np.cross(self._v, operand._v))

    def transform(self, v):
        """指定されたベクトルの方向をx, y, z方向としたベクトルを作成します。"""
        return MathVector(self.dot(v[0]), self.dot(v[1]), self.dot(v[2]))

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
    """数学ベクトル。x, y, z 成分は可変（原コードが直接書き換えるため）。"""

    __slots__ = ("_v",)

    def __init__(self, x=0.0, y=0.0, z=0.0):
        if isinstance(x, MathVector):
            self._v = x._v.copy()
        else:
            self._v = np.array([x, y, z], dtype=float)

    @staticmethod
    def from_mag_lng_lat(mag, lng, lat):
        """極座標から新しいMathVectorオブジェクトを作成します。"""
        return MathVector(mag * math.cos(lng) * math.cos(lat),
                           mag * math.sin(lng) * math.cos(lat),
                           mag * math.sin(lat))

    @property
    def x(self):
        return float(self._v[0])

    @x.setter
    def x(self, value):
        self._v[0] = value

    @property
    def y(self):
        return float(self._v[1])

    @y.setter
    def y(self, value):
        self._v[1] = value

    @property
    def z(self):
        return float(self._v[2])

    @z.setter
    def z(self, value):
        self._v[2] = value

    def get_mag(self):
        """極座標のrを取得します。"""
        return float(np.linalg.norm(self._v))

    def get_lng(self):
        """極座標の経度を取得します。"""
        return math.atan2(self.y, self.x)

    def get_lat(self):
        """極座標の緯度を取得します。"""
        return math.atan2(self.z, math.hypot(self.x, self.y))

    def __eq__(self, other):
        if not isinstance(other, MathVector):
            return NotImplemented
        return bool(np.array_equal(self._v, other._v))

    def plus(self, operand):
        """このベクトルに指定されたベクトルを加算したベクトルを作成します。"""
        result = MathVector()
        result._v = self._v + operand._v
        return result

    def minus(self, operand):
        """このベクトルに指定されたベクトルを減算したベクトルを作成します。"""
        result = MathVector()
        result._v = self._v - operand._v
        return result

    def mult_scalar(self, operand):
        """このベクトルを指定されたスカラー倍したベクトルを作成します。"""
        result = MathVector()
        result._v = self._v * operand
        return result

    def unit_vector(self):
        """このベクトルの単位ベクトルを作成します。"""
        return self.mult_scalar(1.0 / self.get_mag())

    def dot(self, operand):
        """このベクトルの内積を計算します。"""
        return float(np.dot(self._v, operand._v))

    def cross(self, operand):
        """このベクトルの外積を計算します。"""
        result = MathVector()
        result._v = np.cross(self._v, operand._v)
        return result

    def ang(self, v):
        """このベクトルと指定されたベクトルのなす角度を計算します。"""
        cos_theta = self.dot(v) / self.get_mag() / v.get_mag()
        try:
            return math.acos(cos_theta)
        except ValueError:
            # Java の Math.acos は定義域外の入力に対して例外を投げず NaN を返す。
            # NaN を用いた比較(<, >, <=, >=)は常にFalseになる点も含めて挙動を揃える。
            return math.nan

    def project(self, v):
        """指定されたベクトルの方向をz方向としたときに、このベクトルをxy平面に投影したベクトルを作成します。"""
        return self.minus(v.mult_scalar(self.dot(v) / v.dot(v)))

    def reflect(self, v=None):
        """
        v が指定されている場合は、v の方向をz方向としたときにxy平面に関して対称に配置したベクトルを、指定されていない場合は原点に関して対称に配置したベクトルを作成します。
        """
        if v is None:
            result = MathVector()
            result._v = -self._v
            return result
        return self.minus(v.mult_scalar(self.dot(v) * 2 / v.dot(v)))

    def transform(self, v):
        """指定されたベクトルの方向をx, y, z方向としたベクトルを作成します。"""
        return MathVector(self.dot(v[0]), self.dot(v[1]), self.dot(v[2]))

    def __repr__(self):
        return f"[{self.x}, {self.y}, {self.z}]"

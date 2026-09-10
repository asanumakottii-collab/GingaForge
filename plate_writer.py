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
"""原板データをSVGまたは印刷用PDFとして書き出すモジュールです。"""

import math
import os
import tempfile
from abc import ABC, abstractmethod
from enum import Enum

DEFAULT_OUTPUT_DIR = "output"

# SVGのpoints属性にはmmなどの単位を付けられないため、絶対長で記述した円・文字と座標を揃える際にCSS標準の96 dpiへ変換する。
_CSS_PIXELS_PER_MM = 96. / 25.4


def resolve_output_path(output_dir, filename):
    """出力フォルダを作成し、ファイルの出力パスを返します。"""
    output_dir = (output_dir or DEFAULT_OUTPUT_DIR).strip() or DEFAULT_OUTPUT_DIR
    if output_dir == "." or os.path.isabs(filename):
        path = filename
    else:
        path = os.path.join(output_dir, filename)
    dirname = os.path.dirname(path)
    if dirname:
        os.makedirs(dirname, exist_ok=True)
    return path


def categorized_output_dir(output_dir, category):
    """出力ルート配下の分類別フォルダを返します。"""
    output_dir = (output_dir or DEFAULT_OUTPUT_DIR).strip() or DEFAULT_OUTPUT_DIR
    return os.path.join(output_dir, category)


class PlateWriterType(Enum):
    SVG = "SVG"
    PDF = "PDF"


class PlateWriter(ABC):
    @abstractmethod
    def write_star(self, s):
        raise NotImplementedError


    @abstractmethod
    def close(self):
        raise NotImplementedError

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class _PlateWriterBase(PlateWriter):
    """PlateWriterPDF と PlateWriterSVG に共通するレイアウト計算。"""

    WIDTH = 210.  # ISO A4用紙の幅
    HEIGHT = 297. / 2.  # ISO A4用紙の高さの半分（北天・南天を上下に配置）

    def __init__(self, column, row, r, shape, filename_prefix, invert_color, output_dir=DEFAULT_OUTPUT_DIR):
        if column <= 0:
            raise ValueError("横の値が不正です。")
        if row <= 0:
            raise ValueError("縦の値が不正です。")
        self.outs = []
        self.number_of_plate = 0  # 出力用紙1ページ中の原板の数
        self.number_of_page = 0  # 出力用紙のページ数
        self.column = column
        self.row = row
        self.r = (min(self.WIDTH / column, self.HEIGHT / row) - 10.) / 2 if r == 0 else r
        self.shape = shape  # 原板が円形ならTrue, 原板が長方形ならFalse
        self.filename_prefix = filename_prefix
        self.invert_color = invert_color  # Falseなら原板を白で星を黒, Trueなら原板を黒で星を白
        self.output_dir = output_dir

    def _get_cx(self, index):
        """枠の中心のx座標"""
        return self.WIDTH / (self.column * 2) * ((index % (self.column * self.row)) // self.row * 2 + 1)

    def _get_cy(self, dir_, index):
        """枠の中心のy座標"""
        return self.HEIGHT / (self.row * 2) * ((index % self.row) * 2 + 1) + self.HEIGHT * dir_

    def _output_path(self, extension):
        return resolve_output_path(self.output_dir, f"{self.filename_prefix}{self.number_of_page}.{extension}")

    def write_frames(self, number_of_units):
        """上下それぞれに指定個数の原盤枠を、星がない場合も出力します。"""
        if self.outs is None:
            raise IOError("writer is closed")
        if number_of_units > 0:
            self._write_frame(number_of_units - 1)

    def _is_position_in_frame(self, x, y):
        """原盤上の中心座標が、実際に描画する枠の内側かを返します。"""
        if not math.isfinite(x) or not math.isfinite(y):
            return False
        if self.shape:
            return math.hypot(x, y) <= self.r
        return abs(x) <= self.r and abs(y) <= self.r

    def __del__(self):
        if getattr(self, "outs", None) is not None:
            try:
                self.close()
            except Exception:
                pass


class PlateWriterPDF(_PlateWriterBase):
    """ISO A4（210 x 297 mm）の1ページPDFをベクトルで書き出します。

    原盤番号と穴を出力します。
    円形原盤の左上にはSVGと同じ向きマークを出力します。
    """

    def __init__(self, column, row, r, shape, filename_prefix, invert_color,
                 output_dir=DEFAULT_OUTPUT_DIR):
        try:
            from reportlab.lib.units import mm
            from reportlab.pdfgen.canvas import Canvas
        except ImportError as exc:
            raise RuntimeError(
                "PDF出力にはReportLabが必要です。"
                "'pip install -r requirements.txt' を実行してください。"
            ) from exc

        super().__init__(column, row, r, shape, filename_prefix, invert_color, output_dir)
        self._mm = mm
        self._canvas_class = Canvas
        self._page_size = (self.WIDTH * mm, self.HEIGHT * 2 * mm)

    def _x(self, value):
        """mm単位のx座標をPDFのpoint単位に変換します。"""
        return value * self._mm

    def _y(self, value):
        """y軸下向きのmm座標をPDFのy軸上向き座標に変換します。"""
        return (self.HEIGHT * 2 - value) * self._mm

    def _write_frame(self, index):
        page = index // (self.column * self.row)
        while page >= self.number_of_page:
            out = self._canvas_class(
                self._output_path("pdf"),
                pagesize=self._page_size,
                pageCompression=1,
                invariant=1,
            )
            self.outs.append(out)
            self.number_of_page += 1
        while index >= self.number_of_plate:
            out = self.outs[self.number_of_plate // (self.column * self.row)]
            cx = self._get_cx(self.number_of_plate)
            for i in range(2):
                cy = self._get_cy(i, self.number_of_plate)
                label = ("N" if i == 0 else "S") + str(self.number_of_plate)
                out.setFillGray(0)
                out.setFont("Helvetica", self._x(3))
                out.drawCentredString(self._x(cx), self._y(cy - self.r), label)
                out.setLineWidth(self._x(0.1))
                out.setStrokeGray(0)
                if self.shape:
                    out.setFillGray(0)
                    out.circle(
                        self._x(cx), self._y(cy), self._x(self.r),
                        stroke=0 if self.invert_color else 1,
                        fill=1 if self.invert_color else 0,
                    )
                    out.line(
                        self._x(cx - self.r * 3. / 4.), self._y(cy - self.r * 3. / 4.),
                        self._x(cx - self.r), self._y(cy - self.r),
                    )
                    out.line(
                        self._x(cx - self.r * 3. / 4.), self._y(cy - self.r),
                        self._x(cx - self.r), self._y(cy - self.r * 3. / 4.),
                    )
                else:
                    x, y = cx - self.r, cy - self.r
                    w, h = self.r * 2, self.r * 2
                    out.setFillGray(0)
                    out.rect(
                        self._x(x), self._y(y + h), self._x(w), self._x(h),
                        stroke=0 if self.invert_color else 1,
                        fill=1 if self.invert_color else 0,
                    )
            self.number_of_plate += 1

    def write_star(self, s):
        if self.outs is None:
            raise IOError("writer is closed")
        self._write_frame(s.p.index)
        if (not self._is_position_in_frame(s.p.xmm, s.p.ymm)
                or not math.isfinite(s.rmm) or s.rmm <= 0):
            return
        out = self.outs[s.p.index // (self.column * self.row)]
        cx = self._get_cx(s.p.index)
        cy = self._get_cy(s.p.dir, s.p.index)
        out.setFillGray(1 if self.invert_color else 0)
        out.circle(
            self._x(cx + s.p.xmm),
            self._y(cy + s.p.ymm),
            self._x(s.rmm),
            stroke=0,
            fill=1,
        )


    def close(self):
        if self.outs is None:
            raise IOError("writer is closed")
        outs = self.outs
        # 保存失敗後にデストラクタが同じCanvasを再保存しない。
        self.outs = None
        for index, out in enumerate(outs):
            destination = resolve_output_path(
                self.output_dir, f"{self.filename_prefix}{index}.pdf")
            temporary = None
            try:
                # 既存ファイルへの直接書き込みは同期ストレージ等でタイムアウト
                # することがある。同じディレクトリに完成させてから原子的に置換。
                with tempfile.NamedTemporaryFile(
                        mode="wb", dir=os.path.dirname(os.path.abspath(destination)),
                        prefix=".otl-pdf-", suffix=".tmp", delete=False) as stream:
                    temporary = stream.name
                    out.showPage()
                    stream.write(out.getpdfdata())
                os.replace(temporary, destination)
            except OSError as exc:
                raise OSError(
                    exc.errno, f"PDFの保存に失敗しました: {exc.strerror}", destination
                ) from exc
            finally:
                if temporary is not None and os.path.exists(temporary):
                    os.unlink(temporary)


class PlateWriterSVG(_PlateWriterBase):
    def _write_frame(self, index):
        page = index // (self.column * self.row)
        while page >= self.number_of_page:
            out = open(self._output_path("svg"), "w", encoding="utf-8")
            out.write("<?xml version=\"1.0\"?>")
            out.write("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
                      "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">")
            out.write(f"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{self.WIDTH}mm\" "
                      f"height=\"{self.HEIGHT * 2}mm\" version=\"1.1\">")
            self.outs.append(out)
            self.number_of_page += 1
        while index >= self.number_of_plate:
            out = self.outs[self.number_of_plate // (self.column * self.row)]
            cx = self._get_cx(self.number_of_plate)
            for i in range(2):
                cy = self._get_cy(i, self.number_of_plate)
                label = ("N" if i == 0 else "S") + str(self.number_of_plate)
                out.write(f"<text x=\"{cx}mm\" y=\"{cy - self.r}mm\" text-anchor=\"middle\" "
                          f"font-size=\"3mm\" font-familiy=\"Verdana\">{label}</text>")
                if self.shape:
                    out.write(f"<circle cx=\"{cx}mm\" cy=\"{cy}mm\" r=\"{self.r}mm\" "
                              f"fill=\"{'black' if self.invert_color else 'none'}\" "
                              f"stroke=\"{'none' if self.invert_color else 'black'}\" "
                              f"stroke-width=\"0.1mm\" />")
                    out.write(f"<line x1=\"{cx - self.r * 3. / 4.}mm\" y1=\"{cy - self.r * 3. / 4.}mm\" "
                              f"x2=\"{cx - self.r}mm\" y2=\"{cy - self.r}mm\" "
                              f"stroke=\"black\" stroke-width=\"0.1mm\" />")
                    out.write(f"<line x1=\"{cx - self.r * 3. / 4.}mm\" y1=\"{cy - self.r}mm\" "
                              f"x2=\"{cx - self.r}mm\" y2=\"{cy - self.r * 3. / 4.}mm\" "
                              f"stroke=\"black\" stroke-width=\"0.1mm\" />")
                else:
                    out.write(f"<rect x=\"{cx - self.r}mm\" y=\"{cy - self.r}mm\" "
                              f"width=\"{self.r * 2}mm\" height=\"{self.r * 2}mm\" "
                              f"fill=\"{'black' if self.invert_color else 'none'}\" "
                              f"stroke=\"{'white' if self.invert_color else 'black'}\" "
                              f"stroke-width=\"0.1mm\" />")
            self.number_of_plate += 1

    def write_star(self, s):
        if self.outs is None:
            raise IOError("writer is closed")
        self._write_frame(s.p.index)
        if not self._is_position_in_frame(s.p.xmm, s.p.ymm):
            return
        out = self.outs[s.p.index // (self.column * self.row)]
        cx = self._get_cx(s.p.index)
        cy = self._get_cy(s.p.dir, s.p.index)
        out.write(f"<circle cx=\"{cx + s.p.xmm}mm\" cy=\"{cy + s.p.ymm}mm\" r=\"{s.rmm}mm\" "
                  f"fill=\"{'white' if self.invert_color else 'black'}\"/>")


    def close(self):
        if self.outs is None:
            raise IOError("writer is closed")
        for out in self.outs:
            out.write("</svg>")
            out.close()
        self.outs = None

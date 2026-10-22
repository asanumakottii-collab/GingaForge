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
"""原板データをSVGまたは印刷用ベクトルPDFとして書き出すモジュールです。"""

import math
import os
import tempfile


class _PlateWriterBase:
    """PlateWriterPDF と PlateWriterSVG に共通するレイアウト計算。"""

    WIDTH = 210.  # ISO A4用紙の幅
    HEIGHT = 297. / 2.  # ISO A4用紙の高さの半分（上段・下段に配置）

    def __init__(self, column, row, r, filename_prefix, invert_color, output_dir):
        if column <= 0:
            raise ValueError("横の値が不正です。")
        if row <= 0:
            raise ValueError("縦の値が不正です。")
        self.outs = []
        self.number_of_plate = 0  # 出力済みの原板の数
        self.column = column
        self.row = row
        self.r = (min(self.WIDTH / column, self.HEIGHT / row) - 10.) / 2 if r == 0 else r
        self.filename_prefix = filename_prefix
        self.invert_color = invert_color  # Falseなら原板を白で星を黒, Trueなら原板を黒で星を白
        self.output_dir = output_dir

    def _get_cx(self, index):
        """枠の中心のx座標"""
        return self.WIDTH / (self.column * 2) * ((index % (self.column * self.row)) // self.row * 2 + 1)

    def _get_cy(self, dir_, index):
        """枠の中心のy座標"""
        return self.HEIGHT / (self.row * 2) * ((index % self.row) * 2 + 1) + self.HEIGHT * dir_

    def _output_path(self, page, extension):
        os.makedirs(self.output_dir, exist_ok=True)
        return os.path.join(self.output_dir, f"{self.filename_prefix}{page}.{extension}")

    def _write_frame(self, index):
        """指定番号までの原盤枠を上下とも出力し、必要な用紙を開きます。"""
        while index // (self.column * self.row) >= len(self.outs):
            self.outs.append(self._open_page(len(self.outs)))
        while index >= self.number_of_plate:
            out = self.outs[self.number_of_plate // (self.column * self.row)]
            cx = self._get_cx(self.number_of_plate)
            for i in range(2):
                label = ("N" if i == 0 else "S") + str(self.number_of_plate)
                self._draw_frame(out, cx, self._get_cy(i, self.number_of_plate), label)
            self.number_of_plate += 1

    def write_frames(self, number_of_units):
        """上下それぞれに指定個数の原盤枠を、星がない場合も出力します。"""
        if number_of_units > 0:
            self._write_frame(number_of_units - 1)

    def write_hole(self, dir_, index, x, y, r):
        """原盤中心からの位置 (x, y) に半径 r の穴を出力します。枠外は無視します。"""
        self._write_frame(index)
        if not (math.isfinite(x) and math.isfinite(y) and abs(x) <= self.r and abs(y) <= self.r):
            return
        out = self.outs[index // (self.column * self.row)]
        self._draw_hole(out, self._get_cx(index) + x, self._get_cy(dir_, index) + y, r)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class PlateWriterPDF(_PlateWriterBase):
    """ISO A4（210 x 297 mm）の1ページPDFをベクトルで書き出します。"""

    def __init__(self, *args, **kwargs):
        from reportlab.lib.units import mm
        from reportlab.pdfgen.canvas import Canvas

        super().__init__(*args, **kwargs)
        self._mm = mm
        self._canvas_class = Canvas

    def _x(self, value):
        """mm単位のx座標をPDFのpoint単位に変換します。"""
        return value * self._mm

    def _y(self, value):
        """y軸下向きのmm座標をPDFのy軸上向き座標に変換します。"""
        return (self.HEIGHT * 2 - value) * self._mm

    def _open_page(self, page):
        return self._canvas_class(
            self._output_path(page, "pdf"),
            pagesize=(self._x(self.WIDTH), self._x(self.HEIGHT * 2)),
            pageCompression=1,
            invariant=1,
        )

    def _draw_frame(self, out, cx, cy, label):
        out.setFillGray(0)
        out.setFont("Helvetica", self._x(3))
        out.drawCentredString(self._x(cx), self._y(cy - self.r), label)
        out.setLineWidth(self._x(0.1))
        out.setStrokeGray(0)
        out.setFillGray(0)
        x, y, side = cx - self.r, cy - self.r, self.r * 2
        out.rect(
            self._x(x), self._y(y + side), self._x(side), self._x(side),
            stroke=0 if self.invert_color else 1,
            fill=1 if self.invert_color else 0,
        )

    def _draw_hole(self, out, x, y, r):
        out.setFillGray(1 if self.invert_color else 0)
        out.circle(self._x(x), self._y(y), self._x(r), stroke=0, fill=1)

    def close(self):
        for index, out in enumerate(self.outs):
            destination = self._output_path(index, "pdf")
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
        self.outs = []


class PlateWriterSVG(_PlateWriterBase):
    def _open_page(self, page):
        out = open(self._output_path(page, "svg"), "w", encoding="utf-8")
        out.write("<?xml version=\"1.0\"?>")
        out.write("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
                  "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">")
        out.write(f"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{self.WIDTH}mm\" "
                  f"height=\"{self.HEIGHT * 2}mm\" version=\"1.1\">")
        return out

    def _draw_frame(self, out, cx, cy, label):
        out.write(f"<text x=\"{cx}mm\" y=\"{cy - self.r}mm\" text-anchor=\"middle\" "
                  f"font-size=\"3mm\" font-familiy=\"Verdana\">{label}</text>")
        out.write(f"<rect x=\"{cx - self.r}mm\" y=\"{cy - self.r}mm\" "
                  f"width=\"{self.r * 2}mm\" height=\"{self.r * 2}mm\" "
                  f"fill=\"{'black' if self.invert_color else 'none'}\" "
                  f"stroke=\"{'white' if self.invert_color else 'black'}\" "
                  f"stroke-width=\"0.1mm\" />")

    def _draw_hole(self, out, x, y, r):
        out.write(f"<circle cx=\"{x}mm\" cy=\"{y}mm\" r=\"{r}mm\" "
                  f"fill=\"{'white' if self.invert_color else 'black'}\"/>")

    def close(self):
        for out in self.outs:
            out.write("</svg>")
            out.close()
        self.outs = []

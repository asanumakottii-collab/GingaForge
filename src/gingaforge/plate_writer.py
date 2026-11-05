# -*- coding: utf-8 -*-
#
# Orb Transform Library(OTL)
# Copyright (C) 2007,2012,2014 東京大学地文研究会天文部
# Copyright (C) 2026 東京大学地文研究会天文部
# 2026年にGingaForge向けに移行・変更しました。
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
"""原板データをSVGまたは印刷用PDFとして書き出すモジュールです。穴はベクトル、開口率の画像は画素で書き出します。"""

import base64
from io import BytesIO
import math
import os
import tempfile
from xml.sax.saxutils import escape


class _PlateWriterBase:
    """PlateWriterPDF と PlateWriterSVG に共通するレイアウト計算。"""

    WIDTH = 210.  # ISO A4用紙の幅
    HEIGHT = 297. / 2.  # ISO A4用紙の高さの半分（上段・下段に配置）
    LABEL_SIZE = 3.  # 原盤名の文字の大きさ (mm)
    HEADING_SIZE = 2.5  # 列の見出しの文字の大きさの上限 (mm)
    HEADING_TOP = 5.  # 用紙の上端から列の見出しまでの余白 (mm)

    def __init__(self, column, row, r, filename_prefix, invert_color, output_dir, column_headings=((),),
                 plates=None, rows=1):
        """r: 円形の原盤枠の半径。0なら用紙サイズと配置数から自動計算する。

        column_headings: 列ごとの見出しの行。用紙を列の数で横に等分し、各列に上下の原盤を一式ずつ置く。
        plates: 枠を出力する原盤 ``(出力位置, 番号)`` の集合。Noneなら全原盤。
        他の原盤の位置は変えず、ここにない原盤の場所は空ける。原盤が1枚もない用紙は作らない。
        rows: 行の数。全列の原盤を行の数だけ、同じ内容で上から並べる。原盤のない段（N/S）の場所は行に取らない。
        """
        if column <= 0:
            raise ValueError("横の値が不正です。")
        if row <= 0:
            raise ValueError("縦の値が不正です。")
        if not column_headings:
            raise ValueError("列の数が不正です。")
        if rows <= 0:
            raise ValueError("行の数が不正です。")
        self.outs = {}  # 用紙番号から開いた用紙
        self.plates = plates
        self.column_headings = [tuple(lines) for lines in column_headings]
        self.number_of_plate = [0] * len(self.column_headings)  # 列ごとの出力済みの原板の数
        self.column = column
        self.row = row
        self.rows = rows
        self.sections = [0, 1] if plates is None else sorted({dir_ for dir_, _ in plates})  # 原盤のある段
        self.section_height = self.HEIGHT * 2 / max(2, rows * len(self.sections))
        side = min(self.WIDTH / (len(self.column_headings) * column), self.section_height / row)
        self.r = (side - 10.) / 2 if r == 0 else r
        if not 0 < self.r * 2 <= side:
            raise ValueError(f"原盤が用紙に収まりません。plate.frame-diameter を {side:g} mm 以下にするか、"
                             "原盤・列・行の数を減らしてください。")
        self.filename_prefix = filename_prefix
        self.invert_color = invert_color  # Falseなら原板を白で星を黒, Trueなら原板を黒で星を白
        self.output_dir = output_dir
        self.heading_size = self._fit_headings()

    def _fit_headings(self):
        """列の見出しが列の幅と、最上段の原盤名より上に収まる文字の大きさ (mm)。"""
        lines = max(map(len, self.column_headings))
        if not lines:
            return self.HEADING_SIZE
        from reportlab.pdfbase.pdfmetrics import stringWidth

        width = self.WIDTH / len(self.column_headings) - 2.
        widest = max(stringWidth(line, "Helvetica", 1.) for lines in self.column_headings for line in lines)
        height = self._get_cy(self.sections[0], 0) - self.r - self.LABEL_SIZE - self.HEADING_TOP
        size = min(self.HEADING_SIZE, width / widest if widest else self.HEADING_SIZE, height / (1.25 * lines))
        if size < 1.:
            raise ValueError("列の見出し（列で変えた設定）が原盤と重なります。plate.frame-diameter を小さくするか、"
                             "列で変える設定を減らしてください。")
        return size

    def _get_cx(self, index, column=0):
        """枠の中心のx座標"""
        width = self.WIDTH / len(self.column_headings)
        return width * column + width / (self.column * 2) * ((index % (self.column * self.row)) // self.row * 2 + 1)

    def _get_cy(self, dir_, index, row=0):
        """行 row の枠の中心のy座標。1行なら、原盤のない段があっても上半分がN、下半分がSのまま。"""
        section = dir_ if self.rows == 1 else row * len(self.sections) + self.sections.index(dir_)
        height = self.section_height
        return height / (self.row * 2) * ((index % self.row) * 2 + 1) + height * section

    def _output_path(self, page, extension):
        os.makedirs(self.output_dir, exist_ok=True)
        return os.path.join(self.output_dir, f"{self.filename_prefix}{page}.{extension}")

    def _page(self, index):
        """原盤 index の用紙を返します。初めてなら開いて、列の見出しを出力します。"""
        page = index // (self.column * self.row)
        if page not in self.outs:
            out = self.outs[page] = self._open_page(page)
            for i, lines in enumerate(self.column_headings):
                cx = self.WIDTH / len(self.column_headings) * (i + .5)
                for j, line in enumerate(lines):
                    self._draw_text(out, cx, self.HEADING_TOP + self.heading_size * (1 + 1.25 * j),
                                    line, self.heading_size)
        return self.outs[page]

    def _write_frame(self, index, column=0):
        """指定番号までの原盤枠を上下・全行とも出力し、必要な用紙を開きます。"""
        while index >= self.number_of_plate[column]:
            number = self.number_of_plate[column]
            cx = self._get_cx(number, column)
            for i in range(2):
                if self.plates is not None and (i, number) not in self.plates:
                    continue
                label = ("N" if i == 0 else "S") + str(number)
                if len(self.column_headings) > 1:
                    label += f" #{column + 1}"
                for row in range(self.rows):
                    self._draw_frame(self._page(number), cx, self._get_cy(i, number, row), label)
            self.number_of_plate[column] += 1

    def write_frames(self, number_of_units):
        """全列の上下それぞれに指定個数の原盤枠を、星がない場合も出力します。"""
        if number_of_units > 0:
            for column in range(len(self.column_headings)):
                self._write_frame(number_of_units - 1, column)

    def inside_frame(self, x, y, margin=0.):
        """原盤中心からの位置 (x, y) が、枠の縁から margin 以上内側にあるか。"""
        if not (math.isfinite(x) and math.isfinite(y)):
            return False
        return math.hypot(x, y) + margin <= self.r

    def write_hole(self, dir_, index, x, y, r, column=0):
        """全行の列 column の原盤中心からの位置 (x, y) に半径 r の穴を出力します。枠外は無視します。"""
        self._write_frame(index, column)
        if not self.inside_frame(x, y):
            return
        for row in range(self.rows):
            self._draw_hole(self._page(index), self._get_cx(index, column) + x, self._get_cy(dir_, index, row) + y, r)

    def write_image(self, dir_, index, pixels, left, top, width, height, column=0):
        """全行の列 column の原盤に、開口率の画像 pixels（uint8、255が光を通す）を枠の円で切り抜いて出力します。

        left・top は原盤中心から画像の左上まで、width・height は画像の大きさ (mm)。
        白地に黒穴（color.invert = no）では、光を通す所を黒にするため値を反転します。
        """
        from PIL import Image

        self._write_frame(index, column)
        out = self._page(index)
        image = Image.fromarray(pixels if self.invert_color else 255 - pixels, "L")
        cx = self._get_cx(index, column)
        for row in range(self.rows):
            cy = self._get_cy(dir_, index, row)
            self._draw_image(out, cx, cy, image, cx + left, cy + top, width, height)

    def __enter__(self):
        return self

    def __exit__(self, exc_type, exc_val, exc_tb):
        self.close()


class PlateWriterPDF(_PlateWriterBase):
    """ISO A4（210 x 297 mm）のPDFをベクトルで書き出します。"""

    PLATE_FORM = "plate"  # mirrored で、両ページが共有する描画内容の名前

    def __init__(self, *args, mirrored=False, **kwargs):
        """mirrored: Trueなら2ページ目に、裏面印刷用に文字も含めて左右反転した同じ内容を追加する。"""
        from reportlab.lib.units import mm
        from reportlab.pdfgen.canvas import Canvas

        super().__init__(*args, **kwargs)
        self._mm = mm
        self._canvas_class = Canvas
        self.mirrored = mirrored

    def _x(self, value):
        """mm単位のx座標をPDFのpoint単位に変換します。"""
        return value * self._mm

    def _y(self, value):
        """y軸下向きのmm座標をPDFのy軸上向き座標に変換します。"""
        return (self.HEIGHT * 2 - value) * self._mm

    def _open_page(self, page):
        out = self._canvas_class(
            self._output_path(page, "pdf"),
            pagesize=(self._x(self.WIDTH), self._x(self.HEIGHT * 2)),
            pageCompression=1,
            invariant=1,
        )
        if self.mirrored:
            out.beginForm(self.PLATE_FORM)
        return out

    def _finish_pages(self, out):
        """描画を終えた用紙のページを確定します。mirrored では同じ内容を2ページに置きます。"""
        if self.mirrored:
            out.endForm()
            out.doForm(self.PLATE_FORM)
            out.showPage()
            out.translate(self._x(self.WIDTH), 0)
            out.scale(-1, 1)
            out.doForm(self.PLATE_FORM)
        out.showPage()

    def _draw_text(self, out, cx, y, text, size):
        out.setFillGray(0)
        out.setFont("Helvetica", self._x(size))
        out.drawCentredString(self._x(cx), self._y(y), text)

    def _draw_frame(self, out, cx, cy, label):
        self._draw_text(out, cx, cy - self.r, label, self.LABEL_SIZE)
        out.setLineWidth(self._x(0.1))
        out.setStrokeGray(0)
        out.setFillGray(0)
        stroke, fill = (0, 1) if self.invert_color else (1, 0)
        out.circle(self._x(cx), self._y(cy), self._x(self.r), stroke=stroke, fill=fill)

    def _draw_hole(self, out, x, y, r):
        out.setFillGray(1 if self.invert_color else 0)
        out.circle(self._x(x), self._y(y), self._x(r), stroke=0, fill=1)

    def _draw_image(self, out, cx, cy, image, x, y, width, height):
        from reportlab.lib.utils import ImageReader

        out.saveState()
        frame = out.beginPath()
        frame.circle(self._x(cx), self._y(cy), self._x(self.r))
        out.clipPath(frame, stroke=0, fill=0)
        # DeviceGray・補間なしで埋め込むため、画素がそのままプリンターの点に対応する。
        out.drawImage(ImageReader(image), self._x(x), self._y(y + height), self._x(width), self._x(height))
        out.restoreState()

    def close(self):
        for index, out in self.outs.items():
            destination = self._output_path(index, "pdf")
            temporary = None
            try:
                # 既存ファイルへの直接書き込みは同期ストレージ等でタイムアウト
                # することがある。同じディレクトリに完成させてから原子的に置換。
                with tempfile.NamedTemporaryFile(
                        mode="wb", dir=os.path.dirname(os.path.abspath(destination)),
                        prefix=".otl-pdf-", suffix=".tmp", delete=False) as stream:
                    temporary = stream.name
                    self._finish_pages(out)
                    stream.write(out.getpdfdata())
                os.replace(temporary, destination)
            except OSError as exc:
                raise OSError(
                    exc.errno, f"PDFの保存に失敗しました: {exc.strerror}", destination
                ) from exc
            finally:
                if temporary is not None and os.path.exists(temporary):
                    os.unlink(temporary)
        self.outs = {}


class PlateWriterSVG(_PlateWriterBase):
    def _open_page(self, page):
        out = open(self._output_path(page, "svg"), "w", encoding="utf-8")
        out.write("<?xml version=\"1.0\"?>")
        out.write("<!DOCTYPE svg PUBLIC \"-//W3C//DTD SVG 1.1//EN\" "
                  "\"http://www.w3.org/Graphics/SVG/1.1/DTD/svg11.dtd\">")
        out.write(f"<svg xmlns=\"http://www.w3.org/2000/svg\" width=\"{self.WIDTH}mm\" "
                  f"height=\"{self.HEIGHT * 2}mm\" version=\"1.1\">")
        return out

    def _draw_text(self, out, cx, y, text, size):
        # 文字の幅はPDFと同じHelveticaの字幅で見積もる。
        out.write(f"<text x=\"{cx}mm\" y=\"{y}mm\" text-anchor=\"middle\" font-size=\"{size}mm\" "
                  f"font-family=\"Helvetica, Arial, sans-serif\">{escape(text)}</text>")

    def _draw_frame(self, out, cx, cy, label):
        self._draw_text(out, cx, cy - self.r, label, self.LABEL_SIZE)
        out.write(f"<circle cx=\"{cx}mm\" cy=\"{cy}mm\" r=\"{self.r}mm\" "
                  f"fill=\"{'black' if self.invert_color else 'none'}\" "
                  f"stroke=\"{'white' if self.invert_color else 'black'}\" "
                  f"stroke-width=\"0.1mm\" />")

    def _draw_hole(self, out, x, y, r):
        out.write(f"<circle cx=\"{x}mm\" cy=\"{y}mm\" r=\"{r}mm\" "
                  f"fill=\"{'white' if self.invert_color else 'black'}\"/>")

    def _draw_image(self, out, cx, cy, image, x, y, width, height):
        # 切り抜きの名前は枠の中心で決める（1つの枠に画像は1枚）。
        clip = f"frame-{cx:.4f}-{cy:.4f}"
        png = BytesIO()
        image.save(png, "PNG", optimize=True)
        data = base64.b64encode(png.getvalue()).decode("ascii")
        out.write(f"<clipPath id=\"{clip}\"><circle cx=\"{cx}mm\" cy=\"{cy}mm\" r=\"{self.r}mm\"/></clipPath>"
                  f"<image x=\"{x}mm\" y=\"{y}mm\" width=\"{width}mm\" height=\"{height}mm\" "
                  f"preserveAspectRatio=\"none\" image-rendering=\"optimizeSpeed\" clip-path=\"url(#{clip})\" "
                  f"xmlns:xlink=\"http://www.w3.org/1999/xlink\" xlink:href=\"data:image/png;base64,{data}\"/>")

    def close(self):
        for out in self.outs.values():
            out.write("</svg>")
            out.close()
        self.outs = {}

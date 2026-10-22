"""Mellinger全天画像の緑成分を中心配置の円筒投影へ再投影し、OHPフィルム用の原盤PDFを生成する。"""
from __future__ import annotations

import argparse
import configparser
import math
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from PIL import Image
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas
from .config import PROJECT_ROOT, output_directory
from .galaxy_geometry import REGIONS


def read_config(path):
    parser = configparser.ConfigParser(interpolation=None)
    parser.read_string("[ohp]\n" + Path(path).read_text(encoding="utf-8"))
    return parser["ohp"]


def geometry(cfg):
    p = cfg.getfloat("plate-distance-mm")
    bmax = cfg.getfloat("latitude-limit-deg")
    scale = cfg.getfloat("print-scale")
    dpi = cfg.getint("raster-dpi")
    if not (p > 0 and 0 < bmax < 90 and scale > 0 and 72 <= dpi <= 2400):
        raise ValueError("Invalid geometry or raster resolution")
    if cfg.getfloat("projector-horizontal-mm") != 0 or cfg.getfloat("projector-vertical-mm") != 0:
        raise ValueError("This OHP renderer supports only the centered projector")
    width = 2 * p * math.tan(math.radians(bmax)) * scale
    height = p * math.radians(60) * scale
    if max(width, height) > 130:
        raise ValueError("Image exceeds the available 130 mm area in the OTL layout")
    return p, bmax, scale, dpi, width, height


def plate_to_sky(x_mm, y_mm, center_deg, plate_distance_mm, scale=1.):
    return (center_deg - np.degrees(np.asarray(y_mm) / (plate_distance_mm * scale))) % 360., np.degrees(np.arctan(np.asarray(x_mm) / (plate_distance_mm * scale)))


def celestial_wcs(header):
    h = header.copy()
    # Normalize the legacy keyword only in memory; preserve the source FITS.
    if "RADECSYS" in h:
        h["RADESYS"] = h.pop("RADECSYS")
    return WCS(h, naxis=2)


def bilinear(data, x, y):
    ny, nx = data.shape
    if not (np.all(np.isfinite(x)) and np.all(np.isfinite(y))):
        raise ValueError("Nonfinite WCS coordinates")
    if np.any((y < 0) | (y > ny - 1)):
        raise ValueError("Requested sky lies outside the source image")
    x0, y0 = np.floor(x).astype(int), np.floor(y).astype(int)
    dx, dy = x - x0, y - y0
    xa, xb = x0 % nx, (x0 + 1) % nx
    ya, yb = y0, np.minimum(y0 + 1, ny - 1)
    return ((1 - dy) * ((1 - dx) * data[ya, xa] + dx * data[ya, xb])
            + dy * ((1 - dx) * data[yb, xa] + dx * data[yb, xb]))


def reproject(data, wcs, start, end, cfg):
    p, _, scale, dpi, width, height = geometry(cfg)
    nx, ny = round(width / 25.4 * dpi), round(height / 25.4 * dpi)
    x = (np.arange(nx) + .5) / nx * width - width / 2
    output = np.empty((ny, nx), dtype=np.float32)
    for row in range(0, ny, 128):
        y = (np.arange(row, min(row + 128, ny)) + .5) / ny * height - height / 2
        xx, yy = np.meshgrid(x, y)
        longitude, latitude = plate_to_sky(xx, yy, (start + end) / 2, p, scale)
        sx, sy = wcs.world_to_pixel_values(longitude, latitude)
        output[row:row + len(y)] = bilinear(data, sx, sy)
    if not np.all(np.isfinite(output)):
        raise ValueError("Source contains nonfinite values in the requested region")
    return output


def tone_map(values, cfg):
    black = cfg.getfloat("black-level")
    white = cfg.getfloat("white-level")
    soft = cfg.getfloat("asinh-softening")
    exposure = cfg.getfloat("exposure")
    gamma = cfg.getfloat("display-gamma")
    if not (white > black and soft > 0 and exposure > 0 and gamma > 0):
        raise ValueError("Invalid tone mapping settings")
    z = np.clip((values.astype(np.float64) - black) * exposure / (white - black), 0, 1)
    mode = cfg.get("tone-mode")
    if mode == "asinh":
        z = np.arcsinh(z / soft) / np.arcsinh(1 / soft)
    elif mode != "linear":
        raise ValueError("tone-mode must be asinh or linear")
    return np.rint(255 * z ** (1 / gamma)).astype(np.uint8)


def draw_plate(c, name, region, image, width, height, cy):
    cx, side = 105., 138.5
    c.setFillColorRGB(0, 0, 0)
    c.rect((cx - side / 2) * mm, (cy - side / 2) * mm, side * mm, side * mm, fill=1, stroke=0)
    c.drawImage(ImageReader(image), (cx - width / 2) * mm,
                (cy - height / 2) * mm, width * mm, height * mm, mask=None)
    c.setFillColorRGB(1, 1, 1)
    c.setFont("Helvetica-Bold", 14)
    c.drawString((cx - side / 2 + 5) * mm, (cy + side / 2 - 9) * mm, name)
    c.setFont("Helvetica", 8)
    c.drawString((cx - side / 2 + 5) * mm, (cy + side / 2 - 15) * mm,
                 f"l = {region[0]:.0f} to {region[1]:.0f} deg   |   Mellinger 2009, G")
    c.drawCentredString(cx * mm, (cy + height / 2 + 4) * mm, f"l = {region[1]:.0f} deg")
    c.drawCentredString(cx * mm, (cy - height / 2 - 6) * mm, f"l = {region[0]:.0f} deg")
    c.drawCentredString(cx * mm, (cy - side / 2 + 9) * mm, "LEFT: -b     RIGHT: +b     WHITE: CLEAR FILM")
    c.setStrokeColorRGB(1, 1, 1)
    c.setLineWidth(.25)
    # All registration marks are outside the projected field.
    for sx in [-1, 1]:
        for sy in [-1, 1]:
            px, py = cx + sx * (width / 2 + 2), cy + sy * (height / 2 + 2)
            c.line((px - 1) * mm, py * mm, (px + 1) * mm, py * mm)
            c.line(px * mm, (py - 1) * mm, px * mm, (py + 1) * mm)


def make_master(path, images, cfg):
    _, bmax, scale, _, width, height = geometry(cfg)
    c = canvas.Canvas(str(path), pagesize=A4, pageCompression=1)
    c.setTitle("Mellinger OHP masters - N0 N1 S0 S1")
    c.setAuthor("Source sky image: Axel Mellinger (2009)")
    for names in [("N0", "S0"), ("N1", "S1")]:
        for name, cy in zip(names, [222.75, 74.25]):
            draw_plate(c, name, REGIONS[name], images[name], width, height, cy)
        c.setFillColorRGB(0, 0, 0)
        c.setFont("Helvetica", 6)
        c.drawCentredString(105 * mm, 148 * mm,
            f"PRINT AT 100% / NO FIT-TO-PAGE | scale {scale:g} | field {width:.4f} x {height:.4f} mm | |b| <= {bmax:g} deg")
        c.showPage()
    c.save()


def build(config_path):
    """設定ファイルに従ってOHP原盤PDFを生成し、そのパスを返します。"""
    cfg = read_config(config_path)
    directory = Path(config_path).resolve().parent
    source = directory / cfg["source-fits"]
    output = output_directory(directory, "ohp")
    geometry(cfg)
    data = fits.getdata(source)
    if data.shape != (3, 1800, 3600):
        raise ValueError("Expected the Mellinger 3600 x 1800 RGB cube")
    wcs = celestial_wcs(fits.getheader(source))
    # The middle plane of the RGB cube is green; it is not Gaia's broad G band.
    images = {name: Image.fromarray(tone_map(reproject(data[1], wcs, start, end, cfg), cfg))
              for name, (start, end) in REGIONS.items()}
    output.mkdir(parents=True, exist_ok=True)
    path = output / "mellinger-ohp-plates.pdf"
    make_master(path, images, cfg)
    print("Created", path, flush=True)
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--config", type=Path, default=PROJECT_ROOT / "galaxyconfig.properties")
    build(parser.parse_args(argv).config)
    return 0

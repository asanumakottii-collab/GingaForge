#!/usr/bin/env python3
"""Reproject Mellinger's observed G image to the centered OTL cylinder."""
from __future__ import annotations

import argparse
import configparser
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from PIL import Image, ImageDraw, ImageFont
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas
from galaxy_geometry import REGIONS

ROOT = Path(__file__).resolve().parent


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


def draw_plate(c, name, region, image_path, width, height, cy):
    cx, side = 105., 138.5
    c.setFillColorRGB(0, 0, 0)
    c.rect((cx - side / 2) * mm, (cy - side / 2) * mm, side * mm, side * mm, fill=1, stroke=0)
    c.drawImage(ImageReader(str(image_path)), (cx - width / 2) * mm,
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


def make_master(output, cfg):
    _, bmax, scale, dpi, width, height = geometry(cfg)
    c = canvas.Canvas(str(output / "mellinger-ohp-plates.pdf"), pagesize=A4, pageCompression=1)
    c.setTitle("Mellinger OHP masters - N0 N1 S0 S1")
    c.setAuthor("Source sky image: Axel Mellinger (2009)")
    for names in [("N0", "S0"), ("N1", "S1")]:
        for name, cy in zip(names, [222.75, 74.25]):
            draw_plate(c, name, REGIONS[name], output / f"{name}-ohp.png", width, height, cy)
        c.setFillColorRGB(0, 0, 0)
        c.setFont("Helvetica", 6)
        c.drawCentredString(105 * mm, 148 * mm,
            f"PRINT AT 100% / NO FIT-TO-PAGE | scale {scale:g} | field {width:.4f} x {height:.4f} mm | |b| <= {bmax:g} deg")
        c.showPage()
    c.save()


def make_calibration(output, cfg):
    c = canvas.Canvas(str(output / "ohp-calibration.pdf"), pagesize=A4, pageCompression=1)
    c.setTitle("OHP transmission and print scale check")
    c.setFont("Helvetica-Bold", 17)
    c.drawString(20 * mm, 275 * mm, "OHP transmission / print scale check")
    c.setFont("Helvetica", 10)
    for y, text in zip([263, 256, 249, 242], [
        "Print with exactly the same film, printer, ink and driver settings as the masters.",
        "Print at 100%; disable fit-to-page. Gray code is not calibrated transmission.",
        "0 = maximum black ink; 255 = clear film. Inspect under the projector light.",
        "Film density, projector brightness and viewing adaptation require a trial."]):
        c.drawString(20 * mm, y * mm, text)
    for index, code in enumerate(np.rint(np.linspace(0, 255, 21)).astype(int)):
        column, row = index % 7, index // 7
        x, y = 20 + column * 24, 207 - row * 32
        c.setFillGray(code / 255)
        c.setStrokeGray(0)
        c.rect(x * mm, y * mm, 20 * mm, 20 * mm, fill=1, stroke=1)
        c.setFillGray(0)
        c.drawCentredString((x + 10) * mm, (y - 5) * mm, str(code))
    c.setFillGray(0)
    c.setFont("Helvetica", 10)
    c.drawString(20 * mm, 123 * mm, "50.00 mm horizontal and vertical references (measure after printing)")
    c.setLineWidth(.5)
    c.line(25 * mm, 106 * mm, 75 * mm, 106 * mm)
    c.line(95 * mm, 106 * mm, 95 * mm, 56 * mm)
    for x in [25, 75]:
        c.line(x * mm, 104 * mm, x * mm, 108 * mm)
    for y in [56, 106]:
        c.line(93 * mm, y * mm, 97 * mm, y * mm)
    c.drawString(20 * mm, 38 * mm, "Sky source: Axel Mellinger (2009), PASP 121, 1180-1187.")
    c.drawString(20 * mm, 31 * mm, "Source use conditions: https://milkywaysky.com/licenses.html")
    c.save()


def make_preview(output):
    preview = Image.new("RGB", (1100, 1620), "#171a20")
    draw = ImageDraw.Draw(preview)
    font = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 24)
    small = ImageFont.truetype("/System/Library/Fonts/Supplemental/Arial.ttf", 19)
    draw.text((30, 18), "Mellinger G / OHP masters - PREVIEW ONLY", font=font, fill="white")
    for name, (px, py) in zip(["N0", "N1", "S0", "S1"], [(35, 75), (580, 75), (35, 845), (580, 845)]):
        start, end = REGIONS[name]
        draw.text((px, py), f"{name}   l = {start:.0f}-{end:.0f} deg", font=font, fill="white")
        im = Image.open(output / f"{name}-ohp.png").convert("RGB")
        im.thumbnail((440, 635), Image.Resampling.LANCZOS)
        preview.paste(im, (px + 20, py + 42))
        draw.text((px, py + 690), "Left: b=-20 / Right: b=+20 deg", font=small, fill="#ccd3db")
        draw.text((px, py + 715), "Longitude decreases downward", font=small, fill="#ccd3db")
    preview.save(output / "preview.png")


def build(config_path):
    cfg = read_config(config_path)
    directory = Path(config_path).resolve().parent
    source = directory / cfg["source-fits"]
    output = directory / cfg["output-directory"]
    output.mkdir(parents=True, exist_ok=True)
    _, _, _, dpi, width, height = geometry(cfg)
    header = fits.getheader(source)
    if fits.getdata(source).shape != (3, 1800, 3600):
        raise ValueError("Expected the Mellinger 3600 x 1800 RGB cube")
    # The middle plane of the RGB cube is green; it is not Gaia's broad G band.
    data = fits.getdata(source)[1]
    wcs = celestial_wcs(header)
    stats = {}
    for name, (start, end) in REGIONS.items():
        linear = reproject(data, wcs, start, end, cfg)
        np.save(output / f"{name}-linear-G.npy", linear)
        pixels = tone_map(linear, cfg)
        Image.fromarray(pixels).save(output / f"{name}-ohp.png", dpi=(dpi, dpi))
        stats[name] = {"longitude_deg": [start, end], "latitude_deg": [-cfg.getfloat("latitude-limit-deg"), cfg.getfloat("latitude-limit-deg")],
            "shape_rows_columns": list(linear.shape), "display_zero_fraction": float(np.mean(pixels == 0)),
            "display_white_fraction": float(np.mean(pixels == 255)),
            "linear_percentiles_0_50_90_99_100": np.percentile(linear, [0, 50, 90, 99, 100]).tolist()}
        print(name, stats[name], flush=True)
    make_master(output, cfg)
    make_calibration(output, cfg)
    make_preview(output)
    provenance = {"source_url": "https://galaxy.phy.cmich.edu/~axel/mwpan2/mwpan2_RGB_3600.fits",
        "source_author": "Axel Mellinger", "source_paper_doi": "10.1086/648480",
        "license_information": "https://milkywaysky.com/licenses.html",
        "license_clearance": "Not obtained by this workflow; check with the author for educational screenings.",
        "source_sha256": hashlib.file_digest(source.open("rb"), "sha256").hexdigest(),
        "source_plane_zero_based": 1, "source_sampling_deg": .1,
        "source_intensity_units": "Original FITS values; no BUNIT in source header",
        "projection": "Centered OTL cylinder: x=p*tan(b), y=-p*(l-center); radians; y down on page",
        "sampling": "Bilinear interpolation of intensity at pixel centers; no Jacobian intensity factor",
        "printed_field_width_mm": width, "printed_field_height_mm": height,
        "config": dict(cfg), "regions": stats}
    (output / "provenance.json").write_text(json.dumps(provenance, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    (output / "source-fits-header.txt").write_text(header.tostring(sep="\n", endcard=True, padding=False) + "\n")
    print("Created", output / "mellinger-ohp-plates.pdf", flush=True)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("-f", "--config", type=Path, default=ROOT / "galaxyconfig.properties")
    build(parser.parse_args(argv).config)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

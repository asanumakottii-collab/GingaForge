# Copyright (C) 2026 東京大学地文研究会天文部; GPL-2.0-or-later
"""Mellinger全天画像を中心配置の円筒投影へ再投影し、OHPフィルム用の原盤PDFを生成する。

FITSの青・緑・赤を重み付きで混ぜ、StarNet2で星と天の川の光に分ける。恒星投影機・1等星投影機が
映すHipparcosの星（V等級が上限以下）は除き、それより暗い星の光はぼかして天の川の光に戻す。
"""
from __future__ import annotations

import argparse
import configparser
import gzip
import math
import subprocess
import tempfile
from pathlib import Path

import numpy as np
from astropy.io import fits
from astropy.wcs import WCS
from astropy.wcs.utils import proj_plane_pixel_scales
from PIL import Image
from reportlab.lib.pagesizes import A4
from reportlab.lib.units import mm
from reportlab.lib.utils import ImageReader
from reportlab.pdfgen import canvas
from . import starnet_setup
from .config import PROJECT_ROOT, output_directory
from .galaxy_geometry import REGIONS, lens_distance
from .models import SpherePosition

# StarNet2へ渡す16ビット画像の可逆なasinh変換。上限を超えるのは明るい星の芯だけ。
STARNET_RANGE = 20000.
STARNET_SOFTENING = 50.
# A4に2列×2段で4区画を並べる。枠の一辺は区画の幅（210 / 2 mm）から余白10 mmを引いたもの。
PLATE_LAYOUT = {"N0": (52.5, 222.75), "S0": (157.5, 222.75), "N1": (52.5, 74.25), "S1": (157.5, 74.25)}
FRAME_SIDE = 95.


def read_config(path):
    parser = configparser.ConfigParser(interpolation=None)
    parser.read_string("[ohp]\n" + Path(path).read_text(encoding="utf-8"))
    return parser["ohp"]


def geometry(cfg):
    """``(投影, 銀緯の上限, 倍率, dpi, 画像の幅, 画像の高さ)`` を返す。

    投影は ``(原板とレンズの距離, ドーム中心からレンズ中心までの距離, ドーム半径)``（倍率を掛ける前のmm）。
    OTLと同じく、原板は投影機中心から plate-distance-mm の位置にあり、原板とレンズの距離はレンズの公式から決まる。
    """
    plate = cfg.getfloat("plate-distance-mm")
    focus = cfg.getfloat("focal-length-mm")
    sphere = cfg.getfloat("dome-radius-mm")
    bmax = cfg.getfloat("latitude-limit-deg")
    scale = cfg.getfloat("print-scale")
    dpi = cfg.getint("raster-dpi")
    if not (plate >= 0 and focus > 0 and sphere > plate and 0 < bmax < 90 and scale > 0 and 72 <= dpi <= 2400):
        raise ValueError("Invalid geometry or raster resolution")
    if cfg.getfloat("projector-horizontal-mm") != 0 or cfg.getfloat("projector-vertical-mm") != 0:
        raise ValueError("This OHP renderer supports only the centered projector")
    lens = lens_distance(sphere - plate, focus)
    offset = plate + lens
    # 原盤の端は、レンズ中心から見たドーム上の銀経±30度・銀緯±bmaxの点。
    b = math.radians(bmax)
    width = 2 * lens * sphere * math.sin(b) / (sphere * math.cos(b) - offset) * scale
    height = 2 * lens * math.atan2(sphere * math.sin(math.radians(30)), sphere * math.cos(math.radians(30)) - offset) * scale
    if width > FRAME_SIDE - 16 or height > FRAME_SIDE - 52:
        raise ValueError(f"Image and labels do not fit the {FRAME_SIDE:g} mm frame of the 2 x 2 layout")
    return (lens, offset, sphere), bmax, scale, dpi, width, height


def plate_to_sky(x_mm, y_mm, center_deg, projection, scale=1.):
    """原盤上の点から、レンズ中心を通る光線がドームに当たる点の ``(銀経, 銀緯)``（度）を返す。

    projection は geometry() が返す ``(原板とレンズの距離, ドーム中心からレンズ中心までの距離, ドーム半径)``。
    """
    lens, offset, sphere = projection
    angle = -np.asarray(y_mm) / (lens * scale)
    # 原盤中心の方向・銀経の増える方向・銀北極の方向の成分。
    forward, east, north = np.cos(angle), np.sin(angle), np.asarray(x_mm) / (lens * scale)
    norm = np.sqrt(1 + north * north)
    forward, east, north = forward / norm, east / norm, north / norm
    along = offset * forward
    distance = -along + np.sqrt(along * along + sphere * sphere - offset * offset)
    forward = offset + distance * forward
    east, north = distance * east, distance * north
    return ((center_deg + np.degrees(np.arctan2(east, forward))) % 360.,
            np.degrees(np.arctan2(north, np.hypot(forward, east))))


def celestial_wcs(header):
    h = header.copy()
    # Normalize the legacy keyword only in memory; preserve the source FITS.
    if "RADECSYS" in h:
        h["RADESYS"] = h.pop("RADECSYS")
    return WCS(h, naxis=2)


def channel_weights(cfg):
    weights = np.array([float(w) for w in cfg.get("channel-weights").split(",")])
    if weights.shape != (3,) or not np.all(np.isfinite(weights)) or np.any(weights < 0) or weights.sum() <= 0:
        raise ValueError("channel-weights must be three non-negative numbers (blue, green, red)")
    return weights / weights.sum()


def channel_mix(cube, cfg):
    """FITSの3面を重み付き平均する。面の順は青・緑・赤（アンタレスとシャウラの色で確認）。"""
    return np.tensordot(channel_weights(cfg), cube.astype(np.float64), axes=1)


def gaussian_blur(image, sigma):
    """σ（画素）のガウスぼかし。全天画像なので経度方向は360度でつなぎ、緯度方向は端で折り返す。"""
    if sigma == 0:
        return image.copy()
    r = math.ceil(3 * sigma)
    kernel = np.exp(-.5 * (np.arange(-r, r + 1) / sigma) ** 2)
    kernel /= kernel.sum()
    padded = np.pad(image, ((r, r), (0, 0)), mode="reflect")
    image = sum(k * padded[i:i + image.shape[0]] for i, k in enumerate(kernel))
    padded = np.pad(image, ((0, 0), (r, r)), mode="wrap")
    return sum(k * padded[:, i:i + image.shape[1]] for i, k in enumerate(kernel))


def starnet(image, cfg, accept_license=False):
    """StarNet2で星を除いた画像を返す。可逆なasinh変換で16ビットTIFFにして渡し、出力を元の強度に戻す。"""
    command = starnet_setup.prepare(cfg.get("starnet-command"), accept=accept_license)
    print(f"Separating stars with StarNet2 {starnet_setup.version(command)}: {command}", flush=True)
    scale = np.arcsinh(STARNET_RANGE / STARNET_SOFTENING)
    stretched = np.arcsinh(np.clip(image, 0, STARNET_RANGE) / STARNET_SOFTENING) / scale
    with tempfile.TemporaryDirectory() as directory:
        source, starless = Path(directory) / "input.tif", Path(directory) / "starless.tif"
        Image.fromarray(np.rint(stretched * 65535).astype(np.uint16)).save(source)
        upsample = ["-u"] if cfg.getboolean("starnet-upsample") else []
        subprocess.run([command, "-q", *upsample, "-i", str(source), "-o", str(starless)], check=True)
        result = np.asarray(Image.open(starless), dtype=np.float64) / 65535
    if result.shape != image.shape:
        raise ValueError("StarNet2 returned an image of a different size")
    return STARNET_SOFTENING * np.sinh(result * scale)


def read_hipparcos(path, limit):
    """hip_main.dat（.gzも可）から、V等級が limit 以下の星の (赤経, 赤緯, V) を度で返す。

    OTLの恒星原盤と同じ欄（V: 42–46桁、ICRS・J1991.25の赤経・赤緯: 52–63桁と65–76桁）を読む。
    OTLと同じく、等級か位置が空欄の星は投影されないので含めない。
    """
    opener = gzip.open if Path(path).suffix == ".gz" else open
    stars = []
    with opener(path, "rt", encoding="ascii") as f:
        for line in f:
            vmag, ra, dec = line[41:46], line[51:63], line[64:76]
            if vmag.strip() and ra.strip() and dec.strip() and float(vmag) <= limit:
                stars.append((float(ra), float(dec), float(vmag)))
    return np.array(stars).reshape(-1, 3)


def icrs_to_galactic(ra, dec):
    position = SpherePosition()
    position.radeg, position.dedeg = ra, dec
    return position.to_galactic()


def star_disks(shape, x, y, radius):
    """画素座標 (x, y) を中心とする半径 radius（画素）の円の和集合。経度方向は360度でつなぐ。"""
    mask = np.zeros(shape, bool)
    for cx, cy, r in zip(x, y, radius):
        rows = np.arange(max(0, math.floor(cy - r)), min(shape[0], math.ceil(cy + r) + 1))
        cols = np.arange(math.floor(cx - r), math.ceil(cx + r) + 1)
        inside = (rows[:, None] - cy) ** 2 + (cols[None, :] - cx) ** 2 <= r * r
        mask[np.ix_(rows, cols % shape[1])] |= inside
    return mask


def remove_projected_stars(image, wcs, cfg, config_directory, accept_license=False):
    """恒星投影機・1等星投影機が映す星を除き、それより暗い星の光をぼかして天の川の光に戻す。"""
    limit = cfg.getfloat("star-catalog-limit-mag")
    radius0 = cfg.getfloat("star-mask-radius-deg")
    growth = cfg.getfloat("star-mask-growth-deg-per-mag")
    halo_limit = cfg.getfloat("halo-fill-limit-mag")
    blur = cfg.getfloat("faint-star-blur-deg")
    cap = cfg.getfloat("faint-star-cap")
    weight = cfg.getfloat("faint-star-weight")
    if not (radius0 > 0 and growth >= 0 and blur >= 0 and cap > 0 and 0 <= weight <= 1):
        raise ValueError("Invalid star removal settings")
    catalog = read_hipparcos(config_directory / cfg["star-catalog"], limit)
    longitude, latitude = np.array([icrs_to_galactic(ra, dec) for ra, dec, _ in catalog]).reshape(-1, 2).T
    x, y = wcs.world_to_pixel_values(longitude, latitude)
    pixel = float(np.mean(proj_plane_pixel_scales(wcs)))
    vmag = catalog[:, 2]
    radius = (radius0 + growth * (limit - vmag)) / pixel
    # Star cores above STARNET_RANGE include single-pixel spikes; clipping keeps them out of the faint-star glow.
    clipped = np.clip(image, 0, STARNET_RANGE)
    starless = starnet(clipped, cfg, accept_license)
    keep = (~star_disks(image.shape, x, y, radius)).astype(np.float64)
    # The source was downsampled bilinearly, so a sharp star lands in one pixel at its peak brightness and
    # its pixel value varies tens of times between stars of the same magnitude; the cap keeps that out of the glow.
    faint = np.minimum(clipped - starless, cap) * keep
    # Inside the removed discs the glow takes the local average of the fainter stars around them.
    sigma = blur / pixel
    glow = gaussian_blur(faint, sigma) / np.maximum(gaussian_blur(keep, sigma), 1e-12) if sigma else faint
    # StarNet2 leaves a faint halo around the brightest stars; refill it from the surroundings with a soft edge.
    bright = vmag <= halo_limit
    halo = star_disks(image.shape, x[bright], y[bright], 1.5 * radius[bright]).astype(np.float64)
    if halo.any():
        outside = 1. - halo
        spread = max(1., 1.5 * float(radius[bright].max()) / 2)
        filled = gaussian_blur(starless * outside, spread) / np.maximum(gaussian_blur(outside, spread), 1e-12)
        starless = starless + np.clip(gaussian_blur(halo, 1.), 0, 1) * (filled - starless)
    print(f"Removed {len(vmag)} Hipparcos stars with V <= {limit:g} (all sky); fainter stars blurred by {blur:g} deg", flush=True)
    return starless + weight * glow


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
    projection, _, scale, dpi, width, height = geometry(cfg)
    nx, ny = round(width / 25.4 * dpi), round(height / 25.4 * dpi)
    x = (np.arange(nx) + .5) / nx * width - width / 2
    output = np.empty((ny, nx), dtype=np.float32)
    for row in range(0, ny, 128):
        y = (np.arange(row, min(row + 128, ny)) + .5) / ny * height - height / 2
        xx, yy = np.meshgrid(x, y)
        longitude, latitude = plate_to_sky(xx, yy, (start + end) / 2, projection, scale)
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


def draw_plate(c, name, region, image, width, height, cx, cy, caption):
    side = FRAME_SIDE
    c.setFillColorRGB(0, 0, 0)
    c.rect((cx - side / 2) * mm, (cy - side / 2) * mm, side * mm, side * mm, fill=1, stroke=0)
    c.drawImage(ImageReader(image), (cx - width / 2) * mm,
                (cy - height / 2) * mm, width * mm, height * mm, mask=None)
    c.setFillColorRGB(1, 1, 1)
    c.setFont("Helvetica-Bold", 12)
    c.drawString((cx - side / 2 + 5) * mm, (cy + side / 2 - 8) * mm, name)
    c.setFont("Helvetica", 7)
    c.drawString((cx - side / 2 + 5) * mm, (cy + side / 2 - 13) * mm,
                 f"l = {region[0]:.0f} to {region[1]:.0f} deg")
    c.setFont("Helvetica", 6)
    c.drawString((cx - side / 2 + 5) * mm, (cy + side / 2 - 17) * mm, caption)
    c.setFont("Helvetica", 8)
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


def caption(cfg):
    weights = ":".join(f"{w:g}" for w in channel_weights(cfg))
    text = f"Mellinger 2009, B:G:R = {weights}"
    if cfg.get("star-removal") == "starnet":
        text += f", StarNet2 + HIP V <= {cfg.getfloat('star-catalog-limit-mag'):g} removed"
    return text


def make_master(path, images, cfg):
    (lens, _, _), bmax, scale, _, width, height = geometry(cfg)
    c = canvas.Canvas(str(path), pagesize=A4, pageCompression=1)
    c.setTitle("Mellinger OHP masters - N0 N1 S0 S1")
    c.setAuthor("Source sky image: Axel Mellinger (2009)")
    for name, (cx, cy) in PLATE_LAYOUT.items():
        draw_plate(c, name, REGIONS[name], images[name], width, height, cx, cy, caption(cfg))
    c.setFillColorRGB(0, 0, 0)
    c.setFont("Helvetica", 6)
    c.drawCentredString(105 * mm, 148 * mm,
        f"PRINT AT 100% / NO FIT-TO-PAGE | scale {scale:g} | focal length {cfg.getfloat('focal-length-mm') * scale:g} mm"
        f" | lens-plate {lens * scale:.4f} mm | field {width:.4f} x {height:.4f} mm | |b| <= {bmax:g} deg")
    c.showPage()
    c.save()


def milky_way_image(config_path, cfg, accept_starnet_license=False):
    """再投影前の全天画像（色の混合・星の除去済み）とそのWCSを返す。"""
    source = config_path.parent / cfg["source-fits"]
    data = fits.getdata(source)
    if data.shape != (3, 1800, 3600):
        raise ValueError("Expected the Mellinger 3600 x 1800 RGB cube")
    wcs = celestial_wcs(fits.getheader(source))
    image = channel_mix(data, cfg)
    removal = cfg.get("star-removal")
    if removal == "starnet":
        image = remove_projected_stars(image, wcs, cfg, config_path.parent, accept_starnet_license)
    elif removal != "none":
        raise ValueError("star-removal must be starnet or none")
    return image, wcs


def build(config_path, accept_starnet_license=False):
    """設定ファイルに従ってOHP原盤PDFを生成し、そのパスを返します。"""
    config_path = Path(config_path).resolve()
    cfg = read_config(config_path)
    output = output_directory(config_path)
    geometry(cfg)
    image, wcs = milky_way_image(config_path, cfg, accept_starnet_license)
    images = {name: Image.fromarray(tone_map(reproject(image, wcs, start, end, cfg), cfg))
              for name, (start, end) in REGIONS.items()}
    output.mkdir(parents=True, exist_ok=True)
    path = output / f"{config_path.stem}-plates.pdf"
    make_master(path, images, cfg)
    print("Created", path, flush=True)
    return path


def main(argv=None):
    parser = argparse.ArgumentParser(prog="gingaforge mellinger-ohp", description=__doc__)
    parser.add_argument("-f", "--config", type=Path, default=PROJECT_ROOT / "mellinger-ohp.properties")
    parser.add_argument("--accept-starnet-license", action="store_true",
                        help="StarNet2 の利用条件（同梱の LICENSE.txt）に同意する（対話できない環境用）")
    args = parser.parse_args(argv)
    build(args.config, args.accept_starnet_license)
    return 0

"""Gaia DR3のG帯光量マップを読む。ネットワーク取得は明示指定時だけ行う。

ESA TAPで銀経・銀緯ビンごとに集計した光量をキャッシュし、検証してから読み込む。
集計は恒星数ではなく SUM(10**(-0.4*(G-reference)))。
参考: https://www.cosmos.esa.int/web/gaia-users/archive/programmatic-access
Copyright (C) 2026 東京大学地文研究会天文部; GPL-2.0-or-later
"""

import csv
import gzip
import hashlib
import json
import math
import os
from pathlib import Path
import ssl
import tempfile
import time
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from dataclasses import dataclass

from .config import positive_number
from .models import SpherePosition

TAP_URL = "https://gea.esac.esa.int/tap-server/tap"


@dataclass(frozen=True)
class GaiaSettings:
    path: Path
    bright: float = 7.5
    faint: float = 17.0
    reference: float = 7.5
    sky_step: float = 0.1

    @classmethod
    def from_properties(cls, props, config_directory):
        path = Path(props.get("galaxy.gaia.input", "source/gaia"))
        if not path.is_absolute():
            path = Path(config_directory) / path
        bright = positive_number(props, "galaxy.gaia.bright-limit", 7.5, allow_zero=True)
        faint = positive_number(props, "galaxy.gaia.faint-limit", 17.)
        reference = positive_number(props, "galaxy.gaia.reference-magnitude", 7.5, allow_zero=True)
        step = positive_number(props, "galaxy.gaia.sky-bin-deg", 0.1)
        if not bright < faint <= 21 or reference > 21:
            raise ValueError("Gaia等級は 0 <= bright-limit < faint-limit <= 21 としてください。")
        if not 0.05 <= step <= 1. or not math.isclose(round(60 / step) * step, 60):
            raise ValueError("galaxy.gaia.sky-bin-deg: 60度を等分する0.05〜1度を指定してください。")
        return cls(path, bright, faint, reference, step)

    def provenance(self, regions, latitude_limit):
        return {
            "format_version": 1, "catalogue": "gaiadr3.gaia_source", "band": "G",
            "bright_exclusive": self.bright, "faint_inclusive": self.faint,
            "reference_magnitude": self.reference, "sky_bin_deg": self.sky_step,
            "longitude_regions": [list(r) for r in regions], "latitude_limit": latitude_limit,
        }


def build_queries(settings, regions, latitude_limit):
    """全該当星の光量を天球上で集計。TOP・ランダムな間引きは使わない。"""
    step = settings.sky_step
    lon = f"FLOOR(l / {step:.12g})"
    lat = f"FLOOR((b + 90.0) / {step:.12g})"
    return [
        f"SELECT {lon} AS lon_bin, {lat} AS lat_bin, COUNT(*) AS source_count,\n"
        f"SUM(POWER(10.0, -0.4 * (phot_g_mean_mag - {settings.reference:.12g}))) AS relative_flux\n"
        "FROM gaiadr3.gaia_source\n"
        f"WHERE l >= {start:.12g} AND l < {end:.12g}\n"
        f"AND b >= {-latitude_limit:.12g} AND b <= {latitude_limit:.12g}\n"
        f"AND phot_g_mean_mag > {settings.bright:.12g} AND phot_g_mean_mag <= {settings.faint:.12g}\n"
        "GROUP BY lon_bin, lat_bin"
        for start, end in regions
    ]


def file_sha256(path):
    digest = hashlib.sha256()
    with open(path, "rb") as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def iter_flux_samples(settings, regions, latitude_limit):
    """(ビン中心のICRS位置, 基準等級に対する光量)をストリームで返す。"""
    for lon, lat, flux in _iter_cached_bins(settings, regions, latitude_limit):
        yield SpherePosition.from_galactic(lon, lat), flux


def iter_flux_bins(settings, regions, latitude_limit):
    """(銀経の下端, 上端, 銀緯の下端, 上端, 基準等級に対する光量)を度でストリームで返す。

    キャッシュはビンの中心だけを持つので、集計幅から範囲を戻す。銀緯±latitude_limitで
    切り詰めたビンは、中心が切り詰め後の中点であることから下端・上端を求める。
    """
    half = settings.sky_step / 2
    for lon, lat, flux in _iter_cached_bins(settings, regions, latitude_limit):
        bottom, top = lat - half, lat + half
        if top > latitude_limit:
            bottom, top = 2 * lat - latitude_limit, latitude_limit
        elif bottom < -latitude_limit:
            bottom, top = -latitude_limit, 2 * lat + latitude_limit
        yield lon - half, lon + half, bottom, top, flux


def _iter_cached_bins(settings, regions, latitude_limit):
    """(ビン中心の銀経, 銀緯, 光量)を返す。

    キャッシュはmanifestの条件とチェックサムを照合してから読む。
    """
    path = settings.path
    manifest_path = path / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(
            f"Gaiaキャッシュがありません: {path}\n"
            "同じ設定で --download-gaia を付けて実行してください。"
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    expected = settings.provenance(regions, latitude_limit)
    if manifest.get("selection") != expected or not manifest.get("complete"):
        raise ValueError("Gaiaキャッシュの等級・集計幅・担当領域が設定と一致しません。--download-gaia で再取得してください。")
    files = manifest.get("files", [])
    if [entry["name"] for entry in files] != [f"region-{i}.csv.gz" for i in range(len(regions))]:
        raise ValueError("Gaiaキャッシュの領域ファイル一覧が不正です。")
    for entry in files:
        if file_sha256(path / entry["name"]) != entry["sha256"]:
            raise ValueError(f"Gaiaキャッシュのチェックサムが一致しません: {entry['name']}")
    for entry in files:
        with gzip.open(path / entry["name"], "rt", encoding="utf-8", newline="") as stream:
            reader = csv.DictReader(stream)
            if reader.fieldnames != ["l", "b", "relative_flux", "source_count"]:
                raise ValueError("Gaia集計CSVの列が不正です。")
            for row in reader:
                lon, lat, flux = (float(row[k]) for k in ("l", "b", "relative_flux"))
                if not all(math.isfinite(x) for x in (lon, lat, flux)) or flux <= 0:
                    raise ValueError("Gaia集計CSVに不正な光量があります。")
                if not (0 <= lon < 360 and -latitude_limit <= lat <= latitude_limit):
                    raise ValueError("Gaia集計CSVの座標が領域外です。")
                yield lon, lat, flux


def _request(url, data=None):
    if data is not None:
        data = urllib.parse.urlencode(data).encode("ascii")
    context = ssl.create_default_context()
    # python.org版macOS Pythonで既定CAファイルが未配置の場合も、検証を無効化しない。
    if ssl.get_default_verify_paths().cafile is None and Path("/etc/ssl/cert.pem").is_file():
        context.load_verify_locations("/etc/ssl/cert.pem")
    return urllib.request.urlopen(urllib.request.Request(url, data=data, headers={
        "User-Agent": "Via-Lactea-Generator/1.0",
    }), timeout=60, context=context)


def _download_votable(query, target, maxrec):
    # UWS非同期ジョブ。同期TAPの短い実行時間制限を避ける。
    with _request(TAP_URL + "/async", {
        "REQUEST": "doQuery", "LANG": "ADQL", "FORMAT": "votable_plain",
        "QUERY": query, "MAXREC": str(maxrec), "PHASE": "RUN",
    }) as response:
        job = response.geturl().rstrip("/")
    if "/async/" not in job or not job.startswith(TAP_URL + "/async/"):
        raise RuntimeError("ESA TAPがジョブURLを返しませんでした。")
    print(f"Gaiaジョブ: {job}", flush=True)
    deadline = time.monotonic() + 1800
    last_phase = None
    while time.monotonic() < deadline:
        with _request(job + "/phase") as response:
            phase = response.read().decode().strip()
        if phase != last_phase:
            print(f"  {phase}", flush=True)
            last_phase = phase
        if phase == "COMPLETED":
            break
        if phase in {"ERROR", "ABORTED"}:
            with _request(job + "/error") as response:
                detail = response.read(4000).decode(errors="replace")
            raise RuntimeError(f"Gaia TAP {phase}: {detail}")
        time.sleep(5)
    else:
        raise TimeoutError(f"Gaiaジョブが30分で完了しませんでした: {job}")
    with _request(job + "/results/result") as response, open(target, "wb") as out:
        for chunk in iter(lambda: response.read(1024 * 1024), b""):
            out.write(chunk)
    return job


def _votable_to_csv(source, target, settings, region, latitude_limit):
    """TAPのOVERFLOW/ERRORを検出。途中までのデータを完成品として使わない。"""
    fields, statuses = [], []
    rows = stars = 0
    seen = set()
    with gzip.open(target, "wt", encoding="utf-8", newline="") as out:
        writer = csv.writer(out)
        writer.writerow(["l", "b", "relative_flux", "source_count"])
        for event, element in ET.iterparse(source, events=("end",)):
            tag = element.tag.rsplit("}", 1)[-1]
            if tag == "FIELD":
                fields.append(element.attrib.get("name", "").lower())
            elif tag == "INFO" and element.attrib.get("name") == "QUERY_STATUS":
                status = element.attrib.get("value", "").upper()
                statuses.append(status)
                if status != "OK":
                    raise RuntimeError(f"Gaia TAP {status}: {element.text or ''}")
            elif tag in {"BINARY", "BINARY2", "FITS"}:
                raise RuntimeError("TABLEDATA形式のVOTableが必要です。")
            elif tag == "TR":
                record = dict(zip(fields, [child.text for child in element]))
                lon_value, lat_value = float(record["lon_bin"]), float(record["lat_bin"])
                if not lon_value.is_integer() or not lat_value.is_integer():
                    raise ValueError("Gaia集計ビン番号が整数ではありません。")
                lon_bin, lat_bin = int(lon_value), int(lat_value)
                if (lon_bin, lat_bin) in seen:
                    raise ValueError("Gaia集計セルが重複しています。")
                seen.add((lon_bin, lat_bin))
                lon = (lon_bin + 0.5) * settings.sky_step
                bottom = max(-latitude_limit, lat_bin * settings.sky_step - 90.)
                top = min(latitude_limit, (lat_bin + 1) * settings.sky_step - 90.)
                if bottom > top:
                    raise ValueError("Gaia集計ビンが銀緯範囲外です。")
                lat = (bottom + top) / 2
                flux, count = float(record["relative_flux"]), int(record["source_count"])
                if not (region[0] <= lon < region[1] and -latitude_limit <= lat <= latitude_limit
                        and math.isfinite(flux) and flux > 0 and count > 0):
                    raise ValueError("Gaia TAPから不正な集計セルが返されました。")
                writer.writerow([f"{lon:.10f}", f"{lat:.10f}", repr(flux), count])
                rows += 1
                stars += count
                element.clear()
    if not statuses or any(s != "OK" for s in statuses):
        raise RuntimeError("Gaia TAPの完了ステータスを確認できません。")
    if set(fields) != {"lon_bin", "lat_bin", "relative_flux", "source_count"}:
        raise RuntimeError("Gaia TAPの列が想定と異なります。")
    return rows, stars


def download_gaia(settings, regions, latitude_limit):
    """4領域を別々に取得。条件一致・検証済みの領域は再利用する。"""
    directory = settings.path
    directory.mkdir(parents=True, exist_ok=True)
    queries = build_queries(settings, regions, latitude_limit)
    entries = []
    for index, (query, region) in enumerate(zip(queries, regions)):
        name = f"region-{index}.csv.gz"
        path = directory / name
        meta_path = directory / f"region-{index}.json"
        if path.exists() and meta_path.exists():
            meta = json.loads(meta_path.read_text(encoding="utf-8"))
            if meta.get("query") == query and meta.get("sha256") == file_sha256(path):
                print(f"Gaia領域{index}: 検証済みキャッシュを使用", flush=True)
                entries.append(meta)
                continue
        print(f"Gaia領域{index}: 銀経{region[0]:g}〜{region[1]:g}度を集計・取得", flush=True)
        with tempfile.TemporaryDirectory(dir=directory) as temp:
            xml_path, csv_path = Path(temp) / "result.xml", Path(temp) / name
            maxrec = math.ceil((region[1] - region[0]) / settings.sky_step) * (math.ceil(2 * latitude_limit / settings.sky_step) + 2) + 10
            job = _download_votable(query, xml_path, maxrec)
            rows, stars = _votable_to_csv(xml_path, csv_path, settings, region, latitude_limit)
            meta = {"name": name, "query": query, "job": job, "rows": rows, "stars": stars,
                    "sha256": file_sha256(csv_path)}
            os.replace(csv_path, path)
            meta_path.write_text(json.dumps(meta, indent=2) + "\n", encoding="utf-8")
        entries.append(meta)
        print(f"  {stars:,}天体 → {rows:,}光量セル", flush=True)
    manifest = {"complete": True, "selection": settings.provenance(regions, latitude_limit), "files": entries}
    partial = directory / "manifest.json.part"
    partial.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    os.replace(partial, directory / "manifest.json")

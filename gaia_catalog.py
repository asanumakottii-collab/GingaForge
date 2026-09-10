"""Gaia DR3のG帯光量を読む。ネットワーク取得は明示指定時だけ行う。

個別星CSV、またはESA TAPで集計した光量マップを入力する。
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

from models import SpherePosition

TAP_URL = "https://gea.esac.esa.int/tap-server/tap"


def positive_number(props, key, default, *, allow_zero=False):
    try:
        value = float(props.get(key, default))
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{key}: 数値を指定してください。") from exc
    if not math.isfinite(value) or (value < 0 if allow_zero else value <= 0):
        raise ValueError(f"{key}: 有限の{'0以上' if allow_zero else '正'}の数値が必要です。")
    return value


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


def _read_csv(path):
    opener = gzip.open if str(path).endswith(".gz") else open
    return opener(path, "rt", encoding="utf-8-sig", newline="")


def iter_flux_samples(settings, regions, latitude_limit, stats):
    """(ICRS位置, 基準等級に対する光量, 恒星数)をストリームで返す。

    生CSVは source_id,ra,dec,phot_g_mean_mag が必要。重複IDは拒否する。
    キャッシュはmanifestの条件とチェックサムを照合してから読む。
    """
    path = settings.path
    stats.update(rows=0, stars=0, filtered_stars=0)
    if path.is_dir():
        manifest = json.loads((path / "manifest.json").read_text(encoding="utf-8"))
        expected = settings.provenance(regions, latitude_limit)
        if manifest.get("selection") != expected or not manifest.get("complete"):
            raise ValueError("Gaiaキャッシュの等級・集計幅・担当領域が設定と一致しません。--download-gaia で再取得してください。")
        files = manifest.get("files", [])
        if [entry["name"] for entry in files] != [f"region-{i}.csv.gz" for i in range(len(regions))]:
            raise ValueError("Gaiaキャッシュの領域ファイル一覧が不正です。")
        for entry in files:
            if file_sha256(path / entry["name"]) != entry["sha256"]:
                raise ValueError(f"Gaiaキャッシュのチェックサムが一致しません: {entry['name']}")
        stats["provenance"] = manifest
        for entry in files:
            with _read_csv(path / entry["name"]) as stream:
                reader = csv.DictReader(stream)
                if reader.fieldnames != ["l", "b", "relative_flux", "source_count"]:
                    raise ValueError("Gaia集計CSVの列が不正です。")
                for row in reader:
                    lon, lat, flux = (float(row[k]) for k in ("l", "b", "relative_flux"))
                    count = int(row["source_count"])
                    if not all(math.isfinite(x) for x in (lon, lat, flux)) or flux <= 0 or count <= 0:
                        raise ValueError("Gaia集計CSVに不正な光量・個数があります。")
                    if not (0 <= lon < 360 and -latitude_limit <= lat <= latitude_limit):
                        raise ValueError("Gaia集計CSVの座標が領域外です。")
                    stats["rows"] += 1
                    stats["stars"] += count
                    yield SpherePosition.from_galactic(lon, lat), flux, count
        return
    if not path.is_file():
        raise FileNotFoundError(
            f"Gaiaデータがありません: {path}\n"
            "同じ設定で --download-gaia を付けて実行するか、galaxy.gaia.input にGaia CSVを指定してください。"
        )
    stats["provenance"] = {"file": str(path), "sha256": file_sha256(path), "band": "G",
                           "selection": settings.provenance(regions, latitude_limit),
                           "coverage": "user-supplied; completeness not verified"}
    ids = set()
    with _read_csv(path) as stream:
        reader = csv.DictReader(stream)
        if not {"source_id", "ra", "dec", "phot_g_mean_mag"}.issubset(reader.fieldnames or []):
            raise ValueError("Gaia CSVには source_id,ra,dec,phot_g_mean_mag 列が必要です。")
        for line_number, row in enumerate(reader, 2):
            try:
                source_id = int(row["source_id"])
                ra, dec, mag = (float(row[k]) for k in ("ra", "dec", "phot_g_mean_mag"))
                if source_id <= 0 or not all(math.isfinite(v) for v in (ra, dec, mag)):
                    raise ValueError("非有限値または不正なID")
                if not 0 <= ra < 360 or not -90 <= dec <= 90:
                    raise ValueError("座標範囲外")
                if source_id in ids:
                    raise ValueError(f"source_id重複: {source_id}")
                ids.add(source_id)
            except (ValueError, TypeError, KeyError) as exc:
                raise ValueError(f"{path}:{line_number}: Gaiaレコード不正 ({exc})") from exc
            stats["rows"] += 1
            if not settings.bright < mag <= settings.faint:
                stats["filtered_stars"] += 1
                continue
            position = SpherePosition()
            position.radeg, position.dedeg = ra, dec
            stats["stars"] += 1
            yield position, 10 ** (-0.4 * (mag - settings.reference)), 1


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
    if settings.path.suffix in {".csv", ".gz"}:
        raise ValueError("ダウンロード時のgalaxy.gaia.inputにはキャッシュディレクトリを指定してください。")
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

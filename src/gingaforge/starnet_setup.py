# Copyright (C) 2026 東京大学地文研究会天文部; GPL-2.0-or-later
"""StarNet2（Mikita (Nikita) Misiura）のコマンドライン版を探し、必要なら公式配布元から導入する。

StarNet2 は再配布できないため、リポジトリには含めない。公式の最新版フィードから今のOS向けの
ZIPを取得し、サイズとSHA-256を照合して tools/starnet2/（gitの対象外）に展開する。
利用条件（LICENSE.txt）は全文の場所を示してユーザーの同意を得てから使い、同意した条文のSHA-256を記録する。
参考: https://starnetastro.com/cli-tools/latest-versions/
"""

import argparse
import datetime
import json
import os
import platform
import re
import shutil
import stat
import subprocess
import sys
import tempfile
import zipfile
from pathlib import Path

from .config import PROJECT_ROOT
from .gaia_catalog import _request, file_sha256

FEED_URL = "https://starnetastro.com/cli-tools/latest.json"
FEED_SCHEMA = "starnetastro.cli-tools.latest.v1"
# 利用条件により、ダウンロードは公式サイトに掲載された配布元から直接行う。
OFFICIAL_PREFIXES = ("https://starnetastro.com/", "https://download.starnetastro.com/")
INSTALL_ROOT = PROJECT_ROOT / "tools" / "starnet2"
EXECUTABLES = ("starnet2", "starnet2.exe")
INSTALL_HINT = "gingaforge install-starnet で公式配布元から導入するか、設定の starnet-command にパスを書いてください。"


def platform_key():
    """公式フィードのプラットフォーム名（例: macos-arm64）。"""
    system, machine = platform.system(), platform.machine().lower()
    os_name = {"Darwin": "macos", "Linux": "linux", "Windows": "windows"}.get(system)
    arch = {"arm64": "arm64", "aarch64": "arm64", "x86_64": "x64", "amd64": "x64"}.get(machine)
    if os_name is None or arch is None or (os_name != "macos" and arch != "x64"):
        raise RuntimeError(f"StarNet2 のコマンドライン版がないプラットフォームです: {system} {machine}")
    return f"{os_name}-{arch}"


def latest_archive():
    """公式フィードから、今のOS向けの最新ZIPの情報（url, sha256, size_bytes, file_name, version, build）を返す。"""
    with _request(FEED_URL) as response:
        feed = json.load(response)
    if feed.get("schema") != FEED_SCHEMA:
        raise ValueError(f"StarNet2 の最新版フィードの形式が想定と違います: {feed.get('schema')}")
    entry = feed["tools"]["starnet2"]["latest"].get(platform_key()) or {}
    archive = entry.get("packages", {}).get("archive")
    if archive is None:
        raise ValueError(f"最新版フィードに {platform_key()} 向けのZIPがありません。"
                         "https://starnetastro.com/cli-tools/starnet/ から入手してください。")
    if not archive["url"].startswith(OFFICIAL_PREFIXES):
        raise ValueError(f"公式配布元ではないURLです: {archive['url']}")
    return {**archive, "version": entry["version"], "build": entry["build"]}


def download(archive, target):
    with _request(archive["url"]) as response, open(target, "wb") as out:
        shutil.copyfileobj(response, out, 1024 * 1024)
    if target.stat().st_size != archive["size_bytes"] or file_sha256(target) != archive["sha256"]:
        target.unlink()
        raise ValueError(f"ダウンロードしたZIPのサイズかSHA-256が公式フィードと一致しません: {archive['file_name']}")


def extract(zip_path, directory):
    """ZIPを展開する。展開先の外を指すパスは拒否し、実行権限とシンボリックリンクを保つ。"""
    directory = Path(directory).resolve()
    with zipfile.ZipFile(zip_path) as archive:
        for info in archive.infolist():
            target = (directory / info.filename).resolve()
            mode = info.external_attr >> 16
            link = archive.read(info).decode() if stat.S_ISLNK(mode) else None
            if not target.is_relative_to(directory) or (
                    link is not None and not (target.parent / link).resolve().is_relative_to(directory)):
                raise ValueError(f"ZIPに展開先の外を指すパスがあります: {info.filename}")
            if info.is_dir():
                target.mkdir(parents=True, exist_ok=True)
                continue
            target.parent.mkdir(parents=True, exist_ok=True)
            if link is not None:
                os.symlink(link, target)
                continue
            with archive.open(info) as source, open(target, "wb") as out:
                shutil.copyfileobj(source, out)
            if mode & 0o111:
                target.chmod(0o755)


def find_executable(directory):
    found = [p for name in EXECUTABLES for p in Path(directory).rglob(name) if p.is_file()]
    return min(found, key=lambda p: len(p.parts), default=None)


def version(executable):
    result = subprocess.run([str(executable), "--version"], capture_output=True, text=True)
    match = re.search(r"version:\s*(\S+)", result.stdout + result.stderr)
    return match.group(1) if match else "unknown"


def installed_executable(root=None):
    """tools/starnet2/ に導入済みのStarNet2の実行ファイル。なければ None。"""
    pointer = Path(root or INSTALL_ROOT) / "current.json"
    if pointer.is_file():
        executable = pointer.parent / json.loads(pointer.read_text(encoding="utf-8"))["executable"]
        if executable.is_file():
            return executable
    return None


def locate(command="", root=None):
    """使う実行ファイル。command があればそれ、空なら tools/starnet2/ の導入分、次に PATH 上の starnet2。"""
    if command:
        found = shutil.which(str(Path(command).expanduser()))
        if found is None:
            raise FileNotFoundError(f"starnet-command の StarNet2 が見つかりません: {command}")
        return Path(found)
    found = installed_executable(root) or shutil.which("starnet2")
    return Path(found) if found else None


def license_file(executable):
    """配布物の LICENSE.txt（ZIP版は実行ファイルの隣、インストーラー版は share/doc/starnet2/）。"""
    executable = Path(executable).resolve()
    for candidate in (executable.parent / "LICENSE.txt",
                      executable.parent.parent / "share" / "doc" / "starnet2" / "LICENSE.txt"):
        if candidate.is_file():
            return candidate
    return None


def _accepted(root):
    record = Path(root or INSTALL_ROOT) / "accepted-licenses.json"
    return json.loads(record.read_text(encoding="utf-8")) if record.is_file() else {}


def license_accepted(executable, root=None):
    path = license_file(executable)
    return path is not None and file_sha256(path) in _accepted(root)


def ready(command="", root=None):
    """StarNet2 が見つかり、その利用条件に同意済みか（尋ねずに確かめる）。"""
    try:
        executable = locate(command, root)
    except FileNotFoundError:
        return False
    return executable is not None and license_accepted(executable, root)


def ensure_license(executable, *, accept=False, interactive=None, root=None):
    """利用条件に同意済みか確かめる。未同意なら全文の場所を示し、同意を得て記録する（条文ごとに一度）。"""
    path = license_file(executable)
    if path is None:
        raise FileNotFoundError(f"StarNet2 の利用条件（LICENSE.txt）が実行ファイルの近くにありません: {executable}\n"
                                "公式サイトの配布物を、同梱ファイルごとそのままの構成で使ってください。")
    accepted = _accepted(root)
    digest = file_sha256(path)
    if digest in accepted:
        return
    print(f"星の除去には StarNet2（Mikita (Nikita) Misiura）を別のプロセスとして使います。\n"
          f"利用条件の全文: {path}", flush=True)
    if not accept:
        if not (sys.stdin.isatty() if interactive is None else interactive):
            raise PermissionError("StarNet2 の利用条件への同意が必要です。全文を読み、同意する場合は "
                                  "--accept-starnet-license を付けて実行してください。")
        print(path.read_text(encoding="utf-8", errors="replace"), flush=True)
        if input("StarNet2 の利用条件に同意しますか？ [y/N]: ").strip().lower() not in ("y", "yes"):
            raise PermissionError("StarNet2 の利用条件に同意しなかったため、StarNet2 を使いません。")
    accepted[digest] = {"license": str(path),
                        "accepted_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds")}
    record = Path(root or INSTALL_ROOT) / "accepted-licenses.json"
    record.parent.mkdir(parents=True, exist_ok=True)
    record.write_text(json.dumps(accepted, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def _check_local_zip(path):
    """手で取得したZIPが最新版と同名なら、公式フィードのSHA-256と照合する。"""
    try:
        archive = latest_archive()
    except (OSError, ValueError, RuntimeError) as exc:
        print(f"公式フィードを読めず、ZIPのSHA-256を照合できませんでした: {exc}", flush=True)
        return
    if archive["file_name"] != path.name:
        print(f"最新版（{archive['file_name']}）ではないため、SHA-256を照合していません: {path.name}", flush=True)
    elif file_sha256(path) != archive["sha256"]:
        raise ValueError(f"ZIPのSHA-256が公式フィードと一致しません: {path}")


def install(source=None, *, accept=False, interactive=None, root=None):
    """StarNet2 を tools/starnet2/ に導入し、実行ファイルのパスを返す。

    source を省略すると、公式フィードの最新版ZIPをダウンロードする。公式サイトから手で取得したZIPか、
    それを展開したフォルダを source に渡すこともできる。利用条件に同意しなければ導入しない。
    """
    root = Path(root or INSTALL_ROOT)
    root.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(dir=root) as work:
        work = Path(work)
        if source is None:
            archive = latest_archive()
            print(f"StarNet2 {archive['version']}（{archive['file_name']}、{archive['size_bytes'] / 1e6:.1f} MB）を"
                  f"ダウンロードします: {archive['url']}", flush=True)
            source = work / archive["file_name"]
            download(archive, source)
        else:
            source = Path(source).expanduser().resolve()
            if source.suffix.lower() == ".zip":
                _check_local_zip(source)
        if source.is_dir():
            unpacked = work / "package"
            shutil.copytree(source, unpacked, symlinks=True)
        else:
            unpacked = work / "unzipped"
            extract(source, unpacked)
        executable = find_executable(unpacked)
        if executable is None:
            raise FileNotFoundError(f"StarNet2 の実行ファイルが見つかりません: {source}")
        destination = root / source.name.removesuffix(".zip")
        if destination.exists():
            shutil.rmtree(destination)
        shutil.move(executable.parent, destination)
    executable = destination / executable.name
    try:
        ensure_license(executable, accept=accept, interactive=interactive, root=root)
    except PermissionError:
        shutil.rmtree(destination)
        raise
    (root / "current.json").write_text(json.dumps({
        "executable": str(executable.relative_to(root)), "version": version(executable), "source": str(source),
    }, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
    return executable


def prepare(command="", *, accept=False, interactive=None, root=None):
    """使うStarNet2を決め、利用条件への同意を確かめて実行ファイルのパスを返す。

    見つからない場合、対話できる端末なら公式配布元からの導入を提案する。
    """
    interactive = sys.stdin.isatty() if interactive is None else interactive
    executable = locate(command, root)
    if executable is None:
        if not interactive or input(
                "StarNet2 が見つかりません。公式サイト（starnetastro.com）から最新版をダウンロードして"
                "GingaForge の tools/starnet2/ に導入しますか？ [y/N]: ").strip().lower() not in ("y", "yes"):
            raise FileNotFoundError("StarNet2 が見つかりません。" + INSTALL_HINT)
        return install(accept=accept, interactive=interactive, root=root)
    ensure_license(executable, accept=accept, interactive=interactive, root=root)
    return executable


def main(argv=None):
    parser = argparse.ArgumentParser(
        prog="gingaforge install-starnet",
        description="StarNet2 のコマンドライン版を公式配布元から GingaForge の tools/starnet2/ に導入する。")
    parser.add_argument("--from", dest="source", type=Path,
                        help="公式サイトから手で取得したZIP、またはそれを展開したフォルダ（省略時はダウンロード）")
    parser.add_argument("--accept-starnet-license", action="store_true",
                        help="StarNet2 の利用条件（同梱の LICENSE.txt）に同意する")
    args = parser.parse_args(argv)
    executable = install(args.source, accept=args.accept_starnet_license)
    print(f"StarNet2 {version(executable)} を導入しました: {executable}", flush=True)
    return 0

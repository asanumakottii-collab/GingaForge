# Copyright (C) 2026 東京大学地文研究会天文部; GPL-2.0-or-later
"""GingaForge: 銀河投影機（ぎんとう）用の天の川原盤を生成するライブラリ。

- ``gingaforge.ohp.build(config)``: Mellinger全天画像からOHPフィルム用の原盤PDF
- ``gingaforge.gaia_etch.generate(config, pdf=False)``: Gaia DR3のG帯光量からエッチング用の穴原盤（SVG / PDF）
"""

__version__ = "0.1.0"

from .galaxy_geometry import REGIONS, GalaxyTransformer
from .models import PlatePosition, SpherePosition

__all__ = ["REGIONS", "GalaxyTransformer", "PlatePosition", "SpherePosition", "__version__"]

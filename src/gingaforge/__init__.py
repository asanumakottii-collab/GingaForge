# Copyright (C) 2026 東京大学地文研究会天文部; GPL-2.0-or-later
"""GingaForge: 銀河投影機（ぎんとう）用の天の川原盤を生成するライブラリ。

- ``gingaforge.ohp``: Mellinger全天画像からOHPフィルム用の階調原盤
- ``gingaforge.gaia_etch``: Gaia DR3のG帯光量からエッチング用の穴原盤
"""

__version__ = "0.1.0"

from .galaxy_geometry import REGIONS, GalaxyTransformer
from .models import PlatePosition, PlateStar, SpherePosition

__all__ = ["REGIONS", "GalaxyTransformer", "PlatePosition", "PlateStar", "SpherePosition", "__version__"]

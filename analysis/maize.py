"""玉米模型与产区数据的入口：实际代码在 ../maize-crf-map/（maize_model.py、maize_zones.py），
β 扫描（beta_scan.py）与玉米配方地图共用同一套模型。"""
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "maize-crf-map"))

from maize_model import P, season_setup  # noqa: E402,F401
from maize_zones import COUNTRIES, ENV_ZONE, ZONES  # noqa: E402,F401

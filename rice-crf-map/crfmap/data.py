"""读取 data/ 下的稻作区数据库（UTF-8 CSV，可直接用 Excel 编辑）。

zones.csv    每个稻作区一行：代表站、经纬度、一年几季（single/double/triple）、1–12 月平均气温（°C）
seasons.csv  每个稻季一行：
  est   TPR 移栽 / WSR 水直播（湿直播）/ DSR 旱直播（播后旱长、晚淹水）
  type  IND 常规籼稻 / HYB 杂交籼稻 / JAP 温带粳稻 / TJ 热带粳稻 / AROM 香稻、感光品种
  app   施肥日（移栽日或播种日），MM-DD
  days  施肥日到成熟的天数
  Y     可达目标产量（t/ha 稻谷）
  INS   土壤基础供氮（kg N/ha，整季，≈ 无氮区地上部吸氮量）
  FN    当地常规施氮量（kg N/ha，约数，仅作对照）
countries.csv  ISO 3166 数字代码 → 中文国名（地图着色用）

月均温是各代表站气候常年值的近似值；种植日历参考 RiceAtlas（Laborte et al. 2017）
与各国农业部门公布的常规播栽期。均为模型初值，应以实测数据替换。
"""
import csv
from pathlib import Path

DATA_DIR = Path(__file__).resolve().parent.parent / "data"


def _rows(name):
    with open(DATA_DIR / name, newline="", encoding="utf-8-sig") as f:
        return list(csv.DictReader(f))


def load_countries():
    return {r["iso_num"]: r["name_zh"] for r in _rows("countries.csv")}


def load_zones():
    zones, by_id = [], {}
    for r in _rows("zones.csv"):
        z = dict(id=r["id"], country=r["country"], name=r["name"], station=r["station"],
                 lat=float(r["lat"]), lon=float(r["lon"]),
                 T=[float(r[f"T{m:02d}"]) for m in range(1, 13)], crop=r["crop"], seasons=[])
        zones.append(z)
        by_id[z["id"]] = z
    for r in _rows("seasons.csv"):
        by_id[r["zone_id"]]["seasons"].append(dict(
            name=r["name"], est=r["est"], type=r["type"], app=r["app"], days=int(r["days"]),
            Y=float(r["Y"]), INS=int(r["INS"]), FN=int(r["FN"]), note=r["note"]))
    return zones


ZONES = load_zones()
COUNTRIES = load_countries()

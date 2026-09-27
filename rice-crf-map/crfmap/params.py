"""读取模型参数（默认为包目录上一级的 params.toml）。

对外提供与原 model.py 相同的名字：P、CROP、EST、PRODUCTS、MIN_FERT_DEMAND，
以及 PHEN（生育期推算）、OPT（优化）、SCEN（情景与对照）。
P 是可变字典：S 型产品情景会临时修改 P["beta"]，release() 在调用时读取它。
"""
import os
import tomllib
from pathlib import Path

PARAMS_FILE = Path(os.environ.get("CRFMAP_PARAMS", Path(__file__).resolve().parent.parent / "params.toml"))

with open(PARAMS_FILE, "rb") as _f:
    RAW = tomllib.load(_f)

_r, _p, _n = RAW["release"], RAW["pool"], RAW["product_n"]
# 键的顺序与 data.json 中的 params.P 保持一致
P = dict(
    Ea=_r["Ea"], f_soil=_r["f_soil"], f0=_r["f0"], beta=_r["beta"], paddy_dT=_r["paddy_dT"],
    urea_tau=_p["urea_tau"], urea_loss=_p["urea_loss"], k_loss=_p["k_loss"], eta=_p["eta"],
    buffer=_p["buffer"], early_pool=_p["early_pool"],
    n_cru=_n["n_cru"], n_urea=_n["n_urea"],
)
PRODUCTS = list(range(RAW["products"]["min"], RAW["products"]["max"] + 1, RAW["products"]["step"]))
CROP = {k: dict(nreq=v["nreq"], pPI=v["pPI"], pH=v["pH"], label=v["label"]) for k, v in RAW["crop"].items()}
EST = {k: dict(label=v["label"], dPI=v["dPI"], dH=v["dH"], k_mult=v["k_mult"]) for k, v in RAW["est"].items()}
MIN_FERT_DEMAND = RAW["demand"]["min_fert_demand"]
PHEN = RAW["phenology"]
OPT = RAW["optimize"]
SCEN = RAW["scenarios"]

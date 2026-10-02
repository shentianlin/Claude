"""生成“田间释放模拟器”页面：把模型参数与两张地图的站点气温嵌入 template.html。

    python field-release/build.py   →  field-release/index.html
"""
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
ROOT = HERE.parent
sys.path.insert(0, str(ROOT / "common"))
sys.path.insert(0, str(ROOT / "rice-crf-map"))
sys.path.insert(0, str(ROOT / "wheat-crf-map"))

import field_release as FR                     # noqa: E402
from crfmap.data import COUNTRIES as RC, ZONES as RZ   # noqa: E402
from zones import COUNTRIES as WC, ZONES as WZ          # noqa: E402

RICE_SYS = {"TPR": "rice_flooded", "WSR": "rice_flooded", "DSR": "rice_dsr"}
WHEAT_SYS = lambda z: "wheat_irrigated"


def sites():
    out, seen = [], set()
    for z in RZ:
        s = z["seasons"][0]
        out.append(dict(id="r-" + z["id"], country=RC[z["country"]], name=z["name"], station=z["station"],
                        crop="水稻", season=s["name"], T=z["T"], app=s["app"], days=s["days"],
                        system=RICE_SYS.get(s["est"], "rice_flooded")))
        seen.add(z["station"])
    for z in WZ:
        s = z["seasons"][0]
        m = int(s["app"][:2]); d = int(s["app"][3:])
        hm = int(s["harv"][:2]); hd = int(s["harv"][3:])
        days = ((hm - m) * 30.44 + (hd - d)) % 365
        sysk = "wheat_rainfed" if any(w in s["name"] for w in ("旱", "雨养")) else "wheat_irrigated"
        out.append(dict(id="w-" + z["id"], country=WC[z["country"]], name=z["name"], station=z["station"],
                        crop="小麦", season=s["name"], T=z["T"], app=s["app"], days=int(round(days)),
                        system=sysk))
    return out


def main():
    model = dict(DEFAULTS=FR.DEFAULTS, MOISTURE=FR.MOISTURE, SYSTEMS=FR.SYSTEMS, PLACEMENT=FR.PLACEMENT,
                 TEXTURE=FR.TEXTURE)
    html = (HERE / "template.html").read_text()
    html = html.replace("/*__MODEL__*/null", json.dumps(model, ensure_ascii=False))
    html = html.replace("/*__SITES__*/null", json.dumps(sites(), ensure_ascii=False))
    (HERE / "index.html").write_text(html)
    print("written", HERE / "index.html", len(html))


if __name__ == "__main__":
    main()

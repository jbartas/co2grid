# plot_points_lonx_laty.py
import json
from pathlib import Path
import matplotlib.pyplot as plt

# --- change this to your file path ---
PATH = Path(r"C:\jbartas\gridding\elements.json")

def load_points(path: Path):
    """Return (lats, lons) from a JSON array or NDJSON file."""
    lats, lons = [], []
    colors = []
    try:
        data = json.loads(path.read_text(encoding="utf-8"))
        if isinstance(data, dict):
            for key in ("items", "data", "elements", "points"):
                if key in data and isinstance(data[key], list):
                    data = data[key]
                    break
        if not isinstance(data, list):
            raise ValueError("JSON is not a list of objects")
        rows = data
    except Exception:
        rows = []
        with path.open("r", encoding="utf-8") as f:
            for line in f:
                line = line.strip()
                if line:
                    rows.append(json.loads(line))

    for rec in rows:
        if not isinstance(rec, dict): 
            continue
        lat = rec.get("lat", rec.get("latitude"))
        lon = rec.get("lon", rec.get("lng", rec.get("longitude")))
        co2 = rec.get("co2")
        try:
            lats.append(float(lat))
            lons.append(float(lon))
            colors.append(float(co2))
        except (TypeError, ValueError):
            continue
    return lats, lons, colors

def main():
    lats, lons, colors = load_points(PATH)
    print(f"Loaded {len(lats)} points from {PATH}")

    # X = longitude, Y = latitude  ✅
    plt.figure(figsize=(7, 6))
#    plt.scatter(lons, lats, s=6, alpha=0.6)
    plt.scatter(lons, lats, s=6, alpha=0.7, c=colors)
    plt.colorbar(label="CO₂")

    plt.xlabel("Longitude (deg)")   # X
    plt.ylabel("Latitude (deg)")    # Y
    plt.gca().set_aspect("equal")   # keep degrees square
    # Optional, if your data spans the globe:
    # plt.xlim(-180, 180); plt.ylim(-90, 90)
    plt.title("Points from JSON (X = lon, Y = lat)")
    plt.tight_layout()
    plt.savefig("points_scatter_lonx_laty.png", dpi=150)
    plt.show()

if __name__ == "__main__":
    main()

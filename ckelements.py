
import numpy as np
from math import sqrt

import json
from degridding import element, elementfile 

json_data = []   # global list of Element instances
elements = []  # local list

# file to check, from degridding.py
#filestr = "C:\\jbartas\\gridding\\elements.json"

with open(elementfile, "r", encoding="utf-8") as f:
    json_data = json.load(f)  # elements is now a list of {"lat":..., "lon":..., "co2":...} dicts

max_lat = -90
min_lat = 90
max_lon = 0
min_lon = 180

for d in json_data:
    # Cast to float in case the JSON stores strings
    obj = element(float(d["lat"]), float(d["lon"]), float(d["co2"]), float(d["u"]), float(d["v"]) )
    if(min_lat > obj.lat):
        min_lat = obj.lat
    if(max_lat < obj.lat):
        max_lat = obj.lat
    if(min_lon > obj.lon):
        min_lon = obj.lon
    if(max_lon < obj.lon):
        max_lon = obj.lon

    elements.append(obj)

print(len(json_data), "loaded")
print("latitude  max {} min {}".format(max_lat, min_lat))
print("longitude max {} min {}".format(max_lon, min_lon))





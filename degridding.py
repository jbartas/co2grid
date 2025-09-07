import numpy as np
from math import radians, sin, cos, sqrt, asin

def lldist(lat1, lon1, lat2, lon2):
    """
    Spherical coordinate distance between two points on Earth (km).
    """
    R = 6371.0  # Earth radius in km
    phi1, phi2 = radians(lat1), radians(lat2)
    d_phi = phi2 - phi1
    d_lam = radians(lon2 - lon1)
    a = sin(d_phi/2)**2 + cos(phi1)*cos(phi2)*sin(d_lam/2)**2
    return 2 * R * asin(min(1, sqrt(a)))


def latlon_to_km_vector(lat1, lon1, lat2, lon2):
    """
    Approximate vector (dy, dx) in km between two lat/lon points.
    """
    R = 6371.0
    phi1, phi2 = radians(lat1), radians(lat2)
    lam1, lam2 = radians(lon1), radians(lon2)
    dphi = phi2 - phi1
    dlam = lam2 - lam1
    meanphi = 0.5 * (phi1 + phi2)
    dy = R * dphi
    dx = R * cos(meanphi) * dlam
    return dy, dx

def grid_segment(lats, lons, co2s, us, vs,
                 ov,
                 #i_start, ov,
                 rows=8,
                 min_dist_samples=16,
                 threshold_multiplier=1.5,
                 invalid_max=2):
    """
    Reshape a list into rows x Ncol grid.
    Returns:
      lat_grid, lon_grid, co2_grid, u_grid, v_grid        shape (rows, Ncol)
      path_vecs, col_vecs                                 shape (Ncol, 2)
    """

    n_pts = len(lats)

    # compute local spacing threshold from the first gaps
    M = min(min_dist_samples, n_pts-1)
    # find the distance between points
    d0 = [lldist(lats[i], lons[i], lats[i+1], lons[i+1]) for i in range(M)]
    # threshold is some multiplier times the minimum distance
    thresh = threshold_multiplier * min(d0)

    # check where the first point is (not necessarrilly at start of a row)
    # look at first "row" of points
    #i_row = i_start

    # allocate maximum‐possible grid, fill with NaN
    max_cols = int(np.ceil(n_pts/rows)) # Greatest number of columns
    lat_grid = np.full((rows, max_cols), np.nan)
    lon_grid = np.full((rows, max_cols), np.nan)
    co2_grid = np.full((rows, max_cols), np.nan)
    u_grid = np.full((rows, max_cols), np.nan)
    v_grid = np.full((rows, max_cols), np.nan)
    path_vecs  = np.zeros((rows, max_cols, 2))  # (dy, dx) for all points that have data
    col_vecs   = np.zeros((rows, max_cols, 2))

    j_col = 0
    i_row = 0
    #ov_i = 0
    splitFlag = False
    invalidp = 0
    invalidc = 0

    for idx in range(n_pts):
        lat, lon, co2, u, v = lats[idx], lons[idx], co2s[idx], us[idx], vs[idx]

        # distance to previous point - if point is same, reject
        if idx > 0:
            d = lldist(lats[idx-1], lons[idx-1], lat, lon)
        else:
            d = 0.0

        # fill with NaN or real measurement
        if idx>0 and d > thresh:
            # leave as NaN
            invalidp += 1
        else:
            lat_grid[i_row, j_col] = lat
            lon_grid[i_row, j_col] = lon
            co2_grid[i_row, j_col] = co2
            u_grid[i_row, j_col] = u
            v_grid[i_row, j_col] = v
            invalidp = 0

            # path vector: from previous row (same column) to current
            if i_row > 0 and not np.isnan(lat_grid[i_row-1, j_col]):
                path_vecs[i_row, j_col] = latlon_to_km_vector(
                    lat_grid[i_row-1, j_col],
                    lon_grid[i_row-1, j_col],
                    lat, lon
                )

            # column vector: from previous column (same row) to current
            if j_col > 0 and not np.isnan(lat_grid[i_row, j_col-1]):
                col_vecs[i_row, j_col] = latlon_to_km_vector(
                    lat_grid[i_row, j_col-1],
                    lon_grid[i_row, j_col-1],
                    lat, lon
                )

        # Overlaps must be whole numbers of rows, so
        #ov_i = i_row

        # advance to next cell
        i_row += 1
        if i_row == rows:
            i_row = 0
            if invalidp > rows:
              invalidc += 1
            # if too many invalid points, split the data
            if invalidc > invalid_max:
              # Create new segment
              # We know we are at the end of a whole number of rows.
              new_seg = grid_segment(lats[idx-ov:], lons[idx-ov:], co2s[idx-ov:], us[idx-ov:], vs[idx-ov:], ov,
                                    rows, min_dist_samples, threshold_multiplier, invalid_max)
              # Split
              splitFlag = True
            else:
              j_col += 1

    # Find the offset correction
    # 1) average col vec by item in row

    # 2) if one seems to be much larger on average, it is where the starting offset is

    # 3) recompute path vector

    # slice away unused columns
    used_cols = j_col + (1 if i_row>0 else 0)

    # single output
    seg = {
            "path": path_vecs [:, :used_cols, :],
            "col": col_vecs  [:, :used_cols, :],
            "lat": lat_grid[:, :used_cols],
            "lon": lon_grid[:, :used_cols],
            "co2": co2_grid[:, :used_cols],
            "u": u_grid[:, :used_cols],
            "v": v_grid[:, :used_cols]
          }

    if splitFlag == True:
      return [seg] + [new_seg] #, ov_i_n
    else:
      return [seg] #, ov_i

def gridding_with_gaps(lats, lons, co2s, us, vs,
                       i_start=0,
                       rows=8,
                       min_dist_samples=16,
                       threshold_multiplier=1.5,
                       chunk_size=100,
                       overlap_rows=5,
                       invalid_max=16):
    """
    Separate the data into grids, sliced if too long or if the distance between points becomes too great (e.g. area with no data)
    """

    n = len(lats)

    # Split into 100 row blocks
    splits = list(range(chunk_size, n, chunk_size)) + [n]

    segments = []
    start = 0
    ov_pts = overlap_rows * rows

    for end in splits:
        seg_lats = lats[start:end]
        seg_lons = lons[start:end]
        seg_co2 = co2s[start:end]
        seg_us = us[start:end]
        seg_vs = vs[start:end]

        segs, i_start = grid_segment(
                                    seg_lats, seg_lons, seg_co2, seg_us, seg_vs,
                                    i_start, ov_pts,
                                    rows, min_dist_samples, threshold_multiplier, invalid_max
                                    )

        segments.append(segs)

        # next segment considering overlap
        start = end - ov_pts

    return segments

Testing the function

filestr = '/content/drive/MyDrive/CO2_absorption_project/Code/Pirouette/jan19.h5'

# Import h5py to read h5 files
import h5py
import numpy as np

with h5py.File(filestr, 'r') as file:
  co2_full = file['/RetrievalResults/xco2'][()];
  lat_full = file['/RetrievalGeometry/retrieval_latitude'][()];
  lon_full = file['/RetrievalGeometry/retrieval_longitude'][()];
  u_full = file['/RetrievalResults/wind_speed_u_met'][()]; # Eastward wind
  v_full = file['/RetrievalResults/wind_speed_v_met'][()]; # Northward wind
  qflag = file['/RetrievalResults/outcome_flag'][()];

# Only retain HQ points
co2_full = co2_full[qflag==1]
lat_full = lat_full[qflag==1]
lon_full = lon_full[qflag==1]
u_full = u_full[qflag==1]
v_full = v_full[qflag==1]

# Need to manually put in first grid value!!

# Grid
segs = gridding_with_gaps(lat_full, lon_full, co2_full, u_full, v_full, i_start = ???)
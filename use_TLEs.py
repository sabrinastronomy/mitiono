import numpy as np
import ephem
import time
import datetime
from beam_map_sbf import GetBeamMap
from matplotlib import pyplot as plt
import cartopy.crs as ccrs
from urllib import request

def read_latest_tles(url='https://celestrak.com/NORAD/elements/orbcomm.txt'):
    """
    This script reads a file with TLEs of satellites (two line elements) which help us calculate where the satellite
    will be at a given time.
    :param url:
    :return:
    """
    raw = request.urlopen(url)
    lines = [line.decode('utf-8').strip() for line in raw.readlines()]
    nsat = len(lines)//3 # three elements per satellite
    tles = [None]*nsat
    for i in range(nsat):
        tles[i] = lines[3*i:3*(i+1)]
    return tles

def ctime2mjd(tt=None):
    """
    This script converts ctime to mjd.
    """
    if tt is None:
        tt = time.time()
    jd = tt/86400+2440587.5
    return jd-2415020

def check_sat_in_beam(utc_time):
    alt_visible_DRAO, az_visible_DRAO, UTC_times = get_sat_at_DRAO(utc_time)
    alt_az_coord = np.transpose([alt_visible_DRAO, az_visible_DRAO])
    alt_beam = np.full(len(alt_visible_DRAO), 81.41)
    az_beam = np.full(len(alt_visible_DRAO), 180)
    alt_az_beam = np.transpose([alt_beam, az_beam])
    bool_tups = np.abs(alt_az_coord-alt_az_beam) < (2, 2)
    # print(np.abs(alt_az_coord-alt_az_beam))
    # print(bool_tups)
    for b in bool_tups:
        if b[0] and b[1]:
            print("There's a sat within 2 deg of the center of the beam.")
            return

def get_sat_at_DRAO(utc_time, sat=""):
    tles = read_latest_tles('https://celestrak.com/NORAD/elements/gnss.txt')  # This contains TLE files from GNSS
    lat_DRAO = '49.3215'
    lon_DRAO = '-119.6200'  # EAST IS NEGATIVE


    observer_DRAO = ephem.Observer()
    # d = datetime.datetime.utcfromtimestamp(time.time())
    d = ephem.date(utc_time)  # USE UTC TIME
    observer_DRAO.lat = lat_DRAO
    observer_DRAO.lon = lon_DRAO
    observer_DRAO.date = d

    curr_visible_DRAO = []
    alt_visible_DRAO = []
    az_visible_DRAO = []
    plt.clf()
    ax = plt.axes(projection=ccrs.PlateCarree())
    ax.stock_img()

    for i, tle in enumerate(tles):
        tle_rec = ephem.readtle(*tle)
        tle_rec.compute(observer_DRAO)
        alt = tle_rec.alt * 180 / np.pi
        az = tle_rec.az * 180 / np.pi
        # lat = tle_rec.sublat * 180 / np.pi
        # lon = tle_rec.sublong * 180 / np.pi

        # plot_data, = ax.plot(lon, lat, 'r.')
        if alt > 0:
            if sat != "" and sat in tle_rec.name:
                print(tle_rec.name)
                return alt, az, utc_time
            curr_visible_DRAO.append(tle_rec.name)
            alt_visible_DRAO.append(alt)
            az_visible_DRAO.append(az)
    return alt_visible_DRAO, az_visible_DRAO, utc_time

if __name__ == "__main__":
    directories = ["/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/parsed_data/"]
    sat_dict = np.load(directories[0] + "june14_satellite_dict_all.npy", allow_pickle=True).item()
    # new = GetBeamMap()
    # new.replace_dictionary(sat_dict)
    # sat = sat_dict["E33"]
    # times_cno, times_elevs, cnos, elevs, az = sat["times"][0], sat["elev_time"][0], sat["cno"][0], sat["elevations"][0], sat["azimuths"][0]
    # times_elevs, cnos, elevs, az = GetBeamMap.match_elevs(times_cno, times_elevs, cnos, elevs, az)
    # WNC = np.full(len(times_elevs), 2214, dtype=float)
    # UTC_times = GetBeamMap.convert_GPStime_wrapper(times_elevs, WNC)
    # print(UTC_times[1])
    # print(elevs[1], az[1])
    # alt, az, utc_time = get_sat_at_DRAO(UTC_times[1], "PRN E33")
    d = datetime.datetime.utcfromtimestamp(time.time())
    check_sat_in_beam(d)


# lat = np.empty(len(tles))
# lon = np.empty(len(tles))
#
# plt.clf()
# ax = plt.axes(projection=ccrs.PlateCarree())
# gps_sat_names = []
# for i, tle in enumerate(tles):
#     tle_rec = ephem.readtle(*tle)
#     name = tle_rec.name
#     name = "GPS {}".format(name[:-1])
#     gps_sat_names.append(name)
# ax.stock_img()
# tt = time.time()
# color_list = ['b', 'g']

# These are the lat/lon of the CHIME Outrigger sites.
# lat_GBO = 38.4314
# lon_GBO = 79.8192





# a = np.where(curr_visible_DRAO == curr_visible_DRAO[100])

# print(curr_visible_DRAO)
# print(alt_visible_DRAO)

# observer_GBO = ephem.Observer()
# observer_GBO.lat = lat_GBO
# observer_GBO.lon = lon_GBO
#
# observer_HC = ephem.Observer()
# observer_HC.lat = lat_HC
# observer_HC.lon = lon_HC

# for j in range(24):
#     # tt += 3600*j
#     # djd = ctime2mjd()
#     observer_DRAO.date = '2022/6/28 18:13'
#     # observer_GBO.date = djd
#     # observer_HC.date = djd
#
#     t1 = time.time()
#     plt.ion()
#
#     curr_visible_GBO = []
#     curr_visible_DRAO = []
#     curr_visible_HC = []
#
#     alt_visible_GBO = []
#     alt_visible_DRAO = []
#     alt_visible_HC = []
#
#     for i, tle in enumerate(tles):
#         tle_rec = ephem.readtle(*tle)
#
#         # tle_rec.compute(observer_GBO)
#         # if tle_rec.alt > 0:
#         #     curr_visible_GBO.append(tle_rec.name)
#         #     alt_visible_GBO.append(tle_rec.alt)
#
#         tle_rec.compute(observer_DRAO)
#         if tle_rec.alt > 0:
#             curr_visible_DRAO.append(tle_rec.name)
#             alt_visible_DRAO.append(tle_rec.alt)
#
#         # tle_rec.compute(observer_HC)
#         # if tle_rec.alt > 0:
#         #     curr_visible_HC.append(tle_rec.name)
#         #     alt_visible_HC.append(tle_rec.alt)
#
#     overlapping_alt_DRAO_GBO = []
#     overlapping_DRAO_GBO = []
#     for i, vis in enumerate(curr_visible_GBO):
#         if alt_visible_GBO[i] > 2 and vis in curr_visible_DRAO:
#             overlapping_DRAO_GBO.append(vis)
#             overlapping_alt_DRAO_GBO.append((alt_visible_GBO[i]))
#
#     overlapping_alt_DRAO_HC = []
#     overlapping_DRAO_HC = []
#
#     for i, vis in enumerate(curr_visible_HC):
#         if alt_visible_HC[i] > 2 and vis in curr_visible_DRAO:
#             overlapping_DRAO_HC.append(vis)
#             overlapping_alt_DRAO_HC.append((alt_visible_HC[i]))
#     print(curr_visible_DRAO)
#     # print(alt_visible_DRAO)
#     break
# time.sleep(10)
# print('exited')
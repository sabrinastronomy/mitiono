import numpy as np
import ephem
import time
import datetime
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
    with request.urlopen(url) as raw:
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
    curr_visible_DRAO, alt_visible_DRAO, az_visible_DRAO, UTC_times = get_sat_at_DRAO(utc_time)
    alt_az_coord = np.transpose([alt_visible_DRAO, az_visible_DRAO])
    alt_beam = np.full(len(alt_visible_DRAO), 81.41)
    az_beam = np.full(len(alt_visible_DRAO), 180)
    alt_az_beam = np.transpose([alt_beam, az_beam])
    bool_tups = np.abs(alt_az_coord-alt_az_beam) < (5, 5)
    # print(np.abs(alt_az_coord-alt_az_beam))
    # print(bool_tups)
    sat_in_beam = []
    for i, b in enumerate(bool_tups):
        if b[0] and b[1]:
            print(curr_visible_DRAO[i])
            # adding satellite to sat in beam list
            sat_in_beam.append(curr_visible_DRAO[i])
    return sat_in_beam

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
            if sat != "" and sat in tle_rec.name: # getting particular satellite if it's above alt at given UTC_TIME
                print(tle_rec.name)
                return alt, az, utc_time
            curr_visible_DRAO.append(tle_rec.name)
            alt_visible_DRAO.append(alt)
            az_visible_DRAO.append(az)
    return curr_visible_DRAO, alt_visible_DRAO, az_visible_DRAO, utc_time

def get_sat_traj(d_start, minutes, sat_name):
    tles = read_latest_tles('https://celestrak.com/NORAD/elements/gnss.txt')  # This contains TLE files from GNSS
    lat_DRAO = '49.3215'
    lon_DRAO = '-119.6200'  # EAST IS NEGATIVE
    alt_visible_DRAO = []
    az_visible_DRAO = []
    for min in np.arange(minutes, step=1):
        new_time = d_start + datetime.timedelta(minutes=int(min))
        observer_DRAO = ephem.Observer()
        d = ephem.date(new_time)  # USE UTC TIME
        observer_DRAO.lat = lat_DRAO
        observer_DRAO.lon = lon_DRAO
        observer_DRAO.date = d
        for i, tle in enumerate(tles):
            tle_rec = ephem.readtle(*tle)
            if tle_rec.name == sat_name:
                tle_rec.compute(observer_DRAO)
                alt = tle_rec.alt * 180 / np.pi
                az = tle_rec.az * 180 / np.pi
                alt_visible_DRAO.append(alt)
                az_visible_DRAO.append(az)
                break
    return alt_visible_DRAO, az_visible_DRAO

if __name__ == "__main__":
    d = datetime.datetime(2022, 8, 29, 16, 0)
    total_minutes = (1 * 8 * 60) # minutes in the next 14 days
    vis_sats = []
    min_sat_vis = []
    minutes = np.arange(total_minutes, step=1)
    # alt, az = get_sat_traj(d, total_minutes, 'GPS BIIR-11 (PRN 19)')
    # import pandas as pd
    # new_times = []
    # for min in minutes:
    #     new_times.append(d + datetime.timedelta(minutes=int(min)))
    # df = pd.DataFrame({'sat_name': 'GPS BIIR-11 (PRN 19)',
    #                     'UTC time': new_times,
    #                     'alt_at_DRAO': alt,
    #                     'az_at_DRAO': az})
    # df.to_csv("locs_prn19.csv")

    # for min in minutes:
    #     print(min)
    #     new_time = d + datetime.timedelta(minutes=int(min))
    #     sats = check_sat_in_beam(new_time)
    #     if len(sats) > 0:
    #         vis_sats.append(sats)
    #         min_sat_vis.append(new_time)
    import pandas as pd
    # df = pd.DataFrame({'date': min_sat_vis,
    #                    'visible satellites': vis_sats})
    # df.to_csv("8_28_precise_beam.csv")

    df = pd.load_csv("8_28_precise_beam.csv")


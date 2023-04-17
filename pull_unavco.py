import numpy as np
import urllib.request
import georinex as gr
import os

def get_unavco_name(year=2022,day=80,station='drao', hour=None, min=None):
    year=repr(year)
    assert(len(year)==4)
    day=f'{day:03d}'
    dir = f'cacsb.nrcan.gc.ca/gps/data/hrdata/22{day}/22d/{hour}/'
    # dir='data.unavco.org/archive/gnss/rinex/obs/'+year+'/'+day+'/'
    # fname=station+day+'0.'+year[-2:]+'d.Z'
    fname = f"DRAO00CAN_R_{year}{day}{hour}{min}_15M_01S_MO.crx.gz"
    url='https://'+dir+fname
    return url,fname

    
def download_days(year, days, station='drao', outdir='archive', uncompress=True):
    for day in days:
        for hour in np.arange(0, 24, 1):
            hour = str(hour)
            if len(hour) < 2:
                hour = "0" + hour
            print(hour)
            for min in ["00", "15", "30", "45"]:
                url, fname = get_unavco_name(year,day,station, hour=hour, min=min)
                try:
                    print("downloading "+fname)
                    urllib.request.urlretrieve(url, outdir+'/'+fname)
                except:
                    print('error downloading ',url)
                    continue
                if uncompress:
                    os.system('uncompress '+outdir+'/'+fname)
        

station='drao'
year=2022
day=180

url,fname=get_unavco_name(year,day,station)
download_days(2022,[241,242,243,244])

assert(1==0)


#dir='data.unavco.org/archive/gnss/highrate/1-Hz/rinex/'
#dir=dir+repr(year)+'/'+f'{day:03d}/'+station+'/'
#fname=station+f'{day:03d}0.'+f'{year:04d}'[-2:]+'d.Z'
#url='https://'+dir+fname
#print(url)
#https://data.unavco.org/archive/gnss/rinex/obs/2022/180/drao1800.22d.Z

dy = f'{day:03d}'
dir = 'data.unavco.org/archive/gnss/rinex/obs/'+repr(year)+'/'+dy+'/'
fname = station+dy+'0.'+f'{year:04d}'[-2:]+'d.Z'
url = 'https://'+dir+fname
print(url)

outname = 'archive/'+fname
urllib.request.urlretrieve(url, 'archive/'+fname)
dat = gr.load(outname)

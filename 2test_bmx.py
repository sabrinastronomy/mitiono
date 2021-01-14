import bmxdata
import numpy as np
import matplotlib.pyplot as plt
from scipy import signal 
from numpy.fft import rfft

file_path = "/project/s/sievers/sievers/bmx/www.cosmo.bnl.gov/www/bmx/GPS_rings/201107_183400_D2.ring"

data = bmxdata.BMXRingbuffer(file_path)

termchan = data.datad1c2
D2data = [bmxdata.BMXFile(fn.replace('D1','D2')) for fn in flist]
sychan =  np.vstack([d.data['chan3_0'] for d in D2data])
volt = np.hstack([d.data['lj_voltage'] for d in D2data])
ljstate = np.hstack([d.data['lj_diode'] for d in D2data])

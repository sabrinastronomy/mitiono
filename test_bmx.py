import bmxdata
import numpy as np
import matplotlib.pyplot as plt
import scipy
from numpy.fft import rfft

file_path = "/project/s/sievers/sievers/bmx/www.cosmo.bnl.gov/www/bmx/GPS_rings/201107_183400_D2.ring"

data = bmxdata.BMXRingbuffer(file_path)
chan1_data = data.datad0c1[0:200]
len_chan1_data = len(chan1_data)
print(len_chan1_data)
ff = np.fft.rfftfreq(len_chan1_data)
ffted = np.fft.rfft(chan1_data)
print("ff: {}".format(ff))
print("ffted: {}".format(ffted))

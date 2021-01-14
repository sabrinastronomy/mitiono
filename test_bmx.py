import bmxdata
import numpy as np
import matplotlib.pyplot as plt
import scipy
from numpy.fft import rfft

file_path = "/project/s/sievers/sievers/bmx/www.cosmo.bnl.gov/www/bmx/GPS_rings/201107_183400_D2.ring"

data = bmxdata.BMXRingbuffer(file_path).datad0c1
num_samples = len(data)
samp_freq = 1.1e9 #gigasamples
time_per_samples = num_samples/samp_freq

samp_time = 1/samp_freq

chan1_data = data[0:200]
print(len(chan1_data))
len_chan1_data = len(chan1_data)
print(len_chan1_data)
ff = np.fft.rfftfreq(len_chan1_data, samp_time) 
print("ff (pre fftshift): {}".format(ff))
ff = np.fft.fftshift(ff)
ffted = np.fft.rfft(chan1_data)
ffted = np.fft.fftshift(ffted)
plt.plot(ff, ffted)
plt.savefig("fft ex.png")
print("ff: {}".format(ff))
print("ffted: {}".format(ffted))

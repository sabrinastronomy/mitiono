import numpy as np
import bmxdata
#import mkfftw
import scipy
import time
from scipy import fft
from matplotlib import pyplot as plt
import numba as nb
import get_prn_baseband

try:
    print('len(d1) is ',len(d1))
except:
    dat=bmxdata.BMXRingbuffer('201107_024900_D1.ring');tag1='D1';
    d1=dat.datad1c1;tag1=tag1+'d1';



f=1100e6
nsamp=np.int(f/get_prn_baseband.f_prn*1023)*4

mymax=np.zeros(32)
for sat in range(1,33):
    myprn=get_prn_baseband.get_prn_nb(f,nsamp,sat)
    mycorr=fft.irfft(fft.rfft(d1[:nsamp])*np.conj(fft.rfft(myprn)))
    mymax[sat-1]=np.max(np.abs(mycorr))
    print('max value on ',sat,' is ',mymax[sat-1])


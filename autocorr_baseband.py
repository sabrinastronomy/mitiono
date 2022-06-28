"""
This is my sad attempt at trying to do an autocorr with the Airspy data and being confused why the
output data type the Airspy claims.
"""

from scipy import signal
import seaborn as sns
import numpy as np
from matplotlib import pyplot as plt

# Key to end of file
# 0=FLOAT32_IQ, 1=FLOAT32_REAL, 2=INT16_IQ(default), 3=INT16_REAL, 4=U16_REAL, 5=RAW

sns.set_theme(style="whitegrid")

def reader(path2data):
    ## read in file, I'm confused which dtype to use here...
    return np.fromfile(path2data, dtype='int16')

## trying 1=FLOAT32_REAL
path2data = '/Users/sabrinaberger/Library/Mobile Documents/com~apple~CloudDocs/mitiono/D3A(6)_Data/d3a_airspy/airspy_C2_LO_1574MHz_1'
d1 = reader(path2data)

def autocorrelate(data):
    """
    Jon's script that I am confused about how to use
    """
    fsamp = 1.174e9
    nchan = 16384 * 2 # number of samples
    nn = 2 * nchan # two times the number of samples
    x = np.linspace(0, 2 * np.pi, nn) # number of samples over 2pi
    win = 0.5 - 0.5 * np.cos(x) # window function
    nchunk = 10000 # number of chunks to window
    spec = 0
    xspec = 0
    nu0 = fsamp / 1e6 # scaling I think
    nu1 = 1.5 * nu0
    nu = np.linspace(nu0, nu1, nchan + 1)
    sig = 5
    nu_gps = 1176
    filt = np.exp(-0.5 * (nu - nu_gps) ** 2 / sig ** 2)

    ts_back = np.zeros(len(data), dtype='float32')
    d1 = data
    d2 = data
    for i in range(nchunk):
        # windowing stuff and FFTing
        block = win * d1[i * nchan:(i * nchan + nn)]
        tmp = np.fft.rfft(block)
        block2 = win * d2[i * nchan:(i * nchan + nn)]
        tmp2 = np.fft.rfft(block2)
        xspec = xspec + np.abs(tmp * np.conj(tmp2))
        spec = spec + np.real(tmp) ** 2 + np.imag(tmp) ** 2

        tmp = tmp * filt
        tmp = np.fft.irfft(tmp)
        ts_back[i * nchan:(i * nchan + nn)] = ts_back[i * nchan:(i * nchan + nn)] + tmp

    ft1 = np.fft.rfft(ts_back[:10000000])
    acorr = np.fft.irfft(ft1 * np.conj(ft1))

    snr = 60  # this comes from looking at the peak of the autocorr against the std of
    # the autocorr from nearby frequencies

    myoff = 500000;
    ind = np.argmax(acorr[myoff:myoff + 1000000]);
    ind = ind + myoff

    tvec = np.arange(len(acorr)) / fsamp * 1e3  # get the time in microseconds
    plt.clf()
    plt.plot(tvec, acorr)
    plt.xlabel(r'Lag - $\mu$s')
    plt.title('GPS Bandpass Autocorr')
    plt.savefig('gps_autocorr.png')

    plt.axis([999.5, 1000.5, -2.5e6, 2.5e6])
    fac = (snr - 1) / snr
    amp = acorr[ind]
    plt.plot([0, 8000], [amp, amp])
    plt.plot([0, 8000], [amp * fac, amp * fac])
    plt.title('GPS Autocorr Zoom')
    plt.legend(['autocorr', 'Peak amplitude', r'Peak-1$\sigma$'])
    plt.savefig('gps_autocorr_zoom.png')
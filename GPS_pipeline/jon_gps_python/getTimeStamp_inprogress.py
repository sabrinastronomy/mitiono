from scipy import fft
import numpy as np
import matplotlib.pyplot as plt
import get_prn_baseband
import prn
import bmxdata

class TimeStamp:
    """
    This class gets a GPS timestamp from input raw data, i.e., the timestamp of the beginning of a PRN code.
    """
    def __init__(self, sampling_frequency, inputFilename):
        self.sampling_frequency = sampling_frequency
        self.inputFilename = inputFilename # input filename of data to be processed for timestamp
        self.SV = {} # dictionary of GPS satellite numbers with associated taps to generate PRN codes
        self.getChannel1()

    def getChannel1(self):
        data = bmxdata.BMXRingbuffer(self.inputFilename)
        self.dataD1C1 = data.datad1c1
        self.dataD0C1 = data.datad0c1
        return

    def findSatellite(self):
        # Author: Jon Sievers
        nsamp = np.int(self.sampling_frequency / get_prn_baseband.f_prn * 1023) * 4 #????
        mymax = np.zeros(32)
        for sat in range(1, 33):
            myprn = get_prn_baseband.get_prn_nb(self.sampling_frequency, nsamp, sat) # get PRN baseband signal to cross-correlate with
            mycorr = fft.irfft(fft.rfft(self.dataD1C1[:nsamp]) * np.conj(fft.rfft(myprn))) # do the cross-correlation
            mymax[sat - 1] = np.max(np.abs(mycorr))
            print('max value on ', sat, ' is ', mymax[sat - 1])
        self.satNumber = np.argmax(mymax) + 1 # Satellite with highest correlation
        return

    def autocorrelate(self):
        # Author: Jon Sievers
        fsamp = self.sampling_frequency
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

        ts_back = np.zeros(len(self.dataD1C1), dtype='float32')
        d1 = self.dataD0C1
        d2 = self.dataD1C1
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


        f1 = 154 * 10.23e6
        # f1=1450e6
        f2 = f1 - fsamp

        nsamp = 10000
        tvec = np.arange(nsamp) / fsamp

        ts1 = np.exp(2J * np.pi * tvec * f1)
        ts2 = np.exp(2J * np.pi * tvec * f2)
        print(np.std(ts1 - ts2))
        print(f2 / 1e9)

        myk = nsamp * f2 / fsamp
        print(myk)
        myk = int(myk)
        myft = np.fft.fft(ts1)
        print(np.abs(myft)[myk - 20:myk + 20])



if __name__ == "__main__":
    BMXSamplingFrequency = 1100e6 # 1100 MHz
    BMXTest = TimeStamp(BMXSamplingFrequency, "/Users/sabrinaberger/mitiono/GPS_pipeline/201107_174500_D1.ring")
    BMXTest.findSatellite()
    BMXTest.autocorrelate()
    print("Satellite with highest power is: {}".format(BMXTest.satNumber))
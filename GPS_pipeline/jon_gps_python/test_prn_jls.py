import prn
import numpy as np
from matplotlib import pyplot as plt

f = 154 * 10.23e6      # carrier frequency
f_prn = 10.23e6 / 10   # PRN frequency
sample_freq=1100e6


myprn=prn.PRN(24)
myprn=np.asarray(myprn)*2-1


nsamp=100000
tvec=np.arange(nsamp)/sample_freq
carrier=np.sin(2*np.pi*tvec*f)
tvec_prn=np.asarray(np.floor(tvec*f_prn),dtype='int')
tvec_prn=tvec_prn%len(myprn)
ts_prn=myprn[tvec_prn]
gps_sig=ts_prn*carrier
assert(1==0)

# Generate ~1ms of data for a plain carrier, and a GPS signal
signal = []
sine = []
#for i in xrange(10000000):
for i in range(10000):
    t = i*sample_rate
    c = carrier(t)
    sine.append(c)
    signal.append(c * prn(t))

plt.ion()

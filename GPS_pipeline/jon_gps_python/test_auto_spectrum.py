import numpy as np
import bmxdata
from matplotlib import pyplot as plt


try:
    print('len(d1) is ',len(d1))
except:
    dat=bmxdata.BMXRingbuffer('201107_024900_D1.ring');tag1='D1';
    d1=dat.datad0c1;#tag1=tag1+'d1';
    d2=dat.datad1c1;#tag1=tag1+'d1';
    del(dat)




fsamp=1.1e9
nchan=16384
nn=2*nchan
x=np.linspace(0,2*np.pi,nn)
win=0.5-0.5*np.cos(x)
nchunk=10000
spec=0
xspec=0
nu0=fsamp/1e6
nu1=1.5*nu0
nu=np.linspace(nu0,nu1,nchan+1)
sig=5
nu_gps=1176
#nu_gps=1196
filt=np.exp(-0.5*(nu-nu_gps)**2/sig**2)

ts_back=np.zeros(len(d1),dtype='float32')

for i in range(nchunk):
    block=win*d1[i*nchan:(i*nchan+nn)]
    tmp=np.fft.rfft(block)
    block2=win*d2[i*nchan:(i*nchan+nn)]
    tmp2=np.fft.rfft(block2)
    xspec=xspec+np.abs(tmp*np.conj(tmp2))
    spec=spec+np.real(tmp)**2+np.imag(tmp)**2
    
    tmp=tmp*filt
    tmp=np.fft.irfft(tmp)
    ts_back[i*nchan:(i*nchan+nn)]=ts_back[i*nchan:(i*nchan+nn)]+tmp

ft1=np.fft.rfft(ts_back[:10000000])
acorr=np.fft.irfft(ft1*np.conj(ft1))
i0=np.int(1e-3*fsamp)


#plt.clf()
#plt.semilogy(nu,spec)
#plt.plot(nu,xspec)


snr=60 #this comes from looking at the peak of the autocorr against the std of
       #the autocorr from nearby frequencies

myoff=500000;ind=np.argmax(acorr[myoff:myoff+1000000]);ind=ind+myoff

tvec=np.arange(len(acorr))/fsamp*1e6 #get the time in microseconds
plt.clf()
plt.plot(tvec,acorr)
plt.xlabel(r'Lag - $\mu$s')
plt.title('GPS Bandpass Autocorr')
plt.savefig('gps_autocorr.png')

plt.axis([999.5,1000.5,-2.5e6,2.5e6])
fac=(snr-1)/snr
amp=acorr[ind]
plt.plot([0,8000],[amp,amp])
plt.plot([0,8000],[amp*fac,amp*fac])
plt.title('GPS Autocorr Zoom')
plt.legend(['autocorr','Peak amplitude',r'Peak-1$\sigma$'])
plt.savefig('gps_autocorr_zoom.png')


assert(1==0)



f1 = 154 * 10.23e6
#f1=1450e6
f2=f1-fsamp

nsamp=10000
tvec=np.arange(nsamp)/fsamp

ts1=np.exp(2J*np.pi*tvec*f1)
ts2=np.exp(2J*np.pi*tvec*f2)
print(np.std(ts1-ts2))
print(f2/1e9)

myk=nsamp*f2/fsamp
print(myk)
myk=int(myk)
myft=np.fft.fft(ts1)
print(np.abs(myft)[myk-20:myk+20])


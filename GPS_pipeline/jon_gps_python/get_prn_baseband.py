import prn
import numpy as np
import numba as nb



f = 154 * 10.23e6      # carrier frequency
f_prn = 10.23e6 / 10   # C/A PRN frequency
sample_freq=1100e6

def get_prn(freq,nsamp,sat=24):
    myprn=prn.PRN(sat)
    myprn=np.asarray(myprn)*2-1
    tvec=np.arange(nsamp)/freq
    carrier=np.sin(2*np.pi*tvec*f)
    tvec_prn=np.asarray(np.floor(tvec*f_prn),dtype='int')
    tvec_prn=tvec_prn%len(myprn)
    ts_prn=myprn[tvec_prn]
    gps_sig=ts_prn*carrier
    return gps_sig

@nb.njit(parallel=True)
def process_prn_nb(freq,nsamp,myprn):
    vec=np.empty(nsamp)
    finv=1/freq
    fac=2*np.pi*f
    prn_len=len(myprn)
    for i in nb.prange(nsamp):
        t=i*finv
        carrier=np.sin(t*fac)
        t_prn=int(t*f_prn)
        t_prn=t_prn%prn_len
        vec[i]=carrier*myprn[t_prn]
    return vec

def get_prn_nb(freq,nsamp,sat=24):
    myprn=prn.PRN(sat)
    myprn=np.asarray(myprn)*2-1
    return process_prn_nb(freq,int(nsamp),myprn)



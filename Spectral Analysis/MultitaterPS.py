from __future__ import division, print_function, absolute_import
import numpy as np
import os
import pandas as pd
from scipy import signal
from scipy.signal import windows
#from scipy.interpolate import spline
import datetime
import operator
import warnings
from scipy import fftpack, linalg, special
from scipy._lib.six import string_types
import matplotlib.pyplot as plt

def _len_guards(M):
    """Handle small or incorrect window lengths"""
    if int(M) != M or M < 0:
        raise ValueError('Window length M must be a non-negative integer')
    return M <= 1

def _extend(M, sym):
    """Extend window by 1 sample if needed for DFT-even symmetry"""
    if not sym:
        return M + 1, True
    else:
        return M, False

## Code for dpss is from scypy pre-release 

def dpss(M, NW=2.5, Kmax=5, sym=False, norm=None, return_ratios=False):

    if _len_guards(M):
        return np.ones(M)
    if norm is None:
        norm = 'approximate' if Kmax is None else 2
    known_norms = (2, 'approximate', 'subsample')
    if norm not in known_norms:
        raise ValueError('norm must be one of %s, got %s'
                         % (known_norms, norm))
    if Kmax is None:
        singleton = True
        Kmax = 1
    else:
        singleton = False
    Kmax = operator.index(Kmax)
    if not 0 < Kmax <= M:
        raise ValueError('Kmax must be greater than 0 and less than M')
    if NW >= M/2.:
        raise ValueError('NW must be less than M/2.')
    if NW <= 0:
        raise ValueError('NW must be positive')
    M, needs_trunc = _extend(M, sym)
    W = float(NW) / M
    nidx = np.arange(M)

    d = ((M - 1 - 2 * nidx) / 2.) ** 2 * np.cos(2 * np.pi * W)
    e = nidx[1:] * (M - nidx[1:]) / 2.

    w, windows = linalg.eigh_tridiagonal(
        d, e, select='i', select_range=(M - Kmax, M - 1))
    w = w[::-1]
    windows = windows[:, ::-1].T

    fix_even = (windows[::2].sum(axis=1) < 0)
    for i, f in enumerate(fix_even):
        if f:
            windows[2 * i] *= -1    

    thresh = max(1e-7, 1. / M)
    for i, w in enumerate(windows[1::2]):
        if w[w * w > thresh][0] < 0:
            windows[2 * i + 1] *= -1

    if return_ratios:
        dpss_rxx = _fftautocorr(windows)
        r = 4 * W * np.sinc(2 * W * nidx)
        r[0] = 2 * W
        ratios = np.dot(dpss_rxx, r)
        if singleton:
            ratios = ratios[0]
    # Deal with sym and Kmax=None
    if norm != 2:
        windows /= windows.max()
        if M % 2 == 0:
            if norm == 'approximate':
                correction = M**2 / float(M**2 + NW)
            else:
                s = np.fft.rfft(windows[0])
                shift = -(1 - 1./M) * np.arange(1, M//2 + 1)
                s[1:] *= 2 * np.exp(-1j * np.pi * shift)
                correction = M / s.real.sum()
            windows *= correction
    # else we're already l2 normed, so do nothing
    if needs_trunc:
        windows = windows[:, :-1]
    if singleton:
        windows = windows[0]
    return (windows, ratios) if return_ratios else windows

class GetData(object):

    def __init__(  self
                 , datatypes=["S", "B", "nB"]
                 , net_type=None
                 , simulations="all" ):

        self.pre_sim     = "sim-"
        self.datatypes   = datatypes
        self.simulations = simulations
        if net_type in ["equalnet0.05", "equalnet0.1", "betanet", "gammanet"]:
            self.net_type = net_type
        else:
            raise Exception("Please review the net_type provided, it has to be one of the following:\n[equalnet0.05, equalnet0.1, betanet, gammanet]")

        return None

    def getsims(self, files):

        if self.simulations=="all":
            sms=1
            for file in files:
                if "sim-" in file:
                    sms+=1
            sims = np.arange(1, sms, 1)
        elif self.simulations == "last":
            sms=0
            for file in files:
                if "sim-" in file:
                    sms+=1
            sims = np.array([sms])
        else:
            sims = np.asarray(self.simulations, dtype = "int32")

        return sims

    def getraster(self, cell_type, raster_type=None):

        rasters = []
        assert (cell_type in ["pyr", "som", "pv"]), "provide a correct cell type as a string, pyr, som or pv"

        netfolder   = os.getcwd() + '/simulations/' + self.net_type
        files     = os.listdir(netfolder)
        sims      = self.getsims(files)
        path      = "net_results/results_"
        if raster_type == "S":
            file_name = "Sraster_summed.txt"
        elif raster_type == "B":
            file_name = "Braster_summed.txt"
        elif raster_type == "nB":
            file_name = "nBraster_summed.txt"
        else:
            raise Exception("Provide a proper raster_type from [S, B, nB]")
        
        for i in range(len(sims)):
            sim = self.pre_sim + str(sims[i])
            with open(netfolder + "/" + sim + "/" + path + cell_type + "/" + file_name) as rast:
                raster = [[float(num) for num in line.split('\r\n')] for line in rast]
            #raster = np.loadtxt(simdate + "/" + sim + "/" + path + cell_type + "/" + file_name, delimiter='\t')
            rasters.append(np.array(raster).flatten())
        return np.array(rasters)
        
    def gettimes(self):
        times = []
        netfolder = os.getcwd() + '/simulations/' + self.net_type
        files     = os.listdir(netfolder)
        sims      = self.getsims(files)
        path      = "net_results/results_pyr/"

        for i in range(len(sims)):
            sim = self.pre_sim + str(sims[i])    
            with open(netfolder + "/" + sim + "/" + path + "times.txt") as tm:
                time = [[float(num) for num in line.split('\r\n')] for line in tm]
            
            times.append(np.array(time).flatten())
                    
        return np.array(times)
    def __call__(self):
        output = {}
        output["times"]      = self.gettimes()
        output["Srasterpyr"] = self.getraster("pyr", raster_type="S") 
        output["Brasterpyr"] = self.getraster("pyr", raster_type="B") 
        output["nBrasterpyr"]= self.getraster("pyr", raster_type="nB") 
        #output["Srastersom"] = self.getraster("som", raster_type="S") 
        #output["Brastersom"] = self.getraster("som", raster_type="B")     
        #output["nBrastersom"]= self.getraster("som", raster_type="nB")         
        #output["Srasterpv"]  = self.getraster("pv", raster_type="S") 
        #output["Brasterpv"]  = self.getraster("pv", raster_type="B") 
        #output["nBrasterpv"]  = self.getraster("pv", raster_type="nB") 
        return output

# Check the data.

class MTPS(object):

    def __init__(  self
                 , data                # vector of vector
                 , datatimes           # vector
                 , timestep=0.0001
                 , lowpassfreq=200.    # in HZ
                 , highpassfreq=2.     # in Hz
                 , lowpass_order=4     # integer
                 , highpass_order=4    # integer
                 , discard_initial=0.0 # in seconds
                 , discard_last=0.     # in seconds
                 , windowsize=500      # integer
                 , windowstep=50       # integer   
                 , smooth=False        # Bool
                 , smoothinterval=1    # integer
                 , padsize=256         # int and power of 2
                 , multitaper=True
                 , lastfreq=200.
                 ):  
        self.data     = data
        self.datat    = datatimes
        self.trials   = len(data)
        self.timestep = timestep
        self.fs       = 1/timestep
        self.lpassf   = lowpassfreq
        self.hpassf   = highpassfreq
        self.lorder   = lowpass_order
        self.horder   = highpass_order
        self.throw0   = discard_initial
        self.throwf   = discard_last
        self.wsize    = windowsize
        self.wstep    = windowstep
        self.wdata    = []
        self.PS       = []
        self.PS_sig   = None
        self.freqs    = None
        self.smooth   = smooth
        self.smhint   = smoothinterval
        self.padsize  = padsize
        self.lastfreq = lastfreq

    def select_data(self):
   
        t0        = self.throw0 
        tf        = self.datat[-1] - self.throwf 
        idxt0     = int(t0 / self.timestep) 
        idxtf     = int(tf / self.timestep)
        self.data = list(self.data)
        for d in range(len(self.data)):
            self.data[d] = self.data[d][idxt0:idxtf]
        self.datat = self.datat[idxt0:idxtf]


    def butterfly_lowpass(self):

        nyq           = 0.5 * self.fs
        cutoff        = self.lpassf
        normal_cutoff = cutoff / nyq
        b, a          = signal.butter(self.lorder, normal_cutoff, btype='low', analog=False)
        self.data     = signal.filtfilt(b, a, self.data)


    def butterfly_highpass(self):

        nyq           = 0.5 * self.fs
        cutoff        = self.hpassf
        normal_cutoff = cutoff / nyq
        b, a          = signal.butter(self.horder, normal_cutoff, btype='high', analog=False)
        self.data     = signal.filtfilt(b, a, self.data)


    def windowing(self):

        steps         = len(self.data[0])
        totalwindows  = int( (steps - self.wsize)/self.wstep )
        for i in range(self.trials):
            wdata = []
            for w in range(totalwindows):
                wpos0     = w * self.wstep
                wposf     = wpos0 + self.wsize
                prewindow = self.data[i][wpos0:wposf]
                prewindow = prewindow - np.mean(prewindow)
                wdata.append(prewindow/np.max(prewindow))
            self.wdata.append(wdata)


    def paddata(self):

        if self.padsize == 0:
            pass
        else:
            while len(self.wdata[0][0]) > self.padsize: # pad always to powers of 2
                self.padsize *= 2
        halfpad = (self.padsize - len(self.wdata[0][0]))/2
        halfpad = int(halfpad)
        for i in range(len(self.wdata)):
            for w in range(len(self.wdata[0])):
                predata = np.append(np.zeros(halfpad), np.array(self.wdata[i][w]))
                self.wdata[i][w] = np.append(predata, np.zeros(halfpad))

    def slepianfunc(self, NW=3, k=5):

        M   = self.padsize
        win = dpss(M, NW=NW, Kmax=k, return_ratios=False)
        return win
        

    def smoothing(self, data):

        for j in range(self.smhint, len(data) - self.smhint ):
            point = np.average(data[(j-self.smhint):(j+self.smhint)])
            data[j] = point

        return data

    def PowerSpectrum(self, NW, k, allfreqs=False):

        slepians = self.slepianfunc(NW=NW, k=k)
        K        = len(slepians)
        pre_PS = []
        pre_time_freq = []
        for tr in range(len(self.wdata)):
            PS_tr = [] 
            for w in range(len(self.wdata[0])):
                ps = []
                for k in range(K):
                    preps = np.abs( np.fft.fft(self.wdata[tr][w]*slepians[k]) )**2 
                    ps.append(preps)
                ps = np.average(ps, axis=0)
                PS_tr.append(ps)
            pre_time_freq.append(PS_tr)
            pre_PS.append(np.average(PS_tr, axis=0))
        
        time_freq = np.average(pre_time_freq, axis=0)
        time_freq_sig = np.std(pre_time_freq, axis=0)
        PS_sig = np.std(pre_PS, axis=0)
        PS     = np.average(pre_PS, axis=0)
        if self.padsize == 0:
            f_len = self.wsize
        else:
            f_len = self.padsize
        freqs  = np.fft.fftfreq(f_len, self.timestep)
        idx    = np.argsort(freqs)

        pfreqs = idx[int((len(idx)/2)):]

        if allfreqs:
            freqs = freqs[idx]
            PS    = PS[idx]
        else:
            freqs = freqs[pfreqs]
            PS    = PS[pfreqs]

        if self.lastfreq == None:
            freqlim = len(freqs)
        else:
            freqlim = np.where(np.floor(freqs) >= float(self.lastfreq))[0][0]
        self.tf    = np.asarray(time_freq)[:,:freqlim].T
        self.tf_sig = np.asarray(time_freq_sig)[:,:freqlim].T
        self.PS_sig = PS_sig[:freqlim]
        self.PS    = PS[:freqlim]
        self.freqs = freqs[:freqlim] 

    def __call__(self, NW=2.5, k=5):

        self.select_data()
        self.butterfly_lowpass()
        self.windowing()
        self.paddata()
        self.PowerSpectrum(NW, k)

def plot_ps(rasterPS):
    my_dpi = 96
    plt.figure(figsize=(600/my_dpi, 300/my_dpi), dpi=my_dpi)
    ps_plot, = plt.plot(rasterPS.freqs, rasterPS.PS, color='black')
    plt.xlabel('frequency (Hz)', fontsize=15)
    plt.tick_params(axis="x", labelsize=15)
    plt.ylabel('PSD', fontsize=15)
    peak_idx = np.argmax(rasterPS.PS)
    peak_freq_S = rasterPS.freqs[peak_idx]
    x_ticks = np.linspace(0,100, num=6)
    x_ticks = np.append(x_ticks, peak_freq_S)
    plt.plot(np.ones(len(rasterPS.freqs))*rasterPS.freqs[peak_idx], np.linspace(0,rasterPS.PS[peak_idx],num=len(rasterPS.freqs)), '--', c='black')
    print('peak freq = ' + '%0.1f' % peak_freq_S)
    plt.tick_params(axis="y", labelsize=15)
    plt.fill_between(rasterPS.freqs, rasterPS.PS+rasterPS.PS_sig, rasterPS.PS-rasterPS.PS_sig, facecolor='orange', alpha=0.5)
    plt.tight_layout()
    #plt.savefig("/Users/pau/Desktop/MSc_CNS/Thesis/CNS_2019_BCN/Poster/Figures/gammaPS.png")
    return ps_plot

def plot_tf(rasterPS, times):
    my_dpi = 96
    plt.figure(figsize=(1200/my_dpi, 500/my_dpi), dpi=my_dpi)
    times_f = np.linspace(times[0], times[-1], num=rasterPS.tf.shape[1])/1000
    xv, yv = np.meshgrid(times_f, rasterPS.freqs)
    #peak_freqs_idx = np.argmax(rasterPS.tf, axis=0)
    #peak_freqs_S = rasterPS.freqs[peak_freqs_idx]
    plt.plot(5.5*np.ones((len(times_f),)), np.linspace(rasterPS.freqs[0], rasterPS.freqs[-1], num=len(times_f)), 'k', '--',  linewidth=2)
    plt.ylim(0.5, 10)
    tf_plot = plt.contourf(xv, yv, rasterPS.tf)
    plt.xlabel('time (s))', fontsize=15)
    plt.tick_params(axis="x", labelsize=15)
    plt.ylabel('frequency (Hz)', fontsize=15)
    plt.tick_params(axis="y", labelsize=15)
    #plt.plot(times_f, peak_freqs_S, 'k',  linewidth=1)
    plt.yticks([0, 20, 40, 60, 80], [0, 20, 40, 60, 80])
    plt.tick_params(axis="y", labelsize=15)
    plt.tight_layout()
    #plt.xlim((1000.0,times_f[-1]))
    return tf_plot



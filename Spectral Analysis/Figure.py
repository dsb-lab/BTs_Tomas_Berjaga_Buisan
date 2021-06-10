#!/usr/bin/env python2
# -*- coding: utf-8 -*-
#import os
import MultitaterPS as mtps
import matplotlib.pyplot as plt
import numpy as np
hfont = {'fontname':'Helvetica'}
my_dpi= 96


#Please consider that to run this file, you must have first saved some network simulations.
#The file is prepared for the analysis of Multitaper Power Spectrum for 10 networks for different AIS T-type Conductances.
#The only change that is needed to be functional is to change the Directories.

'''
with open("SOM NETWORK/g=0/I=250/PSC_t_Inject_250_SOM_1.txt") as rast:
    raster1 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=0/I=250/PSC_t_Inject_250_SOM_2.txt") as rast:
    raster2 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=0/I=250/PSC_t_Inject_250_SOM_3.txt") as rast:
    raster3 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=0/I=250/PSC_t_Inject_250_SOM_4.txt") as rast:
    raster4 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=0/I=250/PSC_t_Inject_250_SOM_5.txt") as rast:
    raster5 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=0/I=250/PSC_t_Inject_250_SOM_6.txt") as rast:
    raster6 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=0/I=250/PSC_t_Inject_250_SOM_7.txt") as rast:
    raster7 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=0/I=250/PSC_t_Inject_250_SOM_8.txt") as rast:
    raster8 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=0/I=250/PSC_t_Inject_250_SOM_9.txt") as rast:
    raster9 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=0/I=250/PSC_t_Inject_250_SOM_10.txt") as rast:
    raster10 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]

raster_full=[raster1,raster2,raster3,raster4,raster5,raster6,raster7,raster8,raster9,raster10]
x = np.arange(0.5, 3+0.0001,  0.0001)

rasterBpyrPS = mtps.MTPS(np.asarray(raster_full), x, timestep=0.0001, windowsize=1000, windowstep=10, highpassfreq=1.0, lowpassfreq=100.0, lastfreq=100, padsize=4048)
rasterBpyrPS(NW=2.5, k=5)

with open("SOM NETWORK/g=25/I=250/PSC_t_Inject_250_SOM_1.txt") as rast:
    raster1 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=25/I=250/PSC_t_Inject_250_SOM_2.txt") as rast:
    raster2 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=25/I=250/PSC_t_Inject_250_SOM_3.txt") as rast:
    raster3 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=25/I=250/PSC_t_Inject_250_SOM_4.txt") as rast:
    raster4 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=25/I=250/PSC_t_Inject_250_SOM_5.txt") as rast:
    raster5 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=25/I=250/PSC_t_Inject_250_SOM_6.txt") as rast:
    raster6 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=25/I=250/PSC_t_Inject_250_SOM_7.txt") as rast:
    raster7 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=25/I=250/PSC_t_Inject_250_SOM_8.txt") as rast:
    raster8 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=25/I=250/PSC_t_Inject_250_SOM_9.txt") as rast:
    raster9 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=25/I=250/PSC_t_Inject_250_SOM_10.txt") as rast:
    raster10 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]

raster_full=[raster1,raster2,raster3,raster4,raster5,raster6,raster7,raster8,raster9,raster10]
x = np.arange(0.5, 3+0.0001,  0.0001)

rasterBpyrPS2 = mtps.MTPS(np.asarray(raster_full), x, timestep=0.0001, windowsize=1000, windowstep=10, highpassfreq=1.0, lowpassfreq=100.0, lastfreq=100, padsize=4048)
rasterBpyrPS2(NW=2.5, k=5)

with open("SOM NETWORK/g=75/I=250/PSC_t_Inject_250_SOM_1.txt") as rast:
    raster1 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=75/I=250/PSC_t_Inject_250_SOM_2.txt") as rast:
    raster2 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=75/I=250/PSC_t_Inject_250_SOM_3.txt") as rast:
    raster3 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=75/I=250/PSC_t_Inject_250_SOM_4.txt") as rast:
    raster4 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=75/I=250/PSC_t_Inject_250_SOM_5.txt") as rast:
    raster5 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=75/I=250/PSC_t_Inject_250_SOM_6.txt") as rast:
    raster6 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=75/I=250/PSC_t_Inject_250_SOM_7.txt") as rast:
    raster7 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=75/I=250/PSC_t_Inject_250_SOM_8.txt") as rast:
    raster8 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=75/I=250/PSC_t_Inject_250_SOM_9.txt") as rast:
    raster9 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=75/I=250/PSC_t_Inject_250_SOM_10.txt") as rast:
    raster10 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]

raster_full=[raster1,raster2,raster3,raster4,raster5,raster6,raster7,raster8,raster9,raster10]
x = np.arange(0.5, 3+0.0001,  0.0001)

rasterBpyrPS3 = mtps.MTPS(np.asarray(raster_full), x, timestep=0.0001, windowsize=1000, windowstep=10, highpassfreq=1.0, lowpassfreq=100.0, lastfreq=100, padsize=4048)
rasterBpyrPS3(NW=2.5, k=5)

with open("SOM NETWORK/g=150/I=250/PSC_t_Inject_250_SOM_1.txt") as rast:
    raster1 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=150/I=250/PSC_t_Inject_250_SOM_2.txt") as rast:
    raster2 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=150/I=250/PSC_t_Inject_250_SOM_3.txt") as rast:
    raster3 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=150/I=250/PSC_t_Inject_250_SOM_4.txt") as rast:
    raster4 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=150/I=250/PSC_t_Inject_250_SOM_5.txt") as rast:
    raster5 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=150/I=250/PSC_t_Inject_250_SOM_6.txt") as rast:
    raster6 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=150/I=250/PSC_t_Inject_250_SOM_7.txt") as rast:
    raster7 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=150/I=250/PSC_t_Inject_250_SOM_8.txt") as rast:
    raster8 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=150/I=250/PSC_t_Inject_250_SOM_9.txt") as rast:
    raster9 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]
with open("SOM NETWORK/g=150/I=250/PSC_t_Inject_250_SOM_10.txt") as rast:
    raster10 = np.transpose([[float(num) for num in line.split('\r\n')] for line in rast])[0]

raster_full=[raster1,raster2,raster3,raster4,raster5,raster6,raster7,raster8,raster9,raster10]
x = np.arange(0.5, 3+0.0001,  0.0001)

rasterBpyrPS4 = mtps.MTPS(np.asarray(raster_full), x, timestep=0.0001, windowsize=1000, windowstep=10, highpassfreq=1.0, lowpassfreq=100.0, lastfreq=100, padsize=4048)
rasterBpyrPS4(NW=2.5, k=5)

plt.figure(figsize=(1200/my_dpi, 600/my_dpi), dpi=my_dpi)
plt.plot(rasterBpyrPS.freqs, rasterBpyrPS.PS, color='orange',label='$g_T$ = 0 nS')
plt.plot(rasterBpyrPS2.freqs, rasterBpyrPS2.PS, color='red',label='$g_T$ = 25 nS')
plt.plot(rasterBpyrPS3.freqs, rasterBpyrPS3.PS, color='blue',label='$g_T$ = 75 nS')
plt.plot(rasterBpyrPS4.freqs, rasterBpyrPS4.PS, color='green',label='$g_T$ = 150 nS')
plt.xlabel('Frequency (Hz)', fontsize=20)
plt.tick_params(axis="x", labelsize=20)
plt.ylabel('PSD', fontsize=20)
plt.legend(loc='best',prop={'size': 12})
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
#plt.title('Burst events power spectrum', fontsize=12)
peak_idx1 = np.argmax(rasterBpyrPS.PS)
peak_freq_B1 = rasterBpyrPS.freqs[peak_idx1]
print("Low peak at: %.1f Hz" % rasterBpyrPS.freqs[peak_idx1])

peak_idx2 = np.argmax(rasterBpyrPS.PS[0:int(round(len(rasterBpyrPS.PS)/4))])
peak_freq_B2 = rasterBpyrPS.freqs[peak_idx2]
print("High peak at: %.1f Hz" % rasterBpyrPS.freqs[peak_idx2])

x_ticks = np.linspace(0,100, num=6)
x_ticks = np.append(x_ticks, peak_freq_B1)
#plt.plot(np.ones(len(rasterBpyrPS.freqs))*rasterBpyrPS.freqs[peak_idx1], np.linspace(0,rasterBpyrPS.PS[peak_idx1],num=len(rasterBpyrPS.freqs)), '--', c='black')
#plt.plot(np.ones(len(rasterBpyrPS.freqs))*rasterBpyrPS.freqs[peak_idx2], np.linspace(0,rasterBpyrPS.PS[peak_idx2],num=len(rasterBpyrPS.freqs)), '--', c='black')

plt.tick_params(axis="y", labelsize=15)
plt.fill_between(rasterBpyrPS.freqs, rasterBpyrPS.PS+rasterBpyrPS.PS_sig, rasterBpyrPS.PS-rasterBpyrPS.PS_sig, facecolor='orange', alpha=0.5)

plt.tick_params(axis="y", labelsize=15)
plt.fill_between(rasterBpyrPS2.freqs, rasterBpyrPS2.PS+rasterBpyrPS2.PS_sig, rasterBpyrPS2.PS-rasterBpyrPS2.PS_sig, facecolor='red', alpha=0.5)

plt.tick_params(axis="y", labelsize=15)
plt.fill_between(rasterBpyrPS3.freqs, rasterBpyrPS3.PS+rasterBpyrPS3.PS_sig, rasterBpyrPS3.PS-rasterBpyrPS3.PS_sig, facecolor='blue', alpha=0.5)

plt.tick_params(axis="y", labelsize=15)
plt.fill_between(rasterBpyrPS4.freqs, rasterBpyrPS4.PS+rasterBpyrPS4.PS_sig, rasterBpyrPS4.PS-rasterBpyrPS4.PS_sig, facecolor='green', alpha=0.5)


#plt.savefig("figure7/figure7A.png", dpi=my_dpi)
#plt.show()

my_dpi = 96
plt.figure(figsize=(1200/my_dpi, 500/my_dpi), dpi=my_dpi)
times_f = np.linspace(x[0], x[-1], num=rasterBpyrPS4.tf.shape[1])
xv, yv = np.meshgrid(times_f, rasterBpyrPS4.freqs)
#peak_freqs_idx = np.argmax(rasterPS.tf, axis=0)
#peak_freqs_S = rasterPS.freqs[peak_freqs_idx]
#plt.plot(5.5*np.ones((len(times_f),)), np.linspace(rasterBpyrPS.freqs[0], rasterBpyrPS.freqs[-1], num=len(times_f)), 'k', '--',  linewidth=2)
plt.ylim(0.5, 10)
plt.xlim(0.5, 3)
tf_plot = plt.contourf(xv, yv, rasterBpyrPS4.tf/np.max(rasterBpyrPS4.tf))
plt.xlabel('Time (s)', fontsize=25)
plt.tick_params(axis="x", labelsize=15)
plt.ylabel('Frequency (Hz)', fontsize=25)
plt.tick_params(axis="y", labelsize=15)
plt.colorbar(tf_plot)

#plt.plot(times_f, peak_freqs_S, 'k',  linewidth=1)
plt.yticks([0, 20, 40, 60, 80], [0, 20, 40, 60, 80])
plt.tick_params(axis="y", labelsize=15)
plt.gca().spines['top'].set_visible(False)
plt.gca().spines['right'].set_visible(False)
plt.gca().spines['left'].set_visible(False)
plt.gca().spines['bottom'].set_visible(False)
plt.tight_layout()
'''

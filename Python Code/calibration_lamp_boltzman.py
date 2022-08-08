#%%
from cmath import log
from NIST_Database_webscraping.NIST_Lines import Lines_LIBS
import xarray as xr 
from sifparser.sifparser import SifParser, FindPeaks
import numpy as np
from scipy.signal import find_peaks
import matplotlib.pyplot as plt
import scipy.constants as constant
from scipy.constants import lambda2nu
import pandas as pd
from adjustText import adjust_text
from matplotlib.text import Annotation
from scipy.optimize import curve_fit


c = constant.speed_of_light

def signaltonoise(a, axis=0, ddof=0):
    a = np.asanyarray(a)
    m = a.mean(axis)
    sd = a.std(axis=axis, ddof=ddof)
    return np.where(sd == 0, 0, m/sd)



D = xr.open_dataarray('/Users/briansquires/Documents/LIBS/data/20220802/calibration.nc')
title = 'HgAr Lamp'
# linesAr = Lines_LIBS('Ar', 200, 600, strongLines=True)
linesHg = Lines_LIBS('Hg' ,200, 600, strongLines=True)
# lines = pd.concat([linesHg.data_frame, linesAr.data_frame])
lines = linesHg.data_frame



prominence = 100
peaks, props = find_peaks(D.values.squeeze(), prominence=prominence, width=5)
while len(peaks)>len(lines):
    prominence +=1
    peaks, props = find_peaks(D.values.squeeze(), prominence=prominence, width=5)

wp = [D.Wavelength[i].values for i in peaks]
fig, ax = plt.subplots(figsize=(12,8))
D.plot(ax=ax)
ax.set_title(f'{title}')

import math
def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

log_wI_gA = []
E_k = []
PWL = []
conf=False
PEAKWL = []
for row in lines.iterrows():
    
    wl = row[1]['obs_wl_air(nm)']
    peakwl = find_nearest(wp, wl)
    intensity = D.sel(Wavelength = peakwl, method='nearest').item()

    if np.abs(peakwl-wl)>.2 or intensity<D.max()*.01:
        lines = lines.drop(lines[lines['obs_wl_air(nm)']==wl].index)
        continue
    else:
        PEAKWL.append(peakwl)
    
    l = np.log((10**7/wl)*intensity/row[1]['gA(s^-1)'])
    log_wI_gA.append(l)
    E_k.append(row[1]['Ek(cm-1)'])
    if peakwl in PWL:
        color = 'red'
    else:
        color = 'k'
    if peakwl not in PWL:
        element = row[1]['element']
        sp_num = row[1]['sp_num']
        conf_k = row[1]['conf_k']
        conf_i = row[1]['conf_i']
        offset = D.max()*.1
        if conf == True:
            ax.annotate(f'          {element} {sp_num} \n '+ f'{conf_k}'+ r'$\rightarrow$' + f'{conf_i}' , xy = (peakwl , intensity),
                    xytext=(wl , intensity + 500), color=color,
                    arrowprops=dict(arrowstyle = '->', connectionstyle = 'arc3',facecolor='red'))
        else:
            ax.annotate(f'{element} {sp_num}' , xy = (peakwl , intensity),
                    xytext=(wl , intensity + offset), color = color, 
                    arrowprops=dict(arrowstyle = '->', connectionstyle = 'arc3',facecolor='red'))
            
    PWL.append(peakwl)

DIFF = []
for x,y in zip(PEAKWL,lines['obs_wl_air(nm)']):
    diff = (x - y)
    if np.abs(diff)>.5:
        lines = lines.drop(lines[lines['obs_wl_air(nm)']==y].index)
    else:
        DIFF.append(diff)
lines['obs_wl_air(nm)'] - DIFF
ax2 = ax.twinx()
ax2.scatter(PEAKWL, DIFF, color='red')



k = constant.Boltzmann
print(k)
fig, ax = plt.subplots()
ax.scatter(E_k,log_wI_gA)
ax.set(xlabel = '$E_k (cm^{-1})$',ylabel='$log(\lambda_{ik} I)/g_{ik}A_{ik}$')
fit, params = curve_fit(lambda x, m, b: m*x +b, E_k, log_wI_gA)
e_k = np.linspace(np.min(E_k), np.max(E_k))
ax.plot(e_k, fit[0]*e_k+fit[1], color='red')
ax.annotate(f'T = {-1/0.6950356*fit[0]:.2E}K', xy=(45000,4))


# %%

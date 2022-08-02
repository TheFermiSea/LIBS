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

c = constant.speed_of_light



D = xr.open_dataarray(r'C:\Users\rmb0155\Desktop\Python Code\LIBS\data\20220802\Ni.nc')

from NIST_Database_webscraping.NIST_Lines import Lines_LIBS

lineNi = Lines_LIBS('Ni', D.Wavelength.min().values,D.Wavelength.max().values, strongLines=True)

prominence = 100
peaks, props = find_peaks(D.sel(Delay=1350), prominence=prominence, width=5)
while len(peaks)>len(lineNi.data_frame):
    prominence +=1
    peaks, props = find_peaks(D.sel(Delay=1350), prominence=prominence, width=5)

wp = [D.Wavelength[i].values for i in peaks]
fig, ax = plt.subplots(figsize=(12,8))
D.sel(Delay=1350).plot(ax=ax)

import math
def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

log_wI_gA = []
E_k = []
for row in lineNi.data_frame.iterrows():
    
    wl = row[1]['obs_wl_air(nm)']
    peakwl = find_nearest(wp, wl)
    
    
    intensity = D.sel(Wavelength = peakwl, Delay=1350, method='nearest')
    
    l = np.log(lambda2nu(wl)*intensity.item()/row[1]['gA(s^-1)'])
    log_wI_gA.append(l)
    E_k.append(row[1]['Ek(cm-1)'])

    ax.annotate('Ni '+ str(row[1]['sp_num']) , xy = (wl , intensity),
            xytext=(wl , intensity),
            arrowprops=dict(arrowstyle = '-', connectionstyle = 'arc3',facecolor='red'))
       

plt.scatter(wp, [D.sel(Delay=1350, Wavelength=i) for i in wp], marker='x', color='g')


# %%

# %%

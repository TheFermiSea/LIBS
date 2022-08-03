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


c = constant.speed_of_light



D = xr.open_dataarray('/Users/briansquires/Documents/LIBS/data/20220801/Al2O3@ZnO.nc')

lineAl = Lines_LIBS('Al', D.Wavelength.min().values,D.Wavelength.max().values, strongLines=True)
lineZn = Lines_LIBS('Zn', D.Wavelength.min().values,D.Wavelength.max().values, strongLines=True)
lines = pd.concat([lineAl.data_frame, lineZn.data_frame])


prominence = 100
peaks, props = find_peaks(D.sel(Delay=1350), prominence=prominence, width=5)
while len(peaks)>len(lines):
    prominence +=1
    peaks, props = find_peaks(D.sel(Delay=1350), prominence=prominence, width=5)

wp = [D.Wavelength[i].values for i in peaks]
fig, ax = plt.subplots(figsize=(12,8))
D.sel(Delay=1350).plot(ax=ax)
ax.set_title(r'$Al_2O_3$ (2%) $ZnO$ (98%) ')

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
for row in lines.iterrows():
    
    wl = row[1]['obs_wl_air(nm)']
    peakwl = find_nearest(wp, wl)
    intensity = D.sel(Wavelength = peakwl, Delay=1350, method='nearest').item()
    
    l = np.log(10000000/wl*intensity/row[1]['gA(s^-1)'])
    log_wI_gA.append(l)
    E_k.append(row[1]['Ek(cm-1)'])
    if peakwl not in PWL:
        element = row[1]['element']
        sp_num = row[1]['sp_num']
        conf_k = row[1]['conf_k']
        conf_i = row[1]['conf_i']
        ax.annotate(f'          {element} {sp_num} \n '+ f'{conf_k}'+ r'$\rightarrow$' + f'{conf_i}' , xy = (wl , intensity),
                xytext=(wl-3 , intensity + 500),
                arrowprops=dict(arrowstyle = '->', connectionstyle = 'arc3',facecolor='red'))
    PWL.append(peakwl)



# %%

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import least_squares

def lorentzian( x, x0, a, gam ):
    return a * gam**2 / ( gam**2 + ( x - x0 )**2)

def multi_lorentz( x, params ):
    off = params[0]
    paramsRest = params[1:]
    assert not ( len( paramsRest ) % 3 )
    return off + sum( [ lorentzian( x, *paramsRest[ i : i+3 ] ) for i in range( 0, len( paramsRest ), 3 ) ] )

def res_multi_lorentz( params, xData, yData ):
    diff = [ multi_lorentz( x, params ) - y for x, y in zip( xData, yData ) ]
    return diff

xData = D.Wavelength.values
yData = D.sel(Delay=1350).values
yData = yData / max(yData)

generalWidth = 1

yDataLoc = yData
startValues = [ max( yData ) ]
bounds=( [0], [max(yData)*2])
counter = 0

while max( yDataLoc ) - min( yDataLoc ) > .1:
    counter += 1
    if counter > len(lines): ### max 20 peak...emergency break to avoid infinite loop
        break
    maxP = np.argmax( yDataLoc )
    maxY = yData[ maxP ]
    x0 = xData[ maxP ]
    x0 = find_nearest(lines['obs_wl_air(nm)'].values, x0)
    startValues += [ x0, maxY - min( yDataLoc ), generalWidth ]
    [bounds[0].append(i) for i in [x0-1, (maxY - min( yDataLoc ))/2, .01]]
    [bounds[1].append(i) for i in [x0+1, (maxY - min( yDataLoc))*2, 1]]
    res_lsq = least_squares( res_multi_lorentz, startValues, args=( xData, yData ), bounds=bounds )
    yDataLoc = [ y - multi_lorentz( x, res_lsq['x'] ) for x,y in zip( xData, yData ) ]
    print(counter)
#%%  
fig, ax = plt.subplots()
D.sel(Delay=1350).plot(ax=ax)
ax.plot(xData, max(D.sel(Delay=1350).values)*multi_lorentz(xData, res_lsq['x']))
intensity = [D.sel(Wavelength=i, Delay=1350,method='nearest') for i in res_lsq['x'][1::3]]
ax.scatter(res_lsq['x'][1::3], intensity)

wl_lsq = sorted(res_lsq['x'][1::3])
lines.sort_values('obs_wl_air(nm)')
lines = lines.reset_index()
for row in lines.iterrows():
    obs_wl = row[1]['obs_wl_air(nm)']
    # wl = find_nearest(wl_lsq, obs_wl)
    wl = lines['obs_wl_air(nm)'][row[0]]
    element = row[1]['element']
    sp_num = row[1]['sp_num']
    ax.annotate(f'          {element} {sp_num} ' , xy = (wl, D.sel(Delay=1350, Wavelength=wl, method='nearest')),
                xytext=(wl -3, D.sel(Delay=1350, Wavelength=wl, method='nearest')+500),
                arrowprops=dict(arrowstyle = '->', connectionstyle = 'arc3',facecolor='red'))
annotations = [child for child in ax.get_children() if isinstance(child, Annotation)]
adjust_text(annotations,
            autoalign=True, only_move ={'points' :'y', 'text': 'y', 'objects': 'y'} , ha='center', va='bottom'
           )

  # %%

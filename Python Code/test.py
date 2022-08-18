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
from matplotlib.text import Annotation


c = constant.speed_of_light



D = xr.open_dataarray(r'C:\Users\rmb0155\Desktop\Python Code\LIBS\data\20220809\graidient4_pos=4.nc')
title = 'Ni/Fe'
# linesAr = Lines_LIBS('Ar', 200, 600, strongLines=True)
wmin = D.Wavelength.min().item()
wmax = D.Wavelength.max().item()

linesFe = Lines_LIBS('Fe' ,wmin, wmax, strongLines=True)
linesNi = Lines_LIBS('Ni', wmin, wmax, strongLines=True)
# lines = pd.concat([linesHg.data_frame, linesAr.data_frame])
lines = pd.concat([linesFe.data_frame, linesNi.data_frame])
lines = lines[lines['gA(s^-1)']>10**8]
lines = lines.sort_values('gA(s^-1)')



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
conf=False
PEAKWL = []
for row in lines.iterrows():
    
    wl = row[1]['obs_wl_air(nm)']
    peakwl = find_nearest(wp, wl)
    idx = wp.index(peakwl)
    if props['prominences'][idx]>5000:
        
        
        intensity = D.sel(Wavelength = peakwl, method='nearest').item()
        
        l = np.log(10000000/wl*intensity/row[1]['gA(s^-1)'])
        log_wI_gA.append(l)
        E_k.append(row[1]['Ek(cm-1)'])
        if peakwl in PEAKWL:
            color = 'red'
        else:
            color = 'k'
        if peakwl not in PEAKWL:
            element = row[1]['element']
            sp_num = row[1]['sp_num']
            conf_k = row[1]['conf_k']
            conf_i = row[1]['conf_i']
            offset = D.max()*.1
            if conf == True:
                ax.annotate(f'          {element} {sp_num} \n '+ f'{conf_k}'+ r'$\rightarrow$' + f'{conf_i}' , xy = (peakwl , intensity),
                        xytext=(wl , intensity + offset), color=color,
                        arrowprops=dict(arrowstyle = '->', connectionstyle = 'arc3',facecolor='red'))
            else:
                ax.annotate(f'{element} {sp_num}' , xy = (peakwl , intensity),
                        xytext=(wl , intensity + offset), color = color, 
                        arrowprops=dict(arrowstyle = '->', connectionstyle = 'arc3',facecolor='red'))
            PEAKWL.append(peakwl)
    
fig.savefig(r'C:\Users\rmb0155\Desktop\Python Code\LIBS\data\20220809\gradient4_lines.png')


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

def multi_lorentz_fit(D, numpeaks):
    

    xData = D.Wavelength.values
    yData = D.values.squeeze()
    scaling = max(yData)
    yData = yData / scaling
    prom =.01
    peaks, props = find_peaks(yData, prominence=prom)
    while len(peaks)>numpeaks:
        prom+=.01
        peaks, props = find_peaks(yData, prominence=prom)

    

    generalWidth = 1

    yDataLoc = yData
    startValues = [ max( yData ) ]
    bounds=( [0], [max(yData)*2])
    counter = 0

    while max( yDataLoc ) - min( yDataLoc ) > .01:
        counter += 1
        if counter >numpeaks: ### max 20 peak...emergency break to avoid infinite loop
            break
        maxP = np.argmax( yDataLoc )
        maxY = yData[ maxP ]
        x0 = xData[ maxP ]
        x0 = find_nearest(D.Wavelength[peaks], x0)
        startValues += [ x0, maxY - min( yDataLoc ), generalWidth ]
        [bounds[0].append(i) for i in [x0-1, (maxY - min( yDataLoc ))/2, .01]]
        [bounds[1].append(i) for i in [x0+1, (maxY - min( yDataLoc))*2, 1]]
        res_lsq = least_squares( res_multi_lorentz, startValues, args=( xData, yData ), bounds=bounds )
        yDataLoc = [ y - multi_lorentz( x, res_lsq['x'] ) for x,y in zip( xData, yData ) ]
    
    fits = res_lsq['x']
    offset = fits[0]
    params = list(zip(fits[1::3], fits[2::3]*scaling, fits[3::3]))
    params.sort(key=lambda i:i[0])
    res_lsq_sorted = [offset, *[i for sub in params for i in sub]]
    
    return res_lsq_sorted
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
  
D = xr.open_dataarray('/Users/briansquires/Documents/LIBS/data/20220802/calibration.nc')
title = 'HgAr Lamp'
# linesAr = Lines_LIBS('Ar', 200, 600, strongLines=True)
lines = Lines_LIBS('Hg' ,200, 600, strongLines=True)
# lines = pd.concat([linesHg.data_frame, linesAr.data_frame])
lines = lines.data_frame



prominence = 100
peaks, props = find_peaks(D.values.squeeze(), prominence=prominence, width=5)
while len(peaks)>len(lines):
    prominence +=1
    peaks, props = find_peaks(D.values.squeeze(), prominence=prominence, width=5)

wp = [D.Wavelength[i].values for i in peaks]

DIFF = []
for x,y in zip(PEAKWL,lines['obs_wl_air(nm)']):
    diff = (x - y)
    if np.abs(diff)>1:
        diff = np.nan
    DIFF.append(diff)

#%%

expda = xr.open_dataarray('/Users/briansquires/Documents/LIBS/data/20220802/Ni.nc')

normda = xr.DataArray(dims=expda.dims, coords=expda.coords)
for i in expda.Delay:
    normda.loc[dict(Delay=i)] = expda.sel(Delay=i)/expda.sel(Delay=i).max()
    
normda.plot()
    
# %%

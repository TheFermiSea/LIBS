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
import math
from tqdm import tqdm
from scipy.optimize import curve_fit

k_0 =  0.6950356 #cm^-1/K

def get_lines(elements):
    lineList = []
    for i in elements:
        line = Lines_LIBS(i ,400, 550, strongLines=True)
        lineList.append(line.data_frame)
    return pd.concat(lineList)
    
def peak_identifier(D:xr.DataArray.__getitem__, lines, plot=False, conf=False):
    Dnorm = D/D.max()
    wmin = D.Wavelength.min().item()
    wmax = D.Wavelength.max().item()
    lines = lines[(lines['obs_wl_air(nm)']>wmin) & (lines['obs_wl_air(nm)']<wmax)]
    lines = lines[lines['gA(s^-1)']>5*10**7]
    lines = lines.sort_values('gA(s^-1)')

    prominence = .001
    peaks, props = find_peaks(Dnorm.values.squeeze(), prominence=prominence, width=5)
    while len(peaks)>len(lines):
        prominence +=.001
        peaks, props = find_peaks(Dnorm.values.squeeze(), prominence=prominence, width=5)

    wp = sorted([D.Wavelength[i].values for i in peaks])
    if plot==True:
        fig, ax = plt.subplots(figsize=(12,8))
        D.plot(ax=ax)
        

    def find_nearest(array, value):
        array = np.asarray(array)
        idx = (np.abs(array - value)).argmin()
        return array[idx]
        
    PEAKWL = []
    
    Found_Peaks = []
    
    for row in tqdm(lines.iterrows(),leave=False):
        wl = row[1]['obs_wl_air(nm)']
        peakwl = find_nearest(wp, wl)
        idx = wp.index(peakwl)
        peakerr=peakwl-wl
        if props['prominences'][idx]>.01 and np.abs(peakerr)<.2:
            if peakwl not in PEAKWL:
                intensity = D.sel(Wavelength = peakwl, method='nearest').item()
                l = np.log(10000000/wl*intensity/row[1]['gA(s^-1)'])
                element = row[1]['element']
                sp_num = row[1]['sp_num']
                conf_k = row[1]['conf_k']
                conf_i = row[1]['conf_i']
                df = pd.DataFrame({'element':element, 
                                'sp_num':sp_num, 
                                'conf_k':conf_k, 
                                'conf_i':conf_i, 
                                'peakwl': peakwl,
                                'peakerr': peakwl-wl,
                                'obs_wl_air(nm)': wl,
                                'peakintens': D.sel(Wavelength=peakwl).item(),
                                'E_k': float(row[1]['Ek(cm-1)']),
                                'log_wI_gA': l},
                                index=[0])
                Found_Peaks.append(df)
                offset = D.max()*.1
                if plot == True:
                    if conf == True:
                        ax.annotate(f'          {element} {sp_num-1} \n '+ f'{conf_k}'+ r'$\rightarrow$' + f'{conf_i}' , xy = (peakwl , intensity),
                                xytext=(wl , intensity + offset), color='k',
                                arrowprops=dict(arrowstyle = '->', connectionstyle = 'arc3',facecolor='red'))
                    else:
                        ax.annotate(f'{element} {sp_num-1}' , xy = (peakwl , intensity),
                                xytext=(wl , intensity + offset), color = 'k', 
                                arrowprops=dict(arrowstyle = '->', connectionstyle = 'arc3',facecolor='red'))
                PEAKWL.append(peakwl)
    Found_Peaks = pd.concat(Found_Peaks, ignore_index = True).drop_duplicates().to_xarray()
    Found_Peaks['E_k'].attrs['units'] = 'cm^-1'
    Found_Peaks['peakwl'].attrs['units'] = 'nm'
    return Found_Peaks

def boltzmann(da, element, plot=False):
    d = da.where(da['element'] == element).dropna('index')
    E_k = d['E_k']
    log_wI_gA = d['log_wI_gA']
    def line(x,m,b):
        return m*x+b
    fit, params = curve_fit(line, E_k, log_wI_gA)  
    temp = -1/(k_0*fit[0])
    if plot==True:
        x = np.linspace(min(E_k), max(E_k),100)
        fig, ax = plt.subplots()
        d.where(d['element']==element).dropna('index').plot.scatter(x='E_k', y='log_wI_gA', label=element, ax=ax)
        ax.plot(x, line(x, *fit), '--')
        ax.plot([], [], ' ', label=f'T = {temp:.2e}K')
        k_0_eV = 8.617333262e-5
        ax.plot([], [], ' ', label=f'T = {temp*k_0_eV:.2e}eV')
        plt.legend()
    return temp
      

#%%

                

lines = get_lines(['Fe', 'Ni'])
#%%
D = xr.open_dataarray('/Users/briansquires/Documents/LIBS/data/20220817/1200gmm/1_clean.nc')

FP = []
for i in tqdm(D.Position):
    FP.append(peak_identifier(D.sel(Position=i), lines, plot=True))
position_da = xr.combine_nested(FP, concat_dim='Position').assign_coords({'Position': D.Position}).dropna(dim='index')
position_da.to_netcdf('/Users/briansquires/Documents/LIBS/data/20220817/1200gmm/position_da_1.nc')

def FeNiRatio(position_da):
    Fe = position_da.where(position_da['element']=='Fe')['peakintens'].sum('index')
    Ni = position_da.where(position_da['element']=='Ni')['peakintens'].sum('index')
    
    

# %%
fp = peak_identifier(D.isel(Position=0), ['Fe', 'Ni'], plot=True)
# %%

    
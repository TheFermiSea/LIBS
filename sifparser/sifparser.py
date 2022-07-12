#%%
from xmlrpc.server import SimpleXMLRPCDispatcher
import sif_reader
from typing import List
import xarray as xr 
from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
from glob import glob
from scipy.optimize import curve_fit
from SimulatedLIBS import simulation


# class SifParser():
    
    # def __init__(self):
    #     pass

def SifParser(filename: str,
            normalize=False) -> xr.DataArray.__getitem__:
    '''
    **Assumes calibration of detector was already preformed and stored in firmware
    Input:      .sif data file as string
    Returns:    xr.DataArray object with modified attributes
    
        opens .sif file in an xarray object and modifies the following:
            -reassign 'width' dimension to use calibration and change to wavelength
            -drop calibration coordinate (replaced by Wavelength)
            -divide intensity by exposure time to normalize to counts/second.  Add additional normalization here
            -assign units: Time->s, Wavelength->nm, Intensity-> counts/s
            -assign filename to attributes for tracking input files
            
    TODO: 
        add error checking
        add code to verify dimensionality matching for concatonation
    '''
    da = sif_reader.xr_open(filename)
    da = (da
            .assign_coords({'Wavelength': ('width', da.calibration.data)})
            .set_index(width='Wavelength')
            .rename({'width':'Wavelength'})
            .drop_vars('calibration')
    )
    if normalize==False:
        da['Wavelength'].attrs['units'] = 'nm'
        da['Time'].attrs['units'] = 's'
        da.name = 'Intensity'
        da.attrs['units'] = 'arb. units'
    elif normalize==True:
        da['Wavelength'].attrs['units'] = 'nm'
        da['Time'].attrs['units'] = 's'
        da.name = 'Intensity'
        da.attrs['units'] = 'counts/s'
        danorm = da/da.ExposureTime
        danorm.attrs = da.attrs
        da = danorm
        del danorm
    return da

def FindPeaks(da:xr.DataArray,
              numpeaks:int,
              width:int=5, 
              time:float=0.0, 
              plot=False,
              title=None,
              filename='',
              sim=False,
              elements=[],
              percentages=[])->None:
    '''
    Peak finding algorithm based on the number of expected peaks.
    Peak pixel values are saved in da.attrs['PeakPixels']
    Peak properties are saved in da.attrs['PeakProps'] as a dictionary
    Peak wavelength values are saved in da.attrs['PeakWavelengths']
    if plot=True, an annotated plot is generated
    '''
    Prominence = 0
    peaks, props = find_peaks(da.sel(Time=time).values[0], 
                              prominence=Prominence, 
                              width=width)
    while len(peaks) > numpeaks:
        Prominence += 1
        peaks, props = find_peaks(da.sel(Time=time).values[0], 
                                  prominence=Prominence, 
                                  width=width)
    da.attrs['PeakPixels'] = peaks
    da.attrs['PeakProps'] = props
    da.attrs['PeakWavelengths'] = da.Wavelength[peaks].values
    if plot==True:
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot(111)
        da.plot(ax=ax)
        da.sel(Time=time, Wavelength=da.Wavelength[peaks]).plot(marker='x', linestyle='')
        for peak in peaks:
            ax.annotate(f'{da.Wavelength[peak].values:.2f}nm', 
                        xy=(da.Wavelength[peak],
                            da.sel(Wavelength=da.Wavelength[peak])))
        ax.set_title(f'{title}')
        if sim==True:
            ax2 = ax.twinx()
            simda = get_sim_data(elements, percentages)
            simda.sel(Wavelength=slice(da.Wavelength.min(),da.Wavelength.max())).plot(ax=ax2,color='red')
        plt.savefig(f'{filename}')
    return da 


def DataSetGenerator(filelist: List[str], 
                    dimension='files',
                    normalize=False,
                    numpeaks:int=1,
                    width:int=5,
                    plot=False)->xr.Dataset.__getitem__:
    '''
    Input:      
            filelist : .sif data file as string
            dimension : dimension along which to concatonate
            
    Returns:    concatonated xr.Dataset object 
        
    TODO 
    add *args for dimension
    '''
    dataset = []
    for file in filelist:
        if normalize==True:
            da = SifParser(file,normalize=True)
        else:
            da = SifParser(file,normalize=False)
        # da = FindPeaks(da, numpeaks=numpeaks, width=width, plot=plot)
        dataset.append(da)
    ds = xr.concat(dataset,dim = 'Position')
    # ds = ds.to_array()
    return ds
    

def get_sim_data(elements : List,
             percentages: List,  
             Te=1.0, 
             Ne=10**17):
    libs = simulation.SimulatedLIBS(Te=Te, Ne=Ne, elements=elements,percentages=percentages,
                                    resolution=500,low_w=200,upper_w=1000,max_ion_charge=3)
    spectrum = libs.get_raw_spectrum()
    spectrum = spectrum.values.T.astype(float)
    da = xr.DataArray(spectrum[1],coords = {'Wavelength':spectrum[0]})
    da.name = 'Intensity'
    da.attrs['units'] = 'arb. units'
    da.Wavelength.attrs['units'] = 'nm'
    return da

def p2nm(da, pix):
    def line(x, m, b):
        return m*x+b 
    pixels = np.linspace(0, len(da.Wavelength), len(da.Wavelength))
    fit, props = curve_fit(line, da.Wavelength,pixels)
    return pix/fit[0]

def merge_spectrum(filelist, name):
    D = []
    for i in filelist:
        da = SifParser(i)
        D.append(da)
    DA = xr.merge(D)
    DA = DA.to_array()
    DA['name'] = name
    return DA
        
#%%
air = glob('/Users/briansquires/Documents/LIBS/data/20220706/ZnO/air_*.sif')
argon = glob('/Users/briansquires/Documents/LIBS/data/20220706/ZnO/argon_*.sif')
airDA = merge_spectrum(air,'ZnO Air')
argonDA = merge_spectrum(argon, 'ZnO Argon')
simda = get_sim_data(['Zn'], [100])
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot()
# ax2 = ax.twinx()
# simda.sel(Wavelength=slice(300,650)).plot(ax=ax2, color = 'red')
airDA.plot(ax=ax, label = 'Air')
argonDA.plot(ax=ax, label = 'Argon')
ax.legend()
ax.set_title('Zn')
fig.savefig('/Users/briansquires/Documents/LIBS/data/20220706/ZnO/ZnO.png')

# %%

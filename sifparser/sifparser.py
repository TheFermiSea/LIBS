#%%
import sif_reader
from typing import List
import xarray as xr 
from scipy.signal import find_peaks
import numpy as np
import matplotlib.pyplot as plt
from glob import glob


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
              plot=False)->None:
    '''
    Peak finding algorithm based on the number of expected peaks.
    Peak pixel values are saved in da.attrs['PeakPixels']
    Peak properties are saved in da.attrs['PeakProps'] as a dictionary
    Peak wavelength values are saved in da.attrs['PeakWavelengths']
    if plot=True, an annotated plot is generated
    '''
    Prominence = 0
    peaks, props = find_peaks(da.sel(Time=time).values[0], prominence=Prominence, width=width)
    while len(peaks) > numpeaks:
        Prominence += 1
        peaks, props = find_peaks(da.sel(Time=time).values[0], prominence=Prominence, width=width)
    da.attrs['PeakPixels'] = peaks
    da.attrs['PeakProps'] = props
    da.attrs['PeakWavelengths'] = da.Wavelength[peaks].values
    if plot==True:
        fig = plt.figure(figsize=(10,8))
        ax = fig.add_subplot(111)
        da.plot(ax=ax)
        da.sel(Time=time, Wavelength=da.Wavelength[peaks]).plot(marker='x', linestyle='')
        for peak in peaks:
            ax.annotate(f'{da.Wavelength[peak].values:.2f}nm', xy=(da.Wavelength[peak],da.sel(Wavelength=da.Wavelength[peak])))
    
    # return da 

def DataSetGenerator(filelist: List[str], 
                    dimension='files',
                    normalize=False,
                    numpeaks:int=10,
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
        da = FindPeaks(da, numpeaks=numpeaks, width=width, plot=plot)
        dataset.append(da)
        ds = xr.concat(dataset, dim='files')
        # ds = xr.merge(dataset, compat='override')
        return ds
    
        

#%%       Testbed

# da = SifParser('/Users/briansquires/Documents/LIBS/data/20211026/1us4095mcp1.sif')
# FindPeaks(da, 6, 1, plot=True)
filelist=glob('../data/20211026/*.sif')
ds = DataSetGenerator(filelist)

    
        
    
# %%


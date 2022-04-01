#%%
import sif_reader
from typing import List
import xarray as xr 
from scipy.signal import find_peaks
import numpy as np

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
        peak finding algo
        add con
        
        
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
        da.attrs['units'] = 'counts/s'
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



def DataSetGenerator(filelist: List[str], 
                     dimensions='files',
                    normalize=False)->xr.Dataset.__getitem__:
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
        da, danorm = SifParser(file)
        dataset.append(da)
    # ds = xr.concat(dataset, dim=dimension)
    ds = xr.merge(dataset, compat='override')
    return ds

def FindPeaks(da, numpeaks, width, time=0):
    Prominence = 0
    peaks, props = find_peaks(da.sel(Time=time).values[0], prominence=Prominence, width=width)
    while len(peaks) > numpeaks:
        Prominence += 1
        peaks, props = find_peaks(da.sel(Time=time).values[0], prominence=Prominence, width=width)
    da.attrs['PeakPixels'] = peaks
    da.attrs['PeakProps'] = props
    da.attrs['PeakWavelengths'] = da.Wavelength[peaks]
    
    return da 
        
    
    



#%%       Testbed

da = SifParser('/Users/briansquires/Documents/LIBS/data/20211026/1us4095mcp1.sif')
da = FindPeaks(da, 2, 1)
    
        
    
# da = SifParser('/Users/briansquires/Downloads/10_26_2021/1us4095mcp1.sif')

# %%

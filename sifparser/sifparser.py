#%%
import sif_reader
from typing import List
import xarray as xr 

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
    )/da.attrs['ExposureTime']
    da['Wavelength'].attrs['units'] = 'nm'
    da['Time'].attrs['units'] = 's'
    da.name = 'Intensity'
    da.attrs['units'] = 'counts/s'
    da.attrs['filename'] = filename
    return da

def DataSetGenerator(filelist: List[str], 
                     dimension='files')->xr.Dataset.__getitem__:
    '''
    Input:      
            filelist : .sif data file as string
            dimension : dimension along which to concatonate
            
    Returns:    concatonated xr.Dataset object 
        
    '''
    dataset = []
    for file in filelist:
        da = SifParser(file)
        dataset.append(da)
    ds = xr.concat(dataset, dim=dimension)
    return ds


#%%       Testbed
from glob import glob 

filelist = glob('*.sif')

ds = DataSetGenerator(filelist)
    
        
    
# da = SifParser('/Users/briansquires/Downloads/10_26_2021/1us4095mcp1.sif')

# %%

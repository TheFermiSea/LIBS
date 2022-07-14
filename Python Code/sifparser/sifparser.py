#%%
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
# path = '/Users/briansquires/Documents/LIBS/data/20220713/Argon+Green_IR_test'
# dirs = glob(path + '/*')
# elements = [i.split('/')[-1] for i in dirs]

# for element in elements:
#     filelist = glob(path + '/' + element +'/*.sif')
#     fig = plt.figure(figsize=(12,8))
#     ax = fig.add_subplot()
#     for file in filelist:
#         da = SifParser(file)
#         da.plot(ax=ax, label = file.split('/')[-1].split('.')[0])
#     ax.legend()
#     ax2 = ax.twinx()
#     ax.set_title(element)
#     sim = get_sim_data([element],[100])
#     sim.sel(Wavelength=slice(280,430)).plot(color='black', label='simulation')
#     ax2.legend(loc='upper left')
#     fig.savefig(path + '/' + element + '/plot.png')
# %%
da = SifParser('/Users/briansquires/Documents/LIBS/data/20220713/gradient1.sif')
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot()
lns1 = da.plot(ax = ax, label = 'Experiment')
simfe = get_sim_data(['Fe'],[100])
simni = get_sim_data(['Ni'], [100])
ax2 = ax.twinx()
ax3 = ax.twinx()
lns2 = simfe.sel(Wavelength=slice(278,430)).plot(color='red', ax=ax2, label='Simulation Fe')
lns3 = simni.sel(Wavelength=slice(270,430)).plot(color='green',ax=ax3, label='Simulation Ni')
ax.set_title('Gradient Fe/Ni Sample')
leg = lns1 + lns2 + lns3
labs = [l.get_label() for l in leg]
ax.legend(leg, labs, loc=0)
fig.savefig('/Users/briansquires/Documents/LIBS/data/20220713/gradient.png')
# %%
da = SifParser('/Users/briansquires/Documents/LIBS/data/20220713/calibrationlamp350.sif')

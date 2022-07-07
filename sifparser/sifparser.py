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


from itertools import takewhile, count

list( takewhile( lambda x: x < 10, count(0) ) ) # [ 0, 1, 2, 3, 4, 5, 6, 7, 8, 9 ]


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
        da = FindPeaks(da, numpeaks=numpeaks, width=width, plot=plot)
        dataset.append(da)
    ds = xr.concat(dataset,dim = 'Position')
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
AlFiles = glob('/Users/briansquires/Documents/LIBS/data/20220705/Al/air_*.sif')
DA = merge_spectrum(AlFiles)

#%%       Testbed
filepath = '/Users/briansquires/Documents/LIBS/data/20220705/'
elementpaths = glob(filepath +'*/')
elements = [i.split('/')[-2] for i in elementpaths]
elements.sort()
for element in elements:
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot()
    airfiles = glob(filepath+f'{element}/air_*.sif')
    name = file.split('/')[-1].split('.')[0]
    da = merge_spectrum(airfiles, name)
    plot1 = da.plot(ax=ax, label='Air')
    argonfiles = glob(filepath+f'{element}/argon_*.sif')
    da = merge_spectrum(argonfiles, name)
    plot2 = da.plot(ax=ax, label='Argon')
    ax2 = ax.twinx()
    if element not in ['BN', 'Inconel', 'TiO2']:
        simda = get_sim_data([element], [100])
        plot3 = simda.sel(Wavelength = slice(da.Wavelength.min(), da.Wavelength.max())).plot(ax=ax2, label='Simulation', color='red')
    elif element == 'BN':
        simdaB = get_sim_data(['B'], [100])
        simdaN = get_sim_data(['N'], [100])
        plot3 = simdaB.sel(Wavelength = slice(da.Wavelength.min(), da.Wavelength.max())).plot(ax=ax2, label='Simulation B', color='red')
        plot4 = simdaN.sel(Wavelength = slice(da.Wavelength.min(), da.Wavelength.max())).plot(ax=ax2, label='Simulation N', color='green')
    elif element == 'Inconel':
        simdaNi = get_sim_data(['Ni'], [100])
        simdaCr = get_sim_data(['Cr'], [100])
        plot3 = simdaNi.sel(Wavelength = slice(da.Wavelength.min(), da.Wavelength.max())).plot(ax=ax2, label='Simulation Ni', color='red')
        plot4 = simdaCr.sel(Wavelength = slice(da.Wavelength.min(), da.Wavelength.max())).plot(ax=ax2, label='Simulation Cr', color='green')
    elif element == 'TiO2':
        simdaTi = get_sim_data(['Ti'], [100])
        simdaO = get_sim_data(['O'], [100])
        plot3 = simdaTi.sel(Wavelength = slice(da.Wavelength.min(), da.Wavelength.max())).plot(ax=ax2, label='Simulation Ti', color='red')
        plot4 = simdaO.sel(Wavelength = slice(da.Wavelength.min(), da.Wavelength.max())).plot(ax=ax2, label='Simulation O', color='green')

    ax.set_title(element)
    if 'plot4' in locals():
        plots = plot1+plot2+plot3 + plot4
        del plot4
    else:
        plots = plot1 + plot2 + plot3 
    labs = [l.get_label() for l in plots]
    ax2.legend(plots, labs, loc=0)
    fig.savefig(filepath+f'{element}/{element}.png')
   
    
        


# # filelist=glob('/Users/briansquires/Downloads/20220615/Fe_0-2500ns/*.asc')
# # filelist.sort()
# # delay = np.linspace(0,2500,25)
# # data = np.asarray([np.loadtxt(i).T[1] for i in filelist])
# # wavelength = np.loadtxt(filelist[0]).T[0]
# # fig = plt.figure(figsize=(8,12))
# # ax = fig.add_subplot()
# # da = xr.DataArray(data, coords = {'Delay':delay, 'Wavelength':wavelength})
# # for d in da.Delay:
# #     DA = da.sel(Delay=delay)
# #     DA.plot(ax=ax)
# # %%

# filelist = glob('/Users/briansquires/Documents/LIBS/data/20220701/gradient_sample/*.sif')
# filelist.sort()
# pos = [float(i.split('.sif')[0].split('/')[-1]) for i in filelist]
# DA = DataSetGenerator(filelist, dimension='Position')
# DA = DA.assign_coords({'Position': pos})
# DA.plot()



# %%

#%%
import xarray as xr
import numpy as np
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
from SimulatedLIBS import simulation
from typing import List
from itertools import product
import concurrent

def get_sim_data(elements : List,
             percentages: List,  
             Te=1.0, 
             Ne=10**17,
             resolution=2000):
    libs = simulation.SimulatedLIBS(Te=Te, Ne=Ne, elements=elements,percentages=percentages,
                                    resolution=resolution,low_w=200,upper_w=1000,max_ion_charge=3)
    libs.interpolate(resolution=.1)
    spectrum = libs.get_interpolated_spectrum()
    spectrum = spectrum.values.T.astype(float)
    da = xr.DataArray(spectrum[1],coords = {'Wavelength':spectrum[0]})
    da.name = 'Intensity'
    da.attrs['units'] = 'arb. units'
    da.Wavelength.attrs['units'] = 'nm'
    return da

def get_simulated_data_array(da, elements,percentages, Nes, Tes, Res):
    for i in tqdm(product(da.Ne,da.Te,da.Resolution), total=len(list(product(da.Ne,da.Te,da.Resolution)))):
        ne, te, res = i
        try:
            da.loc[dict(Ne=ne, Te=te, Resolution=res)] = get_sim_data(elements,percentages,Te=te.item(), Ne=ne.item(), resolution=res.item())
        except IndexError:
            continue
        
def populate_data_array(da, elements,percentages, ne, te):
    try:
        da.loc[dict(Ne=ne, Te=te)] = get_sim_data(elements,percentages,Te=te.item(), Ne=ne.item())
    except:
        pass
    
def get_simulated_data_array_mp(da, elements,percentages, Nes, Tes):
    executor = concurrent.futures.ProcessPoolExecutor(10)
    
    static_inputs = [da, elements, percentages]
    iter = [product(da.Ne,da.Te,da.Resolution)] 
    futures = [executor.submit(populate_data_array, static_inputs+list(i)) for i in iter]
    concurrent.futures.wait(futures)

            
#%%


elements = ['Fe']
percentages=[100.0]
Nes = 10*np.logspace(10.0, 18.0, 9, dtype=np.float64)
Tes = np.arange(.3,1,.1)
Res = [5000]

testda = get_sim_data(elements,percentages,Te=Tes[0], Ne=Nes[0], resolution=Res[0])
wavelengths = testda.Wavelength

simda = xr.DataArray(dims=['Ne','Te','Resolution','Wavelength'],
                coords = {'Ne':Nes, 'Te':Tes, 'Resolution':Res,'Wavelength':wavelengths})

get_simulated_data_array(simda, elements, percentages,Nes,Tes, Res)


def mean_square_error(simda, expda, ne, te):
    sim = simda.sel(Ne=ne, Te=te)
    if np.sum(sim) != 0:
        mse = np.sum(np.square((expda/expda.max()) - (sim/sim.max())))
    else:
        mse=np.nan
    RMSE = np.sqrt(mse)
    return np.sum(RMSE)

expda = xr.open_dataarray('/Users/briansquires/Documents/LIBS/data/20220802/Fe.nc')
expda = expda.sel(Delay=1350)
simda = simda.dropna(dim='Ne')

mseda = da = xr.DataArray(dims=['Ne','Te'],
                  coords = {'Ne':simda.Ne, 'Te':simda.Te})

for ne in tqdm(mseda.Ne,leave=False, desc='Ne', position=0):
    for te in tqdm(mseda.Te,leave=False, desc='Te',position=1):
        mseda.loc[dict(Ne=ne, Te=te)] = mean_square_error(simda, expda, ne, te)

dummyda = mseda.stack(z=['Ne','Te'])
best = dummyda.idxmin('z').item()
Nebest, Tebest = best
mseda.attrs['best'] = best

wmin = expda.Wavelength.min()
wmax = expda.Wavelength.max()
#%%
fig, ax = plt.subplots()
simda.sel(Ne=Nebest, Te=Tebest, method='nearest').sel(Wavelength=slice(wmin,wmax)).plot(ax=ax,color='red')
ax2=ax.twinx()
expda.plot(ax=ax2)
ax2.set_title('')
plt.tight_layout()
mseda = mseda.squeeze()
fig.savefig('/Users/briansquires/Documents/LIBS/data/20220808/simulation_regression/Fe.png')

fig, ax = plt.subplots()
mseda.plot(yscale='log', ax=ax)
        
DA = xr.Dataset({'Experimental': expda, 'Simulated': simda, 'MSE':mseda})




# %%

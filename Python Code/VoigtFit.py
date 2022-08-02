#%%
import xarray as xr 
from scipy.optimize import curve_fit
from scipy.special import voigt_profile
from scipy.signal import find_peaks
from NIST_Database_webscraping.NIST_Lines import Lines_LIBS
import matplotlib.pyplot as plt



def pseudo_voigt(x, x0, sigma, gamma, A, c ):
    return A*voigt_profile(x - x0, sigma, gamma) + c



total_da = xr.open_dataarray('/Users/briansquires/Documents/LIBS/data/20220801/Al2O3@ZnO.nc')

da = total_da.sel(Delay=1350)
lineAl = Lines_LIBS('Al', da.Wavelength.min().values,da.Wavelength.max().values, strongLines=True)
lineZn = Lines_LIBS('Zn', da.Wavelength.min().values,da.Wavelength.max().values, strongLines=True)

fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot()
da.plot(ax=ax)
PEAKSAL = []
for index, row in lineAl.data_frame.iterrows():
    expected_wavelength = row['obs_wl_air(nm)']
    sliced_da = da.sel(Wavelength=slice(expected_wavelength-.5, expected_wavelength+.5))
    prominence = 1
    peak, prop = find_peaks(sliced_da.values, prominence = prominence, width = 1)
    while len(peak)>1:
        prominence += .1
        peak, prop = find_peaks(sliced_da.values, prominence = prominence)
    PEAKSAL.append(peak)
    ax.annotate('Al '+ str(row['sp_num']) , xy = (row['obs_wl_air(nm)'] , sliced_da.isel(Wavelength = peak)),
            xytext=(row['obs_wl_air(nm)'] , sliced_da.isel(Wavelength = peak) + 500),
            arrowprops=dict(arrowstyle = '-', connectionstyle = 'arc3',facecolor='red'))
       
PEAKSZN = []
for index, row in lineZn.data_frame.iterrows():
    expected_wavelength = row['obs_wl_air(nm)']
    sliced_da = da.sel(Wavelength=slice(expected_wavelength-.5, expected_wavelength+.5))
    prominence = 1
    peak, prop = find_peaks(sliced_da.values, prominence = prominence, width = 1)
    while len(peak)>1:
        prominence += .1
        peak, prop = find_peaks(sliced_da.values, prominence = prominence)
    PEAKSZN.append(peak)
    ax.annotate('Zn '+ str(row['sp_num'])  , xy = (row['obs_wl_air(nm)'] , sliced_da.isel(Wavelength = peak)),
            xytext=(row['obs_wl_air(nm)'] , sliced_da.isel(Wavelength = peak) + 500),
            arrowprops=dict(arrowstyle = '-', connectionstyle = 'arc3',facecolor='red'))
       



#%%
FIT = []
PARAMS = []
for i in da1.Delay:

    fit, params = curve_fit(pseudo_voigt,
                            da1.Wavelength, 
                                da1.sel(Delay=i), 
                                p0 = [334.4, 10,10, 1000,1], 
                                bounds=(
                                    (334,0,0,0, 0), 
                                    (335,100,100,10000000,100)
                                    ))
    FIT.append(fit)
    PARAMS.append(params)

# %%

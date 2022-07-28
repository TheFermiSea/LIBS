#%%
from NIST_Database_webscraping.NIST_Lines import Lines_LIBS
import xarray as xr 
from sifparser.sifparser import SifParser, FindPeaks
import numpy as np

line = Lines_LIBS('Al',200,940,strongLines=True,first_sp=False)
df = line.data_frame

da = xr.load_dataarray(r'C:\Users\rmb0155\Desktop\Python Code\LIBS\data\20220726\gen2\Al.nc')
da = FindPeaks(da, 9, title = 'Al', plot=True)


PW = da.attrs['PeakWavelengths']
#%%
for peak in PW:
    p = np.round(peak, 2)
    print(df.loc[np.round(df['obs_wl_air(nm)'],0)==p])
    
# %%

#%%
import xarray as xr
import matplotlib.pyplot as plt
from glob import glob

path = '/Users/briansquires/Documents/LIBS/data/20220726'

gen2 = glob('/Users/briansquires/Documents/LIBS/data/20220726/gen2/*.nc')
elements = [i.split('.')[0].split('/')[-1] for i in gen2]

for element in elements:
    g2 = xr.load_dataarray(path + '/gen2/' + element+'.nc')
    g3 = xr.load_dataarray(path + '/gen3/' + element+'.nc')
    fig = plt.figure(figsize=(12,8))
    ax = fig.add_subplot()
    ranges = [(230,300), (300,400), (400,500), (500,600)]
    for range in ranges:
        fig = plt.figure(figsize=(12,8))
        ax = fig.add_subplot()
        g2.sel(Wavelength=slice(*range)).plot(ax=ax, label='Gen 2')
        g3.sel(Wavelength=slice(*range)).plot(ax=ax, label='Gen 3')
        plt.legend(loc='upper right')
        ax.set_title(f'{element}')
        plt.savefig(path+f'/combined/{element}_{range}')
    

    


# %%
# %%
grad1 = xr.load_dataarray('/Users/briansquires/Documents/LIBS/data/20220721/gradient1.nc')
grad2 = xr.load_dataarray('/Users/briansquires/Documents/LIBS/data/20220721/gradient2.nc')

grad1.plot()
plt.savefig('/Users/briansquires/Documents/LIBS/data/20220721/gradient1.png')

grad2.plot()
plt.savefig('/Users/briansquires/Documents/LIBS/data/20220721/gradient2.png')

# %%

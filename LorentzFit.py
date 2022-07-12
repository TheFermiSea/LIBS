#%%
from scipy.optimize import curve_fit
import xarray as xr
from sifparser.sifparser import SifParser, FindPeaks, get_sim_data
import matplotlib.pyplot as plt
import numpy as np
import scipy.constants

def lorentzian(x, x0, a, gam):
    return a * gam**2 / ( gam**2 + ( x - x0 )**2)
def double_lorentzian( x, x01, a1, gam1, x02, a2, gam2, c):
    return (a1 * gam1**2 / ( gam1**2 + ( x - x01 )**2) + a2 * gam2**2 / ( gam2**2 + ( x - x02 )**2)) +c

# def double_voit( x, x01, a1, gam1, x02, a2, gam2, c)
da = SifParser('/Users/briansquires/Documents/LIBS/data/20220701/Al_1800gmm_395nm_2.sif')

x = da.Wavelength.values
y = da.data.squeeze()
fit, params = curve_fit(double_lorentzian, x, y, p0=[394.35, 1e6, 5, 396, 1e6,5, 1000], maxfev=1000000)

x01, a1, gam1, x02, a2, gam2, c = fit
fig = plt.figure(figsize=(12,8))
ax = fig.add_subplot()
ax.plot(x, lorentzian(x, *fit[0:3]), label = 'Al I $P_{1/2}$')
ax.plot(x, lorentzian(x, *fit[3:-1]), label = 'Al I $P_{3/2}$')
ax.plot(x, double_lorentzian(x,*fit), label = 'Fit', color='red')
ax.annotate(f'{fit[0]:.2f}nm \n FWHM = {2*fit[2]:.2f}nm', xy=[fit[0]+1,da.sel(Wavelength=fit[0], method='nearest')/2])
ax.annotate(f'{fit[3]:.2f}nm \n FWHM = {-2*fit[-2]:.2f}nm', xy=[fit[3]-4,da.sel(Wavelength=fit[3], method='nearest')/2])
ax.annotate(r'$ln \left(\frac{\lambda_{m n} I_{m n}}{h c g_{m} A_{m n}}\right)=-\frac{E_{m}}{k T_{\mathrm{e}}}+\ln \left(\frac{N(T)}{U(T)}\right)$', xy=(385,2e6))
ax.annotate(r'$\Delta \lambda_{1 / 2}=2 w\left(\frac{N_{\mathrm{e}}}{10^{16}}\right) + 3.5 A\left(\frac{N_{\mathrm{e}}}{10^{16}}\right)^{1 / 4} \times\left[1-1.2 N_{\mathrm{D}}^{-1 / 3}\right] \times w\left(\frac{N_{\mathrm{e}}}{10^{16}}\right)$', xy=(382,1.75e6))
da.plot(ax=ax, color='black', label='Experiment')
plt.legend()
ax.set_title('Al 1800gmm 395nm')

h = scipy.constants.physical_constants["Planck constant in eV/Hz"]
c = scipy.constants.physical_constants['speed of light in vacuum']
k = scipy.constants.physical_constants['Boltzmann constant in eV/K']
lambda1 = x01*10**-9
E1 = h[0]*c[0]/lambda1
Ak1 = 9.8e7 #s^-1
lambda2 = x02*10**-9
E2 = h[0]*c[0]/lambda2
Ak2 = 4.99e7 #s^-1


log1 = np.log(np.abs(lambda1*a1*np.pi*gam1)/(h[0]*c[0]*2*Ak1))
log2 = np.log(np.abs(lambda2*a2*np.pi*gam2)/(h[0]*c[0]*2*Ak2))

slope = np.abs((log2 - log1)/(E2-E1))
T = 1/(slope*k[0])

ax.annotate(f'T = {T:.2f} K', xy=(402, 1e6))

fig.savefig('LorentzFit_Al.png')




                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                                           # %%
fit
# %%
fit.curvefit_coefficients.values
# %%

#%%
from scipy.optimize import curve_fit
from scipy.special import voigt_profile

def pseudo_voigt(x, x0, sigma, gamma, A, c ):
    return A*voigt_profile(x - x0, sigma, gamma) + c


# %%
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

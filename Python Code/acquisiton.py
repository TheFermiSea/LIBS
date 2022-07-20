#%%
from turtle import color
from pyAndorSDK3 import AndorSDK3
from pyAndorSpectrograph.spectrograph import ATSpectrograph
import numpy as np
import xarray as xr
import pint
import time
from trigger import Trigger
from pyAndorSpectrograph.spectrograph import ATSpectrograph
from sifparser.sifparser import SifParser, merge_spectrum, get_sim_data
from tqdm import tqdm
import matplotlib.pyplot as plt

class LIBS:
    def __init__(self):
        
        self.spc = ATSpectrograph()
        ret = self.spc.Initialize("")
        shm = self.spc.SetWavelength(0, 350)
        self.spc.SetSlitWidth(0,2,100)
        self.Trigger = Trigger()
        self.Trigger.configure()
        self.sdk3 = AndorSDK3()
        self.cam = self.sdk3.GetCamera(0)
        self.cam.SensorCooling = True
        temps = self.cam.available_options_TemperatureControl
        target_temp = temps[0] #if cam.TemperatureControl != temps[0] else temps[1]
        print("Target temperature = {}C".format(target_temp))
        # cam.TemperatureControl = target_temp

        # waiting for temperature to stabilise
        pbar = tqdm(self.cam.TemperatureStatus != "Stabilised", leave=False)
        while pbar:
            desc = f"Temperature: {self.cam.SensorTemperature:.5f}C Status: '{self.cam.TemperatureStatus}'"
            pbar.set_description(desc)
            if self.cam.TemperatureStatus == "Fault":
                err_str = "Camera faulted when cooling to target temperature"
                raise RuntimeError(err_str)
            time.sleep(5)
            if self.cam.TemperatureStatus == "Stabilised":
                break
            
        ret = self.spc.SetNumberPixels(0,self.cam.SensorWidth)
        ret = self.spc.SetPixelWidth(0,self.cam.PixelWidth)
        (ret, self.calibration) = self.spc.GetCalibration(0, self.cam.SensorWidth)
        
        self.cam_params = {
            'TriggerMode' : 'External',
            'GateMode' : 'DDG',
            'DDGIOCEnable' : True,
            'MCPGain' : 3600,
            'DDGOutputDelay' : 1000000,
            'DDGOutputWidth' : 500000,
            'MetadataEnable' : True,
            'SpuriousNoiseFilter' : True,
            'MCPIntelligate' : True
        }
        
        self.spc_params = {
            'Slit Width' : self.spc.GetSlitWidth(0,2)[1],
            'Grating Position' : self.spc.GetWavelength(0)[1],
            'Grove Density' : self.spc.GetGratingInfo(0,self.spc.GetGrating(0)[1],10)[1],
            'Blaze Wavelength' : self.spc.GetGratingInfo(0,self.spc.GetGrating(0)[1],10)[2]
        }
        
        self.set_cam_params(**self.cam_params)
        self.get_background()
        print('System Initialized')
            
    def set_cam_params(self, **params: dict):
        for k in params.items():
            self.cam.__setattr__(k[0],k[1])
            self.update_params()
            
    def update_params(self):
            self.cam_params = {
                'TriggerMode' : self.cam.TriggerMode,
                'GateMode' : self.cam.GateMode,
                'DDGIOCEnable' : self.cam.DDGIOCEnable,
                'MCPGain' : self.cam.MCPGain,
                'DDGOutputDelay' : self.cam.DDGOutputDelay,
                'DDGOutputWidth' : self.cam.DDGOutputWidth,
                'MetadataEnable' : self.cam.MetadataEnable,
                'SpuriousNoiseFilter' : self.cam.SpuriousNoiseFilter,
                'MCPIntelligate' : self.cam.MCPIntelligate
            }
            self.spc_params = {
                'Slit Width' : self.spc.GetSlitWidth(0,1)[1],
                'Grating Position' : self.spc.GetWavelength(0)[1],
                'Grove Density' : self.spc.GetGratingInfo(0,self.spc.GetGrating(0)[1],10)[1],
                'Blaze Wavelength' : self.spc.GetGratingInfo(0,self.spc.GetGrating(0)[1],10)[2]
            }
        
    def img_to_xarray(self, acq):
        da = xr.DataArray(acq.image).sum('dim_0')
        da = da.assign_coords({'Wavelength': ('dim_1', self.calibration)}).set_index(dim_1='Wavelength').rename({'dim_1':'Wavelength'})
        for k  in acq.metadata.__dict__.keys():
            da.attrs[k] = acq.metadata.__dict__[k]
        for k in self.cam_params.keys():
            da.attrs[k] = self.cam_params[k]
        for k in self.spc_params.keys():
            da.attrs[k] = self.spc_params[k]
        da.name = 'Intensity'
        da.Wavelength.attrs['units'] = 'nm'
        da.attrs['units'] = 'arb. units'
        return da
        
    def get_background(self):
        imgsize = self.cam.ImageSizeBytes
        buf = np.empty((imgsize,), dtype='B')
        self.cam.queue(buf, imgsize)
        self.spc.SetShutter(0,1)
        self.cam.TriggerMode = 'Software'
        self.cam.AcquisitionStart()
        self.cam.SoftwareTrigger()
        acq = self.cam.wait_buffer(1000000)
        self.cam.AcquisitionStop()
        self.update_params()
        self.bkg = self.img_to_xarray(acq)
        self.spc.SetShutter(0,0)
        self.cam.flush()
        self.cam.TriggerMode = 'External'

        
    def get_spectrum(self, bkgsub=True):
        if self.spc.AtZeroOrder(0)[1] == 0 :
            imgsize = self.cam.ImageSizeBytes
            buf = np.empty((imgsize,), dtype='B')
            self.cam.queue(buf, imgsize)
            ret = self.spc.SetShutter(0,1)
            self.cam.AcquisitionStart()
            self.Trigger.single_task()
            acq = self.cam.wait_buffer(10000)
            self.cam.AcquisitionStop()
            self.cam.flush()
            self.spc.SetShutter(0,0)
            self.update_params()
            _da  = self.img_to_xarray(acq)
            _attrs = _da.attrs
            _da, self.bkg = xr.align(_da, self.bkg, join='exact')
            da = _da - self.bkg
            da.attrs = _attrs
            if bkgsub == False:
                return _da
            elif bkgsub == True:
                return da 
            else:
                raise TypeError
        else:
            print('Grating at Zero Order!')
  
    def move_grating(self, position):
        ret = self.spc.SetWavelength(0, position)
        (ret, self.calibration) = self.spc.GetCalibration(0, self.cam.SensorWidth)
        self.get_background()
        
    def plot_spectra(self,
                     da : list, 
                     elements : list[str], 
                     percentages : list[float],
                     sim = False) -> xr.DataArray.__getitem__:
        if type(da) is not list:
            fig = plt.figure(figsize=(12,8))
            ax = fig.add_subplot()
            da.plot(ax=ax)
            SIMDA = []
            if sim == True:
                for element in elements:
                    simda = get_sim_data(elements, percentages)
                    ax2 = ax.twinx()
                    simda.sel(Wavelength = slice(da.Wavelength.min(), da.Wavelength.max())).plot(ax=ax2, label = element, color='red')           
            ax.set_title(str(elements))
            return da
        if type(da) is list and len(da)>1:
            DA = xr.merge(da)
            DA = DA.to_array()
            fig = plt.figure(figsize=(12,8))
            ax = fig.add_subplot()
            DA.plot(ax=ax)
            SIMDA = []
            if sim ==True:
                for element in elements:
                    simda = get_sim_data(elements, percentages)
                    ax2 = ax.twinx()
                    simda.sel(Wavelength = slice(DA.Wavelength.min(), DA.Wavelength.max())).plot(ax=ax2, label = element, color='red')           
            ax.set_title(str(elements))
            return DA
            
        
#%%
        
if __name__ == '__main__':
    LIBS = LIBS()
#%% 
DA = []
LIBS.Trigger.single_task()
for i in [300, 400, 500]:
    LIBS.move_grating(i)
    time.sleep(2)
    da = LIBS.get_spectrum()
    if da.max() < 4e9:
        DA.append(da)
    else:
        while da.max()>4e9:
            da = LIBS.get_spectrum()
        DA.append(da)
LIBS.move_grating(350)

DA = LIBS.plot_spectra(da = DA,elements=['Cu'], percentages=[100], sim=True)







# %%
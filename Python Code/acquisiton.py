#%%
from pyAndorSDK3 import AndorSDK3
from pyAndorSpectrograph.spectrograph import ATSpectrograph
import numpy as np
import xarray as xr
import pint
import time
from trigger import Trigger
from pyAndorSpectrograph.spectrograph import ATSpectrograph
from sifparser.sifparser import SifParser
from tqdm import tqdm

class LIBS:
    def __init__(self):
        
        self.spc = ATSpectrograph()
        ret = self.spc.Initialize("")
        shm = self.spc.SetWavelength(0, 350)
        self.spc.SetSlitWidth(0,0,100)
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
            'MCPGain' : 1000,
            'DDGOutputDelay' : 1000000,
            'DDGOutputWidth' : 500000,
            'MetadataEnable' : True,
            'SpuriousNoiseFilter' : True,
            'MCPIntelligate' : True
        }
        
        self.spc_params = {
            'Slit Width' : self.spc.GetSlitWidth(0,0),
            'Grating Position' : self.spc.GetWavelength(0),
            'Grove Density' : self.spc.GetGratingInfo(0)[1],
            'Blaze Wavelength' : self.spc.GetGratingInfo(0)[2]
        }
        
        self.set_cam_params(**self.cam_params)
        self.get_background()
        print('System Initialized')
            
    def set_cam_params(self, **params: dict):
        for k in params.items():
            self.cam.__setattr__(k[0],k[1])
            self.cam_params[str(k[0])] = k[1]
            
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
        self.bkg = self.img_to_xarray(acq)
        self.spc.SetShutter(0,0)
        self.cam.flush()
        self.cam.TriggerMode = 'External'

        
    def get_spectrum(self):
        if self.spc.AtZeroOrder(0)[1] == 0 :
            imgsize = self.cam.ImageSizeBytes
            buf = np.empty((imgsize,), dtype='B')
            self.cam.queue(buf, imgsize)
            self.spc.SetShutter(0,1)
            self.cam.AcquisitionStart()
            self.Trigger.single_task()
            acq = self.cam.wait_buffer(10000)
            self.cam.AcquisitionStop()
            self.cam.flush()
            self.spc.SetShutter(0,0)
            _da  = self.img_to_xarray(acq)
            _attrs = _da.attrs
            da = _da - self.bkg
            da.attrs = _attrs
            return da 
        else:
            print('Grating at Zero Order!')
  
    def move_grating(self, position):
        ret = self.spc.SetWavelength(0, position)
        (ret, self.calibration) = self.spc.GetCalibration(0, self.cam.SensorWidth)
        

        
if __name__ == '__main__':
    LIBS = LIBS()







# %%

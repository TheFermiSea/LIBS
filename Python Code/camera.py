#%%
from pyAndorSDK3 import AndorSDK3
from pyAndorSpectrograph.spectrograph import ATSpectrograph
import numpy as np
import nidaqmx
import xarray as xr
import pint
import time
from trigger import Trigger
from pyAndorSpectrograph.spectrograph import ATSpectrograph
from sifparser.sifparser import SifParser
from tqdm import tqdm, trange

class LIBS:
    def __init__(self):
        
        self.spc = ATSpectrograph()
        ret = self.spc.Initialize("")
        shm = self.spc.SetWavelength(0, 350)
        self.spc.SetSlitWidth(0,1,100)

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
            time.sleep(1)
            if self.cam.TemperatureStatus == "Stabilised":
                break
            
        self.spc.SetNumberPixels(0,self.cam.SensorWidth)
        self.spc.SetPixelWidth(0,self.cam.PixelWidth)
        (ret, self.calibration) = self.spc.GetCalibration(0, self.cam.SensorWidth)
        
        self.set_cam_params()
        self.get_background()
        print('System Initialized')
            
    def set_cam_params(self, **params):
        self.cam.ExposureTime = 0.1
        self.cam.TriggerMode = 'External'
        self.cam.GateMode = 'DDG'
        self.cam.DDGIOCEnable = True
        self.cam.MCPGain = 1000
        self.cam.DDGOutputDelay = 1000000
        self.cam.DDGOutputWidth = 500000
        self.cam.MetadataEnable = True
        self.cam.SpuriousNoiseFilter = True
        self.cam.MCPIntelligate = True
        
        
    def get_background(self):
        imgsize = self.cam.ImageSizeBytes
        buf = np.empty((imgsize,), dtype='B')
        self.cam.queue(buf, imgsize)
        self.spc.SetShutter(0,0)
        acq = self.cam.acquire()
        self.bkg = xr.DataArray(acq.image).sum('dim_1')
        self.spc.SetShutter(0,0)
        self.cam.flush()
        
    def get_spectrum(self):
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

        da  = xr.DataArray(acq.image).sum('dim_1')
        # da = da - self.bkg
        return da   
    
    def img_to_xarray(self, img):
        da = xr.DataArray(img).sum('dim_1')
        da.assign_coords({'Wavelength': self.calibration}) 
        
if __name__ == '__main__':
    LIBS = LIBS()





# %%

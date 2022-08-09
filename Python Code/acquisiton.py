#%%
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
from scipy.signal import find_peaks


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
            'DDGOutputDelay' : 1350000,
            'DDGOutputWidth' : 1000000,
            'MetadataEnable' : True,
            'SpuriousNoiseFilter' : True,
            'MCPIntelligate' : True
        }
        
        self.spc_params = {
            'Slit Width' : self.spc.GetSlitWidth(0,2)[1],
            'Grating Position' : self.spc.GetWavelength(0)[1],
            'Grove Density' : self.spc.GetGratingInfo(0,self.spc.GetGrating(0)[1],10)[1],
            'Blaze Wavelength' : self.spc.GetGratingInfo(0,self.spc.GetGrating(0)[1],10)[2],
            'Output Flipper Mirror' : self.spc.GetFlipperMirror(0,2)[1]
        }
        
        self.set_cam_params(**self.cam_params)
        self.get_background()
        print('System Initialized')
        for k, v in self.cam_params.items():
            print(k, v)
        for k, v in self.spc_params.items():
            print(k,v)
            
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
                'Blaze Wavelength' : self.spc.GetGratingInfo(0,self.spc.GetGrating(0)[1],10)[2],
                'Output Flipper Mirror' : self.spc.GetFlipperMirror(0,2)[1]
            }
        
    def img_to_xarray(self, acq):
        da = xr.DataArray(acq.image).sum('dim_0')
        if self.spc_params['Output Flipper Mirror'] == 0:
            da = da.assign_coords({'Wavelength': ('dim_1', self.calibration)}).set_index(dim_1='Wavelength').rename({'dim_1':'Wavelength'})
        elif self.spc_params['Output Flipper Mirror'] == 1:
            da = da.assign_coords({'Wavelength': ('dim_1', np.flip(self.calibration))}).set_index(dim_1='Wavelength').rename({'dim_1':'Wavelength'})
        for k in acq.metadata.__dict__.keys():
            da.attrs[k] = acq.metadata.__dict__[k]
        for k in self.cam_params.keys():
            da.attrs[k] = self.cam_params[k]
        for k in self.spc_params.keys():
            da.attrs[k] = self.spc_params[k]
        del da.attrs['_config']
        del da.attrs['_atutil']
        da.name = 'Intensity'
        da.Wavelength.attrs['units'] = 'nm'
        da.attrs['units'] = 'arb. units'
        return da
        
    def get_background(self, shutter=True):
        imgsize = self.cam.ImageSizeBytes
        buf = np.empty((imgsize,), dtype='B')
        self.cam.queue(buf, imgsize)
        if shutter:
            self.spc.SetShutter(0,1)
        else:
            pass
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

        
    def get_spectrum(self, bkgsub=True, shutter=True):
        if self.spc.AtZeroOrder(0)[1] == 0 :
            imgsize = self.cam.ImageSizeBytes
            buf = np.empty((imgsize,), dtype='B')
            self.cam.queue(buf, imgsize)
            if shutter==True:
                ret = self.spc.SetShutter(0,1)
            self.cam.AcquisitionStart()
            self.Trigger.single_task()
            acq = self.cam.wait_buffer(10000)
            self.cam.AcquisitionStop()
            self.cam.flush()
            if shutter==True:
                self.spc.SetShutter(0,0)
            self.update_params()
            _da  = self.img_to_xarray(acq)
            _attrs = _da.attrs
            _da, self.bkg = xr.align(_da, self.bkg, join='exact')
            da = _da.astype('float64') - self.bkg.astype('float64')
            da.attrs = _attrs
            if bkgsub == False:
                return _da
            elif bkgsub == True:
                return da 
            else:
                raise TypeError
        else:
            print('Grating at Zero Order!')
            
    def get_calibration_spectrum(self, bkgsub=True):
        if self.spc.AtZeroOrder(0)[1] == 0 :
            self.set_cam_params(**{'MCPGain':10,'TriggerMode': 'Internal', 'GateMode': 'CW On'})
            imgsize = self.cam.ImageSizeBytes
            buf = np.empty((imgsize,), dtype='B')
            self.cam.queue(buf, imgsize)
            ret = self.spc.SetShutter(0,1)
            acq = self.cam.acquire()
            self.cam.flush()
            self.spc.SetShutter(0,0)
            self.update_params()
            _da  = self.img_to_xarray(acq)
            _attrs = _da.attrs
            _da, self.bkg = xr.align(_da, self.bkg, join='exact')
            da = _da.astype('float64') - self.bkg.astype('float64')
            da.attrs = _attrs
            if bkgsub == False:
                return _da
            elif bkgsub == True:
                return da 
            else:
                raise TypeError
        else:
            print('Grating at Zero Order!')
  
    def move_grating(self, position, calibration_collection = False):
        ret = self.spc.SetWavelength(0, position)
        (ret, self.calibration) = self.spc.GetCalibration(0, self.cam.SensorWidth)
        if calibration_collection:
            self.get_background(shutter=False)
        else:
            self.get_background(shutter=True)
        
    def plot_spectra(self,
                     da : list, 
                     elements : list[str], 
                     percentages : list[float],
                     sim = False) -> xr.DataArray.__getitem__:
        from cycler import cycler
        if type(da) is not list:
            fig = plt.figure(figsize=(12,8))
            ax = fig.add_subplot()
            da.plot(ax=ax)
            SIMDA = []
            if sim == True:
                for element in elements:
                    simda = get_sim_data(elements, percentages)
                    ax2 = ax.twinx()
                    color = iter(plt.cm.rainbow(np.linspace(0, 1, 4)))
                    c = next(color)
                    simda.sel(Wavelength = slice(da.Wavelength.min(), da.Wavelength.max())).plot(ax=ax2, label = element, color=c)           
            ax.set_title(str(elements))
            return da
        if type(da) is list and len(da)>1:
            DA = xr.merge(da)
            DA = DA.to_array()
            fig = plt.figure(figsize=(12,8))
            ax = fig.add_subplot()
            lb1 = DA.plot(ax=ax, label='Experiment')
            SIMDA = []
            if sim ==True:
                for element in elements:
                    simda = get_sim_data(element, percentages)
                    ax2 = ax.twinx()
                    color = iter(plt.cm.rainbow(np.linspace(0, 1, 4)))
                    c = next(color)
                    lb2 = simda.sel(Wavelength = slice(DA.Wavelength.min(), DA.Wavelength.max())).plot(ax=ax2, label = element +' Simulation', color=c)     
                    SIMDA.append(lb2)     
            ax.set_title(str(elements))
            lines_labels = [ax.get_legend_handles_labels() for ax in fig.axes]
            lines, labels = [sum(lol, []) for lol in zip(*lines_labels)]
            ax.legend(lines, labels, loc='upper right')
            return DA
            
        
#%%
        
if __name__ == '__main__':
    LIBS = LIBS()
# %%
DA = []
R = np.arange(0,10,.25)
for i in tqdm(R, leave=False):
    # LIBS.set_cam_params(**{'DDGOutputDelay': i*1000})
    input('')
    DA.append(LIBS.get_spectrum())

D = xr.concat(DA, dim='Position')
D = D.assign_coords({'Position': R})
D.name = 'Fe_Ni_Gradient'

#%%
from NIST_Database_webscraping.NIST_Lines import Lines_LIBS

lines = Lines_LIBS('Cu', D.Wavelength.min().values,D.Wavelength.max().values, strongLines=True)

prominence = 10
peaks, props = find_peaks(D.sel(Delay=1350), prominence=prominence)
while len(peaks)>len(lines.data_frame):
    prominence +=1
    peaks, props = find_peaks(D.sel(Delay=1350), prominence=prominence)
#%%
wp = [D.Wavelength[i].values for i in peaks]
fig, ax = plt.subplots(figsize=(12,8))
D.sel(Delay=1350).plot(ax=ax)

import math
def find_nearest(array,value):
    idx = np.searchsorted(array, value, side="left")
    if idx > 0 and (idx == len(array) or math.fabs(value - array[idx-1]) < math.fabs(value - array[idx])):
        return array[idx-1]
    else:
        return array[idx]

for row in lines.data_frame.iterrows():
    wl = row[1]['obs_wl_air(nm)']
    peakwl = find_nearest(wp, wl)
    
    intensity = D.sel(Wavelength = peakwl, Delay=1350, method='nearest')
    ax.annotate('Cu '+ str(row[1]['sp_num']) , xy = (wl , intensity),
            xytext=(wl , intensity),
            arrowprops=dict(arrowstyle = '-', connectionstyle = 'arc3',facecolor='red'))
       
    

plt.scatter(wp, [D.sel(Delay=1350, Wavelength=i) for i in wp], marker='x', color='g')

# plt.scatter(lineNi.data_frame['obs_wl_air(nm)'], [D.sel(Delay=1350, Wavelength=i, method='nearest') for i in lineNi.data_frame['obs_wl_air(nm)']], marker='o', color='r')






# %%
LIBS.set_cam_params(**{'DDGOutputDelay':1300000})
LIBS.move_grating(220)
# %%



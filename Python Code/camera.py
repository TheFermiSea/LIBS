from pyAndorSDK3 import AndorSDK3
from pyAndorSpectrograph.spectrograph import ATSpectrograph
import numpy as np
import nidaqmx
import xarray as xr
import pint
import time
from testbed import Trigger
from pyAndorSpectrograph.spectrograph import ATSpectrograph

spc = ATSpectrograph()
ret = spc.Initialize("")
shm = spc.SetWavelength(0, 350)

Trigger = Trigger()
Trigger.configure()


sdk3 = AndorSDK3()
cam = sdk3.GetCamera(0)
cam.SensorCooling = True
temps = cam.available_options_TemperatureControl
print(temps)
target_temp = temps[0] #if cam.TemperatureControl != temps[0] else temps[1]
print("Target temperature = {}C".format(target_temp))
# cam.TemperatureControl = target_temp

# waiting for temperature to stabilise
while(cam.TemperatureStatus != "Stabilised"):
    time.sleep(5)
    print("Temperature: {:.5f}C".format(cam.SensorTemperature), end="  ")
    print("Status: '{}'".format(cam.TemperatureStatus))
    if cam.TemperatureStatus == "Fault":
        err_str = "Camera faulted when cooling to target temperature"
        raise RuntimeError(err_str)

cam.ExposureTime = 0.1
cam.TriggerMode = 'External'
cam.GateMode = 'DDG'
cam.DDGIOCEnable = True
cam.MCPGain = 1000
cam.DDGOutputDelay = 1000000
cam.DDGOutputWidth = 500000
cam.MetadataEnable = True
SpuriousNoiseFilter = True
MCPIntelligate = True
imgsize = cam.ImageSizeBytes
buf = np.empty((imgsize,), dtype='B')
cam.queue(buf, imgsize)

cam.AcquisitionStart()
Trigger.single_task()
acq = cam.wait_buffer(1000)
cam.AcquisitionStop()
cam.flush()

da = xr.DataArray(acq.image)
da = da.sum('dim_0')
da.plot()


from tokenize import single_quoted
from pyAndorSDK3 import AndorSDK3
from pyAndorSpectrograph.spectrograph import ATSpectrograph
import numpy as np
import nidaqmx
import xarray as xr
import pint
import time
from nidaqmx.constants import (AcquisitionType, CountDirection, Edge,
    READ_ALL_AVAILABLE, TaskMode, TriggerType)
from nidaqmx.stream_readers import CounterReader


class Trigger(object):
    def __init__(self):
        self.task = nidaqmx.Task()
        
    def configure(self, high_time = .1, low_time = .1, samps_per_chan = 1):
        self.task.co_channels.add_co_pulse_chan_time("Dev1/ctr0", high_time=high_time, low_time=low_time)
        self.task.timing.cfg_implicit_timing(
            sample_mode = AcquisitionType.FINITE,
            samps_per_chan = 1)
        
    def single_task(self):
        tic = time.time()
        self.task.start()
        self.task.wait_until_done()
        self.task.stop()
        toc = time.time()
        print(toc - tic)
        
    def close(self):
        self.task.close()
        
# Trigger = Trigger()
# Trigger.configure()
# print('initialized trigger')
# for i in range(10):
#     Trigger.single_task()
# Trigger.close()


# cam.SensorCooling = True

# # getting legal available options for TemperatureControl feature
# temps = cam.available_options_TemperatureControl
# target_temp = temps[0] #if cam.TemperatureControl != temps[0] else temps[1]
# print("Target temperature = {}C".format(target_temp))
# # cam.TemperatureControl = target_temp

# # waiting for temperature to stabilise
# while(cam.TemperatureStatus != "Stabilised"):
#     time.sleep(5)
#     print("Temperature: {:.5f}C".format(cam.SensorTemperature), end="  ")
#     print("Status: '{}'".format(cam.TemperatureStatus))
#     if cam.TemperatureStatus == "Fault":
#         err_str = "Camera faulted when cooling to target temperature"
#         raise RuntimeError(err_str)

# print("Sensor Temperature now Stabilised and Camera is ready to use")

# cam.CycleMode = "Fixed"
# cam.TriggerMode = 'External'
# cam.ExposureTime = 0.1
# print(cam.Delay)

# acq = cam.acquire(timeout=20000)

# acq.show()


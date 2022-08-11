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
            samps_per_chan = samps_per_chan)
        
    def single_task(self):
        self.task.start()
        self.task.wait_until_done()
        self.task.stop()
        
    def close(self):
        self.task.close()
        


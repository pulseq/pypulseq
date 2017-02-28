# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
"""
TODO
@author: Stefan Kroboth
"""

from Sequence import Sequence
import numpy as np
import mr

# Create a new sequence object
seq = Sequence()

# Define FOV and resolution
fov = 220e-3
Nx = 256
Ny = 256

rf = mr.make_sinc_pulse(20.0*np.pi/180.0, duration=4e-3, slice_thickness=5e-3,
                        apodization=0.5, time_bw_product=4)
gz = rf.gz

# Define other gradients and ADC events
deltak = 1/fov
gx = mr.make_trapezoid('x', flat_area=Nx*deltak, flat_time=6.4e-3)
adc = mr.make_adc(Nx, duration=gx.flat_time, delay=gx.rise_time)
gx_pre = mr.make_trapezoid('x', area=-gx.area/2.0, duration=2e-3)
gz_reph = mr.make_trapezoid('z', area=-gz.area/2.0, duration=2e-3)
phase_areas = (np.arange(Ny)-Ny/2.0)*deltak

# Calculate timing
delay_te = 20e-3 - mr.calc_duration(gx_pre) - mr.calc_duration(gz)/2.0 -\
    mr.calc_duration(gx)/2.0
delay_tr = 100e-3 - mr.calc_duration(gx_pre) - mr.calc_duration(gz) -\
    mr.calc_duration(gx) - delay_te

# Loop over phase encodes and define sequence blocks
for i in range(Ny):
    seq.add_block([rf, gz])
    gyPre = mr.make_trapezoid('y', area=phase_areas[i], duration=2e-3)
    seq.add_block([gx_pre, gyPre, gz_reph])
    seq.add_block(mr.make_delay(delay_te))
    seq.add_block([gx, adc])
    seq.add_block(mr.make_delay(delay_tr))

seq.write('gre.seq')

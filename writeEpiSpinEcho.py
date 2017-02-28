# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
"""
@author Stefan Kroboth
"""

from Sequence import Sequence
import numpy as np
import mr

# Create a new sequence object
seq = Sequence()

# Define FOV and resolution
fov = 220e-3
Nx = 64
Ny = 64

# Set system limits
system = mr.opts(max_grad=32, grad_unit='mT/m', max_slew=130,
                 slew_unit='T/m/s', adc_dead_time=10e-6,
                 rf_dead_time=10e-6)

# Create 90 degree slice selection pulse and gradient
rf = mr.make_sinc_pulse(np.pi/2, system, duration=3e-3, slice_thickness=3e-3,
                        apodization=0.5, time_bw_product=4)
gz = rf.gz

# Define other gradients and ADC events
deltak = 1/fov
k_width = Nx*deltak
readout_time = 3.2e-4
gx = mr.make_trapezoid('x', system, flat_area=k_width, flat_time=readout_time)
# Development side note: Due to collections.namedtuple being immutable,
# it is neccessary to create the gx gradient that is going back explicitly.
# In the matlab solution, the amplitude is reversed in the for loop at the
# bottom of the script. In Python we need to witch between the two gradients.
# This has no influence on the produced Pulseq files, meaning that both
# approaches lead to the same result. (SK)
gx_reverse = mr.make_trapezoid('x', system, flat_area=-k_width,
                               flat_time=readout_time)
adc = mr.make_adc(Nx, system, duration=gx.flat_time, delay=gx.rise_time)

# Pre-phasing gradients
pre_time = 8e-4
gx_pre = mr.make_trapezoid('x', system, area=-gx.area/2-deltak/2,
                           duration=pre_time)
gz_reph = mr.make_trapezoid('z', system, area=-gz.area/2, duration=pre_time)
gy_pre = mr.make_trapezoid('y', system, area=-Ny/2*deltak, duration=pre_time)

# Phase blip in shortest possible time
dur = np.ceil(2*np.sqrt(deltak/system['max_slew'])/10e-6)*10e-6
gy = mr.make_trapezoid('y', system, area=deltak, duration=dur)

# Refocusing pulse with spoiling gradients
rf180 = mr.make_block_pulse(np.pi, 500e-6, system)
gz_spoil = mr.make_trapezoid('z', system, area=gz.area*2, duration=3*pre_time)

# Calculate delay time
TE = 55e-3
duration_to_center = (Nx/2+0.5)*mr.calc_duration(gx) +\
    Ny/2*mr.calc_duration(gy)
delay_te1 = TE/2 - mr.calc_duration(gz)/2 - pre_time -\
    mr.calc_duration(gz_spoil) - mr.calc_duration(rf180)/2
delay_te2 = TE/2 - mr.calc_duration(rf180)/2 - mr.calc_duration(gz_spoil) -\
    duration_to_center

# Define sequence blocks
seq.add_block([rf, gz])
seq.add_block([gx_pre, gy_pre, gz_reph])
seq.add_block(mr.make_delay(delay_te1))
seq.add_block(gz_spoil)
seq.add_block(rf180)
seq.add_block(gz_spoil)
seq.add_block(mr.make_delay(delay_te2))
for i in range(Ny):
    if i % 2 == 0:
        seq.add_block([gx, adc])  # read one line of k-space
    else:
        seq.add_block([gx_reverse, adc])  # read one line of k-space (backwds)
    seq.add_block(gy)  # phase blip
seq.add_block(mr.make_delay(1))

seq.write('epi.seq')

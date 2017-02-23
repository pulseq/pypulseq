# -*- coding: utf-8 -*-
"""
Create a gradient echo sequence and export for execution

The |Sequence| class proviedes functionality to create magnetic
resonance sequences (MRI or NMR) from basic building blocks.

This provides an implementation of the open file format for MR
sequences described here: Http://pulseq.github.io/specification.pdf

This example performs the following steps:
    1. Create a slice selective RF pulse for imaging
    2. Create readout gradient and phase encode strategy
    3. Loop through phase encoding and generate sequence blocks
    4. Write the sequence to an open file format suitable for execution
       on a scanner

@author: Stefan Kroboth
"""
# pylint: disable=invalid-name

from Sequence import Sequence
import numpy as np
import matplotlib.pyplot as plt
import mr

##
# Instantiation and gradient limits
# The system gradient limits can be specified in various units _mT/m_,
# _Hz/cm_, or _Hz/m_. However, the limits will be stored internally in
# units of _Hz/m_ for amplitude and _Hz/m/s_ for slew. Unspecified
# hardware parameters will be assigned default values.
system = mr.opts(max_grad=30, grad_unit='mT/m', max_slew=170,
                 slew_unit='T/m/s')

##
# A new sequence object is created by calling the class constructor
seq = Sequence(system)

##
# Sequence events
# Some sequence parameters are defined using standard Python variables
fov = 220e-3
Nx = 64
Ny = 64
slice_thickness = 5e-3

##
# Slice selection
# Key concepts in the sequence description are *blocks* and *events*.
# Blocks describe a group of events that are executed simultanaeously.
# This hierarchical structure means that one event can be used in
# multiple blocks, a common occurrence in MR sequences, particularly
# in imaging sequences.
#
# First, a slice selective RF pulse (and corresponding slice graidnet)
# can be generated using the |makeSincPulse| function.
flip = 15*np.pi/180
rf = mr.make_sinc_pulse(flip, system, duration=4e-3,
                        slice_thickness=slice_thickness,
                        apodization=0.5, time_bw_product=4)
gz = rf.gz

plt.figure  # Dear pylint, this is not pointless at all
plt.plot(1e3*rf.t, rf.signal)
plt.xlabel(r'$time (ms)$')
plt.ylabel(r'$Signal (Hz)$')
plt.show()

##
# Readout gradient
# To define the remaining encoding gradients we need to calculate the
# $k$-space sampling. The Fourier relationship is
#
# $$\Delta k = \frac{1}{FOV}$$
#
# Therefore the area of the readout gradient is $n\Delta k$.
deltak = 1/fov
k_width = Nx*deltak
readout_time = 6.4e-3
gx = mr.make_trapezoid('x', system, flat_area=k_width, flat_time=readout_time)

adc = mr.make_adc(Nx, duration=gx.flat_time, delay=gx.rise_time)

##
# Phase encoding
# To move the $k$-space trajectory away from 0 prior to the readout a
# prephasing gradient must be used. Furthermore rephasing of the slice
# select gradient is required.
gx_pre = mr.make_trapezoid('x', system, area=-gx.area/2, duration=2e-3)
gz_reph = mr.make_trapezoid('z', system, area=-gz.area/2, duration=2e-3)
phase_areas = (np.arange(Ny)-Ny/2)*deltak

##
# Calculate timing
delay_te = 10e-3 - mr.calc_duration(gx_pre) - mr.calc_duration(rf)/2 - \
    mr.calc_duration(gx)/2
delay_tr = 20e-3 - mr.calc_duration(gx_pre) - mr.calc_duration(rf) - \
    mr.calc_duration(gx) - delay_te
delay_1 = mr.make_delay(delay_te)
delay_2 = mr.make_delay(delay_tr)

##
# Define sequence blocks
# Next, the blocks are put together to form the sequence
for i in range(Ny):
    seq.add_block([rf, gz])
    gy_pre = mr.make_trapezoid('y', system, area=phase_areas[i], duration=2e-3)
    seq.add_block([gx_pre, gy_pre, gz_reph])
    seq.add_block(delay_1)
    seq.add_block([gx, adc])
    seq.add_block(delay_2)

##
# Write to file
# The sequence is written to file in compressed form according to the file
# format specification using the /write/ method.
seq.write('bla2.seq')

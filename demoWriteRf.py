# -*- coding: utf-8 -*-
"""
Create a 2D selective RF pulse for a spin echo sequence

This demo defines an entire MRI sequence in Python to selectively excite
a volume. A slice through this excited volume is then imaged with a
slice-selective refocusing pulse.

This example performs the following steps:
    1. Create a 2D RF pulse and corresponding k-space trajectory
    2. Calculate the gradient waveforms and phase encode strategy
    3. Loop through phase encoding and generate sequence blocks.
    4. Write the sequence to an open file format suitable for execution
       on a scanner

@author: Stefan Kroboth
"""
# pylint: disable=invalid-name

from Sequence import Sequence
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import mr

# Create Sequence object
seq = Sequence()

# Sequence parameters are defined using standard Python variables
fov = 220*1e-3           # Field of view
Nx = 256                 # Imaging resolution x
Ny = 256                 # Imaging resolution y
foe = 200*1e-3           # Field of excitation
target_width = 22.5*1e-3  # Diameter of target excitation pattern
n = 8                    # Number of spiral turns
T = 8*1e-3               # Pulse duration

##
# Excitation k-space
# A inward spiral trajectory for the excitation k-space is defined. The
# field-of-excitation and number of spiral turns defines the maximum
# k-space extent.

k_max = (2*n)/foe/2       # Units of 1/m (not rad/m)
tk = np.linspace(0, T-seq.grad_raster_time, T/seq.grad_raster_time)
kx = k_max*(1-tk/T)*np.cos(2*np.pi*n*tk/T)
ky = k_max*(1-tk/T)*np.sin(2*np.pi*n*tk/T)

plt.figure
plt.plot(kx, ky)
plt.xlabel(r'$k_x (1/m)$')
plt.ylabel(r'$k_y (1/m)$')
plt.show()

##
# RF pulse definition
# The RF pulse is defined closely following Pauly et al, JMR 1989;
# 81:43-56. The target excitation is a Gaussian defined by
#
# $$f(x) = a \exp(-|x|^2/\sigma^2)$$
#
# The equivalent in k-space is calculated with the Fourier transform
#
# $$g(x) = b \exp(-\beta^2 |k|^2)$$
#
# where the width is given by
#
# $$\beta = \frac{2\pi K_{\rm max} \sigma}{2\sqrt{2}}$$

tr = np.arange(0, T-seq.rf_raster_time, seq.rf_raster_time)
f = interp1d(tk, kx, kind='linear', fill_value='extrapolate')
kx_rf = f(tr)
f = interp1d(tk, ky, kind='linear', fill_value='extrapolate')
ry_rf = f(tr)
beta = 2*np.pi*k_max*target_width/2/np.sqrt(2)    # Gaussian width in k-space
signal0 = np.exp(-beta**2*(1-tr/T)**2)*np.sqrt((2*np.pi*n*(1-tr/T))**2+1)

##
# Two RF waveforms are superimposed to excite a replica pattern offset 5cm
# in x and y directions. The shifted pattern is achieved with modulation by
# a complex exponential
signal = signal0*(1+np.exp(-1j*2*np.pi*5e-2*(kx_rf+ry_rf)))

plt.figure(2)
plt.plot(1e3*tr, signal.real, 1e3*tr, signal.imag)
plt.xlabel('t (ms)')
plt.ylabel('Signal (Hz)')
plt.show()

##
# Add gradient ramps to achieve the starting gradient value and moment
# (first k-space point) and likewise ramp the gradients to zero afterwards.
# the RF pulse is also padded with zeros during the ramp times.
out = mr.add_ramps([kx, ky], rf=signal)
kx = out[0]
ky = out[1]
signal = out[2]

##
# The gradient waveforms are calculated based from the k-space trajectory
# using the |traj2grad| function, wich internally calculates the finite
# differences. The arbitrary gradient and RF events are then defined using
# functions in the |mr| toolbox
gx = mr.traj2grad(kx)
gy = mr.traj2grad(ky)

rf = mr.make_arbitrary_rf(signal, 20*np.pi/180)
gx_rf = mr.make_arbitrary_grad('x', gx)
gy_rf = mr.make_arbitrary_grad('y', gy)

##
# Define other gradients and ADC events
deltak = 1/fov
gx = mr.make_trapezoid('x', flat_area=Nx*deltak, flat_time=6.4e-3)
adc = mr.make_adc(Nx, duration=gx.flat_time, delay=gx.rise_time)
gx_pre = mr.make_trapezoid('x', area=-gx.area/2, duration=2e-3)
phase_areas = (np.arange(Ny)-Ny/2)*deltak

##
# Refocusing pulse and spoiling gradients
# the refocusing pulse selects a single slice through the excited volume
rf180 = mr.make_block_pulse(np.pi, duration=1e-3, slice_thickness=5e-3)
gz = rf180.gz
gz_spoil = mr.make_trapezoid('z', area=gx.area, duration=2e-3)

##
# Calculate timing
# Echo time and repetition time are TE=20ms, TR=500ms
delayTE1 = (20e-3)/2 - mr.calc_duration(gz_spoil) - mr.calc_duration(rf180)/2
delayTE2 = delayTE1 - mr.calc_duration(gx_pre) - mr.calc_duration(gx)/2
delayTR = 500e-3 - 20e-3 - mr.calc_duration(rf) - mr.calc_duration(gx)/2

##
# Define sequence blocks
# Loop over phase encodes and define sequence blocks
for i in range(Ny):
    seq.add_block([rf, gx_rf, gy_rf])
    seq.add_block(mr.make_delay(delayTE1))
    seq.add_block(gz_spoil)
    seq.add_block([rf180, gz])
    seq.add_block(gz_spoil)
    seq.add_block(mr.make_delay(delayTE2))
    gyPre = mr.make_trapezoid('y', area=phase_areas[i], duration=2e-3)
    seq.add_block([gx_pre, gyPre])
    seq.add_block([gx, adc])
    seq.add_block(mr.make_delay(delayTR))

##
# Write to file
# The sequence is written to file in compressed form according to the file
# format specification using the /write/ method.
seq.write('bla.seq')

# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
"""
TODO
@author: Stefan Kroboth
"""

from Sequence import Sequence
import numpy as np
from scipy.interpolate import interp1d
import mr

# Create a new sequence object
seq = Sequence()

# Define FOV and resolution
fov = 220e-3
nx = 256
ny = 256

# Field of excitation
foe = 200e-3

# Diameter of target excitation pattern
target_width = 22.5e-3

# Number of spiral turns
n = 8

# Pulse duration
T = 8e-3

# Define spiral k-space trajectory
k_max = (2.0*n)/foe/2.0  # Units of 1/m (not rad/m)
tk = np.linspace(0, T-seq.grad_raster_time, T/seq.grad_raster_time)
kx = k_max*(1-tk/T)*np.cos(2*np.pi*n*tk/T)
ky = k_max*(1-tk/T)*np.sin(2*np.pi*n*tk/T)

# Define RF pulse
tr = np.arange(0, T-seq.rf_raster_time, seq.rf_raster_time)
f = interp1d(tk, kx, kind='linear', fill_value='extrapolate')
kxRf = f(tr)
f = interp1d(tk, ky, kind='linear', fill_value='extrapolate')
kyRf = f(tr)
beta = 2.0*np.pi*k_max*target_width/2.0/np.sqrt(2)  # Gaussian width in k-space
signal0 = np.exp(-beta**2*(1.0-tr/T)**2)*np.sqrt((2.0*np.pi*n*(1.0-tr/T))**2+1)
signal = signal0*(1.0 + np.exp(-1j*2.0*np.pi*5e-2*(kxRf+kyRf)))

# Add gradient ramps
out = mr.add_ramps([kx, ky], rf=signal)
kx = out[0]
ky = out[1]
signal = out[2]

rf = mr.make_arbitrary_rf(signal, 20.0*np.pi/180)
gxRf = mr.make_arbitrary_grad('x', mr.traj2grad(kx))
gyRf = mr.make_arbitrary_grad('y', mr.traj2grad(ky))

# Define other gradients and ADC events
deltak = 1.0/fov
gx = mr.make_trapezoid('x', flat_area=nx*deltak, flat_time=6.4e-3)
adc = mr.make_adc(nx, duration=gx.flat_time, delay=gx.rise_time)
gx_pre = mr.make_trapezoid('x', area=-gx.area/2, duration=2e-3)
phase_areas = (np.arange(ny)-ny/2)*deltak

# Refocusing pulse and spoiling gradients
rf180 = mr.make_block_pulse(np.pi, duration=1e-3, slice_thickness=5e-3)
gz = rf180.gz
gzSpoil = mr.make_trapezoid('z', area=gx.area, duration=2e-3)

# Calculate timing (TE=20ms, TR=500ms)
delay_TE1 = 20e-3/2.0 - mr.calc_duration(gzSpoil) - mr.calc_duration(rf180)/2.0
delay_TE2 = delay_TE1 - mr.calc_duration(gx_pre) - mr.calc_duration(gx)/2.0
delay_TR = 500e-3 - 20e-3 - mr.calc_duration(rf) - mr.calc_duration(gx)/2.0

# Loop over phase encodes and define sequence blocks
for i in range(ny):
    seq.add_block([rf, gxRf, gyRf])
    seq.add_block(mr.make_delay(delay_TE1))
    seq.add_block(gzSpoil)
    seq.add_block([rf180, gz])
    seq.add_block(gzSpoil)
    seq.add_block(mr.make_delay(delay_TE2))
    gyPre = mr.make_trapezoid('y', area=phase_areas[i], duration=2e-3)
    seq.add_block([gx_pre, gyPre])
    seq.add_block([gx, adc])
    seq.add_block(mr.make_delay(delay_TR))

seq.set_definition('Scan_ID', 2068)
seq.set_definition('Recon_Mode', 1)
seq.write('external.seq')

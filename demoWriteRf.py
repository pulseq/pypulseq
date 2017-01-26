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

from Sequence import Sequence
import numpy as np
import matplotlib.pyplot as plt
from scipy.interpolate import interp1d
import mr

# Create Sequence object
seq = Sequence()

# Sequence parameters are defined using standard Python variables
fov = 220*1e-3          # Field of view
Nx = 256                 # Imaging resolution x
Ny = 256                 # Imaging resolution y
foe = 200*1e-3          # Field of excitation
targetWidth = 22.5*1e-3 # Diameter of target excitation pattern
n = 8                    # Number of spiral turns
T = 8*1e-3              # Pulse duration

## Excitation k-space
# A inward spiral trajectory for the excitation k-space is defined. The
# field-of-excitation and number of spiral turns defines the maximum
# k-space extent.

kMax = (2*n)/foe/2       # Units of 1/m (not rad/m)
#tk = np.arange(0,T-seq.gradRasterTime, seq.gradRasterTime) # arange is not reliable for floats!
tk = np.linspace(0,T-seq.gradRasterTime, T/seq.gradRasterTime)
kx = kMax*(1-tk/T)*np.cos(2*np.pi*n*tk/T)
ky = kMax*(1-tk/T)*np.sin(2*np.pi*n*tk/T)

plt.figure
plt.plot(kx, ky)
plt.xlabel(r'$k_x (1/m)$')
plt.ylabel(r'$k_y (1/m)$')
plt.show()

## RF pulse definition
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

tr = np.arange(0,T-seq.rfRasterTime, seq.rfRasterTime)
f = interp1d(tk, kx, kind='linear', fill_value='extrapolate')
kxRf = f(tr)
f = interp1d(tk, ky, kind='linear', fill_value='extrapolate')
kyRf = f(tr)
beta = 2*np.pi*kMax*targetWidth/2/np.sqrt(2)    # Gaussian width in k-space
signal0 = np.exp(-beta**2*(1-tr/T)**2)*np.sqrt((2*np.pi*n*(1-tr/T))**2+1)

## 
# Two RF waveforms are superimposed to excite a replica pattern offset 5cm
# in x and y directions. The shifted pattern is achieved with modulation by
# a complex exponential
signal = signal0*(1+np.exp(-1j*2*np.pi*5e-2*(kxRf+kyRf)))

plt.figure(2)
plt.plot(1e3*tr, signal.real, 1e3*tr, signal.imag)
plt.xlabel('t (ms)')
plt.ylabel('Signal (Hz)')
plt.show()

##
# Add gradient ramps to achieve the starting gradient value and moment 
# (first k-space point) and likewise ramp the gradients to zero afterwards.
# the RF pulse is also padded with zeros during the ramp times.
out = mr.addRamps([kx, ky], rf=signal)
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

rf = mr.makeArbitraryRf(signal, 20*np.pi/180)
gxRf = mr.makeArbitraryGrad('x', gx)
gyRf = mr.makeArbitraryGrad('y', gy)

## 
# Define other gradients and ADC events
deltak = 1/fov
gx = mr.makeTrapezoid('x', flatArea=Nx*deltak, flatTime= 6.4e-3)
adc = mr.makeAdc(Nx, duration=gx.flatTime, delay=gx.riseTime)
gxPre = mr.makeTrapezoid('x', area=-gx.area/2, duration=2e-3)
phaseAreas = (np.arange(Ny)-Ny/2)*deltak

##
# Refocusing pulse and spoiling gradients
# the refocusing pulse selects a single slice through the excited volume
rf180 = mr.makeBlockPulse(np.pi, duration=1e-3, sliceThickness=5e-3)
gz = rf180.gz
gzSpoil = mr.makeTrapezoid('z', area=gx.area, duration=2e-3)

##
# Calculate timing
# Echo time and repetition time are TE=20ms, TR=500ms
delayTE1 = (20e-3)/2 - mr.calcDuration(gzSpoil) - mr.calcDuration(rf180)/2
delayTE2 = delayTE1 - mr.calcDuration(gxPre) - mr.calcDuration(gx)/2
delayTR  = 500e-3 - 20e-3 - mr.calcDuration(rf) - mr.calcDuration(gx)/2

##
# Define sequence blocks
# Loop over phase encodes and define sequence blocks
for i in range(Ny):
    seq.addBlock([rf, gxRf, gyRf])
    seq.addBlock(mr.makeDelay(delayTE1))
    seq.addBlock(gzSpoil)
    seq.addBlock([rf180, gz])
    seq.addBlock(gzSpoil)
    seq.addBlock(mr.makeDelay(delayTE2))
    gyPre = mr.makeTrapezoid('y', area=phaseAreas[i], duration=2e-3)
    seq.addBlock([gxPre, gyPre])
    seq.addBlock([gx, adc])
    seq.addBlock(mr.makeDelay(delayTR))

## Write to file
# The sequence is written to file in compressed form according to the file
# format specification using the /write/ method.
seq.write('bla.seq')

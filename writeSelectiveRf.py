# -*- coding: utf-8 -*-

from Sequence import Sequence
import numpy as np
from scipy.interpolate import interp1d
import mr

# Create a new sequence object
seq = Sequence()

# Define FOV and resolution
fov = 220e-3
Nx = 256
Ny = 256

# Field of excitation
foe = 200e-3

# Diameter of target excitation pattern
targetWidth = 22.5e-3

# Number of spiral turns
n = 8

# Pulse duration
T = 8e-3

# Define spiral k-space trajectory
kMax = (2.0*n)/foe/2.0  # Units of 1/m (not rad/m)
tk = np.linspace(0, T-seq.gradRasterTime, T/seq.gradRasterTime)
kx = kMax*(1-tk/T)*np.cos(2*np.pi*n*tk/T)
ky = kMax*(1-tk/T)*np.sin(2*np.pi*n*tk/T)

# Define RF pulse
tr = np.arange(0, T-seq.rfRasterTime, seq.rfRasterTime)
f = interp1d(tk, kx, kind='linear', fill_value='extrapolate')
kxRf = f(tr)
f = interp1d(tk, ky, kind='linear', fill_value='extrapolate')
kyRf = f(tr)
beta = 2.0*np.pi*kMax*targetWidth/2.0/np.sqrt(2)  # Gaussian width in k-space
signal0 = np.exp(-beta**2*(1.0-tr/T)**2)*np.sqrt((2.0*np.pi*n*(1.0-tr/T))**2+1)
signal = signal0*(1.0 + np.exp(-1j*2.0*np.pi*5e-2*(kxRf+kyRf)))

# Add gradient ramps
out = mr.addRamps([kx, ky], rf=signal)
kx = out[0]
ky = out[1]
signal = out[2]

rf = mr.makeArbitraryRf(signal, 20.0*np.pi/180)
gxRf = mr.makeArbitraryGrad('x', mr.traj2grad(kx))
gyRf = mr.makeArbitraryGrad('y', mr.traj2grad(ky))

# Define other gradients and ADC events
deltak = 1.0/fov
gx = mr.makeTrapezoid('x', flatArea=Nx*deltak, flatTime=6.4e-3)
adc = mr.makeAdc(Nx, duration=gx.flatTime, delay=gx.riseTime)
gxPre = mr.makeTrapezoid('x', area=-gx.area/2, duration=2e-3)
phaseAreas = (np.arange(Ny)-Ny/2)*deltak

# Refocusing pulse and spoiling gradients
rf180 = mr.makeBlockPulse(np.pi, duration=1e-3, sliceThickness=5e-3)
gz = rf180.gz
gzSpoil = mr.makeTrapezoid('z', area=gx.area, duration=2e-3)

# Calculate timing (TE=20ms, TR=500ms)
delayTE1 = 20e-3/2.0 - mr.calcDuration(gzSpoil) - mr.calcDuration(rf180)/2.0
delayTE2 = delayTE1 - mr.calcDuration(gxPre) - mr.calcDuration(gx)/2.0
delayTR = 500e-3 - 20e-3 - mr.calcDuration(rf) - mr.calcDuration(gx)/2.0

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

seq.setDefinition('Scan_ID', 2068)
seq.setDefinition('Recon_Mode', 1)
seq.write('external.seq')

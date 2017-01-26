# -*- coding: utf-8 -*-
"""
TODO
@author: Stefan Kroboth
"""

from Sequence import Sequence
import numpy as np
import matplotlib.pyplot as plt
import mr

# Create a new sequence object
seq = Sequence()

# Define FOV and resolution
fov = 220e-3
Nx = 256
Ny = 256

rf = mr.makeSincPulse(20.0*np.pi/180.0, duration=4e-3, sliceThickness=5e-3, apodization=0.5, timeBwProduct=4)
gz = rf.gz

# Define other gradients and ADC events
deltak = 1/fov
gx = mr.makeTrapezoid('x', flatArea=Nx*deltak, flatTime=6.4e-3)
adc = mr.makeAdc(Nx, duration=gx.flatTime, delay=gx.riseTime)
gxPre = mr.makeTrapezoid('x', area=-gx.area/2.0, duration=2e-3)
gzReph = mr.makeTrapezoid('z', area=-gz.area/2.0, duration=2e-3)
phaseAreas = (np.arange(Ny)-Ny/2.0)*deltak

# Calculate timing
delayTE = 20e-3 - mr.calcDuration(gxPre) - mr.calcDuration(gz)/2.0 - mr.calcDuration(gx)/2.0
delayTR = 100e-3 - mr.calcDuration(gxPre) - mr.calcDuration(gz) - mr.calcDuration(gx) - delayTE

# Loop over phase encodes and define sequence blocks
for i in range(Ny):
    seq.addBlock([rf, gz])
    gyPre = mr.makeTrapezoid('y', area=phaseAreas[i], duration=2e-3)
    seq.addBlock([gxPre, gyPre, gzReph])
    seq.addBlock(mr.makeDelay(delayTE))
    seq.addBlock([gx, adc])
    seq.addBlock(mr.makeDelay(delayTR))

seq.write('gre.seq')

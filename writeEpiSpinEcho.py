# -*- coding: utf-8 -*-
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
system = mr.opts(maxGrad=32, gradUnit='mT/m', maxSlew=130, 
                 slewUnit='T/m/s', adcDeadTime=10e-6, 
                 rfDeadTime=10e-6) 

# Create 90 degree slice selection pulse and gradient
rf = mr.makeSincPulse(np.pi/2, system, duration=3e-3, sliceThickness=3e-3, 
                      apodization=0.5, timeBwProduct=4)
gz = rf.gz

# Define other gradients and ADC events
deltak = 1/fov
kWidth = Nx*deltak
readoutTime = 3.2e-4
gx = mr.makeTrapezoid('x', system, flatArea=kWidth, flatTime=readoutTime)
# Development side note: Due to collections.namedtuple being immutable,
# it is neccessary to create the gx gradient that is going back explicitly.
# In the matlab solution, the amplitude is reversed in the for loop at the
# bottom of the script. In Python we need to witch between the two gradients.
# This has no influence on the produced Pulseq files, meaning that both 
# approaches lead to the same result. (SK)
gxReverse = mr.makeTrapezoid('x', system, flatArea=-kWidth, flatTime=readoutTime)
adc = mr.makeAdc(Nx, system, duration=gx.flatTime, delay=gx.riseTime)

# Pre-phasing gradients
preTime = 8e-4
gxPre  = mr.makeTrapezoid('x', system, area=-gx.area/2-deltak/2, duration=preTime)
gzReph = mr.makeTrapezoid('z', system, area=-gz.area/2,          duration=preTime)
gyPre  = mr.makeTrapezoid('y', system, area=-Ny/2*deltak,        duration=preTime)

# Phase blip in shortest possible time
dur = np.ceil(2*np.sqrt(deltak/system['maxSlew'])/10e-6)*10e-6
gy = mr.makeTrapezoid('y', system, area=deltak, duration=dur)

# Refocusing pulse with spoiling gradients
rf180 = mr.makeBlockPulse(np.pi, 500e-6, system)
gzSpoil = mr.makeTrapezoid('z', system, area=gz.area*2, duration=3*preTime)

# Calculate delay time
TE=55e-3
durationToCenter = (Nx/2+0.5)*mr.calcDuration(gx) + Ny/2*mr.calcDuration(gy)
delayTE1 = TE/2 - mr.calcDuration(gz)/2 - preTime -mr.calcDuration(gzSpoil) - mr.calcDuration(rf180)/2
delayTE2 = TE/2 - mr.calcDuration(rf180)/2 - mr.calcDuration(gzSpoil) - durationToCenter

# Define sequence blocks
seq.addBlock([rf, gz])
seq.addBlock([gxPre, gyPre, gzReph])
seq.addBlock(mr.makeDelay(delayTE1))
seq.addBlock(gzSpoil)
seq.addBlock(rf180)
seq.addBlock(gzSpoil)
seq.addBlock(mr.makeDelay(delayTE2))
for i in range(Ny):
    if i%2==0:
        seq.addBlock([gx, adc]) # read one line of k-space
    else:
        seq.addBlock([gxReverse, adc]) # read one line of k-space (backwards)
    seq.addBlock(gy)             # phase blip
seq.addBlock(mr.makeDelay(1))

seq.write('epi.seq')

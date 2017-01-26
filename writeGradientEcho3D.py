# -*- coding: utf-8 -*-
"""
TODO
@author: Stefan Kroboth
"""

from Sequence import Sequence
import numpy as np
import matplotlib.pyplot as plt
import mr

# Define FOV
fov = np.array((190e-3, 190e-3, 190e-3))

# Define resolution
Nx = 64
Ny = Nx
Nz = Nx

# Define other sequence parameters
Tread = 3.2e-3
Tpre = 3e-3
riseTime = 400e-6
Ndummy = 50

# Define system limits
sys = mr.opts(maxGrad=20, gradUnit='mT/m', riseTime=riseTime, rfDeadTime=10e-6, adcDeadTime=10e-6)

# Create a new sequence object
seq = Sequence(system=sys)

# Define other gradients and ADC events
deltak  = 1.0/fov
gx      = mr.makeTrapezoid('x', system=sys, flatArea=Nx*deltak[0], flatTime=Tread)
gxPre   = mr.makeTrapezoid('x', system=sys, area=-gx.area/2, duration=Tpre)
gxSpoil = mr.makeTrapezoid('x', system=sys, area=gx.area,    duration=Tpre)
rf      = mr.makeBlockPulse(8*np.pi/180, system=sys, duration=0.2e-3) # So far only for timing calculation, will be overwritten later (not necessary, becaus duration is known -- just to keep it more consistent with the Matlab version (SK)
areaY = (np.arange(Ny)-Ny/2.0)*deltak[1]
areaZ = (np.arange(Nz)-Nz/2.0)*deltak[2]

# Calculate timing
TE = 10e-3
TR = 40e-3
delayTE = TE - mr.calcDuration(rf)/2.0 - mr.calcDuration(gxPre) - mr.calcDuration(gx)/2.0
delayTR = TR - mr.calcDuration(rf)     - mr.calcDuration(gxPre) - mr.calcDuration(gx) - mr.calcDuration(gxSpoil) - delayTE

# Drive magnetization to steady state
for iY in range(1,Ndummy+1):
    # RF
    # Create non-selective pulse
    rf = mr.makeBlockPulse(8*np.pi/180, system=sys, duration=0.2e-3, 
            phaseOffset=(np.mod(117.0*(iY**2+iY+2)*np.pi/180,2*np.pi)))
    seq.addBlock(rf)
    # Gradients
    gyPre  = mr.makeTrapezoid('y', system=sys, area= areaY[int(np.floor(Ny/2.0))-1], duration=Tpre)
    gyReph = mr.makeTrapezoid('y', system=sys, area=-areaY[int(np.floor(Ny/2.0))-1], duration=Tpre)
    #flupidu
    seq.addBlock([gxPre, gyPre])
    seq.addBlock(mr.makeDelay(delayTE))
    seq.addBlock(gx)
    seq.addBlock([gyReph, gxSpoil])
    seq.addBlock(mr.makeDelay(delayTR))

# Make trapezoids for inner loop to save computation
gyPre  = []
gyReph = []
rf     = []
adc    = []
for iY in range(1,Ny+1):
    gyPre.append( mr.makeTrapezoid('y', system=sys, area= areaY[iY-1], duration=Tpre))
    gyReph.append(mr.makeTrapezoid('y', system=sys, area=-areaY[iY-1], duration=Tpre))
    phaseOffset = np.mod(117.0*(iY**2+iY+2.0)*np.pi/180.0,2*np.pi)
    rf.append(mr.makeBlockPulse(8*np.pi/180.0, system=sys, duration=0.2e-3, phaseOffset=phaseOffset))
    adc.append(mr.makeAdc(Nx, system=sys, duration=gx.flatTime, delay=gx.riseTime, phaseOffset=phaseOffset))


# Loop over phase encodes and define sequence blocks
for iZ in range(Nz):
    gzPre  = mr.makeTrapezoid('z', system=sys, area= areaZ[iZ], duration=Tpre)
    gzReph = mr.makeTrapezoid('z', system=sys, area=-areaZ[iZ], duration=Tpre)
    for iY in range(Ny):
        # Excitation with RF spoiling
        seq.addBlock(rf[iY])

        # Encoding
        seq.addBlock([gxPre, gyPre[iY], gzPre])
        seq.addBlock(mr.makeDelay(delayTE))
        seq.addBlock([gx, adc[iY]])
        seq.addBlock([gyReph[iY], gzReph, gxSpoil])
        seq.addBlock(mr.makeDelay(delayTR))

# Output in Pulseq format for execution
seq.write('external.seq')

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

from Sequence import Sequence
import numpy as np
import matplotlib.pyplot as plt
import mr

## Instantiation and gradient limits
# The system gradient limits can be specified in various units _mT/m_,
# _Hz/cm_, or _Hz/m_. However, the limits will be stored internally in 
# units of _Hz/m_ for amplitude and _Hz/m/s_ for slew. Unspecified 
# hardware parameters will be assigned default values.
system = mr.opts(maxGrad=30, gradUnit='mT/m', maxSlew=170, slewUnit='T/m/s')

##
# A new sequence object is created by calling the class constructor
seq = Sequence(system)

## Sequence events
# Some sequence parameters are defined using standard Python variables
fov = 220e-3
Nx = 64
Ny = 64
sliceThickness = 5e-3

## Slice selection
# Key concepts in the sequence description are *blocks* and *events*.
# Blocks describe a group of events that are executed simultanaeously. 
# This hierarchical structure means that one event can be used in 
# multiple blocks, a common occurrence in MR sequences, particularly 
# in imaging sequences.
#
# First, a slice selective RF pulse (and corresponding slice graidnet)
# can be generated using the |makeSincPulse| function.
flip = 15*np.pi/180
rf = mr.makeSincPulse(flip, system, duration=4e-3, 
        sliceThickness=sliceThickness, apodization=0.5, timeBwProduct=4)
gz = rf.gz

#plt.figure
#plt.plot(1e3*rf.t, rf.signal)
#plt.xlabel(r'$time (ms)$')
#plt.ylabel(r'$Signal (Hz)$')
#plt.show()

## Readout gradient
# To define the remaining encoding gradients we need to calculate the
# $k$-space sampling. The Fourier relationship is
#
# $$\Delta k = \frac{1}{FOV}$$
#
# Therefore the area of the readout gradient is $n\Delta k$.
deltak = 1/fov
kWidth = Nx*deltak
readoutTime = 6.4e-3
gx = mr.makeTrapezoid('x', system, flatArea=kWidth, flatTime=readoutTime)

adc = mr.makeAdc(Nx, duration=gx.flatTime, delay=gx.riseTime)

## Phase encoding
# To move the $k$-space trajectory away from 0 prior to the readout a 
# prephasing gradient must be used. Furthermore rephasing of the slice
# select gradient is required.
gxPre  = mr.makeTrapezoid('x', system, area=-gx.area/2, duration=2e-3)
gzReph = mr.makeTrapezoid('z', system, area=-gz.area/2, duration=2e-3)
phaseAreas = (np.arange(Ny)-Ny/2)*deltak

## Calculate timing
delayTE = 10e-3 - mr.calcDuration(gxPre) - mr.calcDuration(rf)/2 - mr.calcDuration(gx)/2
delayTR = 20e-3 - mr.calcDuration(gxPre) - mr.calcDuration(rf) - mr.calcDuration(gx) - delayTE
delay1 = mr.makeDelay(delayTE)
delay2 = mr.makeDelay(delayTR)

## Define sequence blocks
# Next, the blocks are put together to form the sequence
for i in range(Ny):
    seq.addBlock([rf, gz])
    gyPre = mr.makeTrapezoid('y', system, area=phaseAreas[i], duration=2e-3)
    seq.addBlock([gxPre, gyPre, gzReph])
    seq.addBlock(delay1)
    seq.addBlock([gx, adc])
    seq.addBlock(delay2)

## Write to file
# The sequence is written to file in compressed form according to the file
# format specification using the /write/ method.
seq.write('bla2.seq')

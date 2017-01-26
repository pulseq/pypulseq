# -*- coding: utf-8 -*-
"""
@author Stefan Kroboth
"""

from Sequence import Sequence
import numpy as np
import mr


# Create new sequence object
seq = Sequence() 

# Define parameters
Nx = 256
Nrep = 1

# Create non-selective pulse
rf = mr.makeBlockPulse(np.pi/2, duration=0.1e-3)

# Define delays and ADC evenets
adc = mr.makeAdc(Nx, duration=3.2e-3)
delayTE = 20e-3
delayTR = 1000e-3
delay1 = mr.makeDelay(delayTE)

# Loop over repetitions and define sequence blocks
for i in range(Nrep):
    seq.addBlock(rf)
    #seq.addBlock(mr.makeDelay(delayTE))
    seq.addBlock(delay1)
    seq.addBlock(adc)
    seq.addBlock(mr.makeDelay(delayTR))

# Write to Pulseq file
seq.write('fid.seq')

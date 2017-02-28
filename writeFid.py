# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
"""
@author Stefan Kroboth
"""

from Sequence import Sequence
import numpy as np
import mr


# Create new sequence object
seq = Sequence()

# Define parameters
nx = 256
nrep = 1

# Create non-selective pulse
rf = mr.make_block_pulse(np.pi/2, duration=0.1e-3)

# Define delays and ADC evenets
adc = mr.make_adc(nx, duration=3.2e-3)
delay_te = 20e-3
delay_tr = 1000e-3
delay1 = mr.make_delay(delay_te)

# Loop over repetitions and define sequence blocks
for i in range(nrep):
    seq.add_block(rf)
    seq.add_block(delay1)
    seq.add_block(adc)
    seq.add_block(mr.make_delay(delay_tr))

# Write to Pulseq file
seq.write('fid.seq')

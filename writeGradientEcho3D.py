# -*- coding: utf-8 -*-
# pylint: disable=invalid-name
"""
TODO
@author: Stefan Kroboth
"""

from Sequence import Sequence
import numpy as np
import mr

# Define FOV
fov = np.array((190e-3, 190e-3, 190e-3))

# Define resolution
nx = 64
ny = nx
nz = nx

# Define other sequence parameters
t_read = 3.2e-3
t_pre = 3e-3
rise_time = 400e-6
n_dummy = 50

# Define system limits
sys = mr.opts(max_grad=20, grad_unit='mT/m', rise_time=rise_time,
              rf_dead_time=10e-6, adc_dead_time=10e-6)

# Create a new sequence object
seq = Sequence(system=sys)

# Define other gradients and ADC events
deltak = 1.0/fov
gx = mr.make_trapezoid('x', system=sys, flat_area=nx*deltak[0],
                       flat_time=t_read)
gx_pre = mr.make_trapezoid('x', system=sys, area=-gx.area/2, duration=t_pre)
gx_spoil = mr.make_trapezoid('x', system=sys, area=gx.area, duration=t_pre)
# So far only for timing calculation, will be overwritten later (not necessary,
# because duration is known -- just to keep it more consistent with the Matlab
# version (SK):
rf = mr.make_block_pulse(8*np.pi/180, system=sys, duration=0.2e-3)
area_y = (np.arange(ny)-ny/2.0)*deltak[1]
area_z = (np.arange(nz)-nz/2.0)*deltak[2]

# Calculate timing
TE = 10e-3
TR = 40e-3
delay_TE = TE - mr.calc_duration(rf)/2.0 - mr.calc_duration(gx_pre) -\
    mr.calc_duration(gx)/2.0
delay_TR = TR - mr.calc_duration(rf) - mr.calc_duration(gx_pre) -\
    mr.calc_duration(gx) - mr.calc_duration(gx_spoil) - delay_TE

# Drive magnetization to steady state
for iY in range(1, n_dummy+1):
    # RF
    # Create non-selective pulse
    rf = mr.make_block_pulse(8*np.pi/180, system=sys, duration=0.2e-3,
                             phase_offset=(np.mod(117.0*(iY**2+iY+2)*np.pi/180,
                                                  2*np.pi)))
    seq.add_block(rf)
    # Gradients
    gy_pre = mr.make_trapezoid('y', system=sys,
                               area=area_y[int(np.floor(ny/2.0))-1],
                               duration=t_pre)
    gy_reph = mr.make_trapezoid('y', system=sys,
                                area=-area_y[int(np.floor(ny/2.0))-1],
                                duration=t_pre)

    seq.add_block([gx_pre, gy_pre])
    seq.add_block(mr.make_delay(delay_TE))
    seq.add_block(gx)
    seq.add_block([gy_reph, gx_spoil])
    seq.add_block(mr.make_delay(delay_TR))

# Make trapezoids for inner loop to save computation
gy_pre = []
gy_reph = []
rf = []
adc = []
for iY in range(1, ny+1):
    gy_pre.append(mr.make_trapezoid('y', system=sys, area=area_y[iY-1],
                                    duration=t_pre))
    gy_reph.append(mr.make_trapezoid('y', system=sys, area=-area_y[iY-1],
                                     duration=t_pre))
    phase_offset = np.mod(117.0*(iY**2+iY+2.0)*np.pi/180.0, 2*np.pi)
    rf.append(mr.make_block_pulse(8*np.pi/180.0, system=sys, duration=0.2e-3,
                                  phase_offset=phase_offset))
    adc.append(mr.make_adc(nx, system=sys, duration=gx.flat_time,
                           delay=gx.rise_time, phase_offset=phase_offset))


# Loop over phase encodes and define sequence blocks
for iZ in range(nz):
    gz_pre = mr.make_trapezoid('z', system=sys, area=area_z[iZ],
                               duration=t_pre)
    gz_reph = mr.make_trapezoid('z', system=sys, area=-area_z[iZ],
                                duration=t_pre)
    for iY in range(ny):
        # Excitation with RF spoiling
        seq.add_block(rf[iY])

        # Encoding
        seq.add_block([gx_pre, gy_pre[iY], gz_pre])
        seq.add_block(mr.make_delay(delay_TE))
        seq.add_block([gx, adc[iY]])
        seq.add_block([gy_reph[iY], gz_reph, gx_spoil])
        seq.add_block(mr.make_delay(delay_TR))

# Output in Pulseq format for execution
seq.write('external.seq')

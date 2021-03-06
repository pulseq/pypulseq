# -*- coding: utf-8 -*-
# pylint: disable=too-many-public-methods
# pylint: disable=invalid-name
# pylint: disable=missing-docstring
"""
Unit tests for the functionality in the mr toolbox

@author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
"""

import unittest as ut
from mr import convert, opts, inside_limits, traj2grad, make_arbitrary_rf
from mr import make_delay
from mr import RFPulse, Gradient, ADC, Delay
from random import random
import numpy as np


class TestConvert(ut.TestCase):
    """Test unit conversion functionality

    Tests every possible conversion defined in the function mr.convert().
    If you add functionality in convert, do not forget to add the corresponding
    testcases here.

    Some conversions may seem pointless, like (num/150.0)*150. This is due to
    potential numerical errors. If we test, we have to reproduce the same
    numerical conditions here as in the convert function.

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def __init__(self, *args, **kwargs):
        super(TestConvert, self).__init__(*args, **kwargs)
        self.num = abs(random())
        self.gamma = 42.57747892e6  # Hz/T

    # Gradient units
    def test_Hz_over_m_to_Hz_over_m(self):
        self.assertEqual(convert(self.num, 'Hz/m', 'Hz/m'), self.num)

    def test_Hz_over_m_to_mT_over_m(self):
        self.assertEqual(convert(self.num, 'Hz/m', 'mT/m'),
                         1e3*self.num/self.gamma)

    def test_Hz_over_m_to_rad_over_ms_over_mm(self):
        self.assertEqual(convert(self.num, 'Hz/m', 'rad/ms/mm'),
                         self.num*2*np.pi*1e-6)

    def test_mT_over_m_to_mT_over_m(self):
        self.assertEqual(convert(self.num, 'mT/m', 'mT/m'), self.num)

    def test_mT_over_m_to_Hz_over_m(self):
        self.assertEqual(convert(self.num, 'mT/m', 'Hz/m'),
                         self.num*1e-3*self.gamma)

    def test_mT_over_m_to_rad_over_ms_over_mm(self):
        self.assertEqual(convert(self.num, 'mT/m', 'rad/ms/mm'),
                         self.num*1e-3*self.gamma*2*np.pi*1e-6)

    def test_rad_over_ms_over_mm_to_rad_over_ms_over_mm(self):
        self.assertEqual(convert(self.num, 'rad/ms/mm', 'rad/ms/mm'), self.num)

    def test_rad_over_ms_over_mm_to_Hz_over_m(self):
        self.assertEqual(convert(self.num, 'rad/ms/mm', 'Hz/m'),
                         self.num*1e6/(2*np.pi))

    def test_rad_over_ms_over_mm_to_mT_over_m(self):
        self.assertEqual(convert(self.num, 'rad/ms/mm', 'mT/m'),
                         1e3*(self.num*1e6/(2*np.pi))/self.gamma)

    # Slew units
    def test_Hz_over_m_over_s_to_Hz_over_m_over_s(self):
        self.assertEqual(convert(self.num, 'Hz/m/s', 'Hz/m/s'), self.num)

    def test_Hz_over_m_over_s_to_mT_over_m_over_ms(self):
        self.assertEqual(convert(self.num, 'Hz/m/s', 'mT/m/ms'),
                         self.num/self.gamma)

    def test_Hz_over_m_over_s_to_T_over_m_over_s(self):
        self.assertEqual(convert(self.num, 'Hz/m/s', 'T/m/s'),
                         self.num/self.gamma)

    def test_Hz_over_m_over_s_to_rad_over_ms_over_mm_over_ms(self):
        self.assertEqual(convert(self.num, 'Hz/m/s', 'rad/ms/mm/ms'),
                         self.num*2*np.pi*1e-9)

    def test_mT_over_m_over_ms_to_Hz_over_m_over_s(self):
        self.assertEqual(convert(self.num, 'mT/m/ms', 'Hz/m/s'),
                         self.num*self.gamma)

    def test_mT_over_m_over_ms_to_mT_over_m_over_ms(self):
        self.assertEqual(convert(self.num, 'mT/m/ms', 'mT/m/ms'), self.num)

    def test_mT_over_m_over_ms_to_T_over_m_over_s(self):
        self.assertEqual(convert(self.num, 'mT/m/ms', 'T/m/s'),
                         (self.num*self.gamma)/self.gamma)

    def test_mT_over_m_over_ms_to_rad_over_ms_over_mm_over_ms(self):
        self.assertEqual(convert(self.num, 'mT/m/ms', 'rad/ms/mm/ms'),
                         self.num*self.gamma*2*np.pi*1e-9)

    def test_T_over_m_over_s_to_Hz_over_m_over_s(self):
        self.assertEqual(convert(self.num, 'T/m/s', 'Hz/m/s'),
                         self.num*self.gamma)

    def test_T_over_m_over_s_to_mT_over_m_over_ms(self):
        self.assertEqual(convert(self.num, 'T/m/s', 'mT/m/ms'),
                         (self.num*self.gamma)/self.gamma)

    def test_T_over_m_over_s_to_T_over_m_over_s(self):
        self.assertEqual(convert(self.num, 'T/m/s', 'T/m/s'),
                         self.num)

    def test_T_over_m_over_s_to_rad_over_ms_over_mm_over_ms(self):
        self.assertEqual(convert(self.num, 'T/m/s', 'rad/ms/mm/ms'),
                         self.num*self.gamma*2*np.pi*1e-9)

    def test_rad_over_ms_over_mm_over_ms_to_Hz_over_m_over_s(self):
        self.assertEqual(convert(self.num, 'rad/ms/mm/ms', 'Hz/m/s'),
                         self.num*1e9/(2*np.pi))

    def test_rad_over_ms_over_mm_over_ms_to_mT_over_m_over_ms(self):
        self.assertEqual(convert(self.num, 'rad/ms/mm/ms', 'mT/m/ms'),
                         self.num*1e9/(2*np.pi)/self.gamma)

    def test_rad_over_ms_over_mm_over_ms_to_T_over_m_over_s(self):
        self.assertEqual(convert(self.num, 'rad/ms/mm/ms', 'T/m/s'),
                         self.num*1e9/(2*np.pi)/self.gamma)

    def test_rad_over_ms_over_mm_over_ms_to_rad_over_ms_over_mm_over_ms(self):
        self.assertEqual(convert(self.num, 'rad/ms/mm/ms', 'rad/ms/mm/ms'),
                         self.num)

    # External Gradients
    def test_1_over_s_to_None(self):
        self.assertEqual(convert(self.num, '1/s', ''), self.num)

    def test_1_over_s_to_1_over_s(self):
        self.assertEqual(convert(self.num, '1/s', '1/s'), self.num)

    def test_1_over_s_to_A(self):
        self.assertEqual(convert(self.num, '1/s', 'A'), 150*self.num)

    def test_1_over_s_to_A_over_s(self):
        self.assertEqual(convert(self.num, '1/s', 'A/s'), 150*self.num)

    def test_None_to_1_over_s(self):
        self.assertEqual(convert(self.num, '', '1/s'), self.num)

    def test_None_to_None(self):
        self.assertEqual(convert(self.num, '', ''), self.num)

    def test_None_to_A(self):
        self.assertEqual(convert(self.num, '', 'A'), 150*self.num)

    def test_None_to_A_over_s(self):
        self.assertEqual(convert(self.num, '', 'A/s'), 150*self.num)

    def test_A_to_1_over_s(self):
        self.assertEqual(convert(self.num, 'A', '1/s'), self.num/150.0)

    def test_A_to_None(self):
        self.assertEqual(convert(self.num, 'A', ''), self.num/150.0)

    def test_A_to_A(self):
        self.assertEqual(convert(self.num, 'A', 'A'), self.num)

    def test_A_to_A_over_s(self):
        self.assertEqual(convert(self.num, 'A', 'A/s'), (self.num/150.0)*150)

    def test_A_over_s_to_1_over_s(self):
        self.assertEqual(convert(self.num, 'A/s', '1/s'), self.num/150.0)

    def test_A_over_s_to_None(self):
        self.assertEqual(convert(self.num, 'A/s', ''), self.num/150.0)

    def test_A_over_s_to_A(self):
        self.assertEqual(convert(self.num, 'A/s', 'A'), (self.num/150.0)*150)

    def test_A_over_s_to_A_over_s(self):
        self.assertEqual(convert(self.num, 'A/s', 'A/s'), self.num)


class TestRamps(ut.TestCase):
    """Test add_ramp(), calc_ramp(),...

    TODO: implement

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def __init__(self, *args, **kwargs):
        super(TestRamps, self).__init__(*args, **kwargs)


class TestOpts(ut.TestCase):
    """Test opt struct creation functionality

    Tests opt() function of mr.py

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def __init__(self, *args, **kwargs):
        super(TestOpts, self).__init__(*args, **kwargs)
        self.num = abs(random())
        self.gamma = 42.57747892e6  # Hz/T

    def test_opts_all_None(self):
        opt = {'grad_unit': 'Hz/m',
               'slew_unit': 'Hz/m/s',
               'max_grad': 40*1e-3*self.gamma,
               'max_slew': 170*self.gamma,
               'rise_time': None,
               'rf_dead_time': 0,
               'rf_ringdown_time': 0,
               'adc_dead_time': 0,
               'rf_raster_time': 1e-6,
               'grad_raster_time': 10e-6}
        self.assertEqual(opts(), opt)

    def test_opts_max_grad_rand(self):
        opt = {'grad_unit': 'Hz/m',
               'slew_unit': 'Hz/m/s',
               'max_grad': self.num*1e-3*self.gamma,
               'max_slew': 170*self.gamma,
               'rise_time': None,
               'rf_dead_time': 0,
               'rf_ringdown_time': 0,
               'adc_dead_time': 0,
               'rf_raster_time': 1e-6,
               'grad_raster_time': 10e-6}
        self.assertEqual(opts(max_grad=self.num), opt)

    def test_opts_max_slew_rand(self):
        opt = {'grad_unit': 'Hz/m',
               'slew_unit': 'Hz/m/s',
               'max_grad': 40*1e-3*self.gamma,
               'max_slew': self.num*self.gamma,
               'rise_time': None,
               'rf_dead_time': 0,
               'rf_ringdown_time': 0,
               'adc_dead_time': 0,
               'rf_raster_time': 1e-6,
               'grad_raster_time': 10e-6}
        self.assertEqual(opts(max_slew=self.num), opt)

    def test_opts_rise_time_rand(self):
        opt = {'grad_unit': 'Hz/m',
               'slew_unit': 'Hz/m/s',
               'max_grad': 40*1e-3*self.gamma,
               'max_slew': None,
               'rise_time': self.num,
               'rf_dead_time': 0,
               'rf_ringdown_time': 0,
               'adc_dead_time': 0,
               'rf_raster_time': 1e-6,
               'grad_raster_time': 10e-6}
        self.assertEqual(opts(rise_time=self.num), opt)

    def test_opts_rf_dead_time_rand(self):
        opt = {'grad_unit': 'Hz/m',
               'slew_unit': 'Hz/m/s',
               'max_grad': 40*1e-3*self.gamma,
               'max_slew': 170*self.gamma,
               'rise_time': None,
               'rf_dead_time': self.num,
               'rf_ringdown_time': 0,
               'adc_dead_time': 0,
               'rf_raster_time': 1e-6,
               'grad_raster_time': 10e-6}
        self.assertEqual(opts(rf_dead_time=self.num), opt)

    def test_opts_rf_ringdown_time_rand(self):
        opt = {'grad_unit': 'Hz/m',
               'slew_unit': 'Hz/m/s',
               'max_grad': 40*1e-3*self.gamma,
               'max_slew': 170*self.gamma,
               'rise_time': None,
               'rf_dead_time': 0,
               'rf_ringdown_time': self.num,
               'adc_dead_time': 0,
               'rf_raster_time': 1e-6,
               'grad_raster_time': 10e-6}
        self.assertEqual(opts(rf_ringdown_time=self.num), opt)

    def test_opts_adc_dead_time_rand(self):
        opt = {'grad_unit': 'Hz/m',
               'slew_unit': 'Hz/m/s',
               'max_grad': 40*1e-3*self.gamma,
               'max_slew': 170*self.gamma,
               'rise_time': None,
               'rf_dead_time': 0,
               'rf_ringdown_time': 0,
               'adc_dead_time': self.num,
               'rf_raster_time': 1e-6,
               'grad_raster_time': 10e-6}
        self.assertEqual(opts(adc_dead_time=self.num), opt)

    def test_opts_rf_raster_time_rand(self):
        opt = {'grad_unit': 'Hz/m',
               'slew_unit': 'Hz/m/s',
               'max_grad': 40*1e-3*self.gamma,
               'max_slew': 170*self.gamma,
               'rise_time': None,
               'rf_dead_time': 0,
               'rf_ringdown_time': 0,
               'adc_dead_time': 0,
               'rf_raster_time': self.num,
               'grad_raster_time': 10e-6}
        self.assertEqual(opts(rf_raster_time=self.num), opt)

    def test_opts_grad_raster_time_rand(self):
        opt = {'grad_unit': 'Hz/m',
               'slew_unit': 'Hz/m/s',
               'max_grad': 40*1e-3*self.gamma,
               'max_slew': 170*self.gamma,
               'rise_time': None,
               'rf_dead_time': 0,
               'rf_ringdown_time': 0,
               'adc_dead_time': 0,
               'rf_raster_time': 1e-6,
               'grad_raster_time': self.num}
        self.assertEqual(opts(grad_raster_time=self.num), opt)

    def test_opts_external_grads(self):
        opt = {'grad_unit': '',
               'slew_unit': '',
               'max_grad': self.num/150.0,
               'max_slew': self.num/150.0,
               'rise_time': None,
               'rf_dead_time': 0,
               'rf_ringdown_time': 0,
               'adc_dead_time': 0,
               'rf_raster_time': 1e-6,
               'grad_raster_time': 10e-6}
        self.assertEqual(opts(slew_unit='A/s', grad_unit='A',
                              max_grad=self.num, max_slew=self.num), opt)


class TestInsideLimits(ut.TestCase):
    """Test inside_limits functionality

    Tests inside_limits() function of mr.py

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def test_grad_and_slew_inside_limits(self):
        self.assertTrue(inside_limits([0, 0.2], [0, -0.3], 1, 1))

    def test_grad_outside_limits(self):
        self.assertFalse(inside_limits([0, 4], [0, 0.2], 1, 7))

    def test_slew_outside_limits(self):
        self.assertFalse(inside_limits([0, 0.4], [6, 0.2], 3, 1))

    def test_grad_and_slew_outside_limits(self):
        self.assertFalse(inside_limits([8, 7.2], [2, -9.3], 4, 3.5))


class TestTraj2Grad(ut.TestCase):
    """Test trajectory to gradient conversion functionality

    Tests traj2grad() function of mr.py
    TODO: more sophisticated tests may be necessary!

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def test_basic_traj(self):
        traj = np.array([[0, 1], [1, 0]], dtype=np.float64)
        grad_raster_time = 10e-6
        self.assertTrue((traj2grad(traj) ==
                         np.array([[1, 0], [-1, 0]])/grad_raster_time).all())

    def test_grad_raster_time_provided(self):
        traj = np.array([[0, 1], [1, 0]], dtype=np.float64)
        grad_raster_time = 1e-6
        self.assertTrue((traj2grad(traj, grad_raster_time=grad_raster_time) ==
                         np.array([[1, 0], [-1, 0]])/grad_raster_time).all())

    def test_sys_provided(self):
        traj = np.array([[0, 1], [1, 0]], dtype=np.float64)
        grad_raster_time = 10e-6
        sys = opts()
        self.assertTrue((traj2grad(traj, system=sys) ==
                         np.array([[1, 0], [-1, 0]])/grad_raster_time).all())

    def test_modified_sys_provided(self):
        traj = np.array([[0, 1], [1, 0]], dtype=np.float64)
        grad_raster_time = 1e-3
        sys = opts(grad_raster_time=1e-3)
        self.assertTrue((traj2grad(traj, system=sys) ==
                         np.array([[1, 0], [-1, 0]])/grad_raster_time).all())

    def test_modified_sys_and_grad_raster_time_provided(self):
        traj = np.array([[0, 1], [1, 0]], dtype=np.float64)
        grad_raster_time = 1e-4
        sys = opts(grad_raster_time=1e-3)
        self.assertTrue((traj2grad(traj, system=sys,
                                   grad_raster_time=grad_raster_time) ==
                         np.array([[1, 0], [-1, 0]])/grad_raster_time).all())


class TestMakeArbitraryRF(ut.TestCase):
    """Test make_arbitrary_rf()

    TODO: more sophisticated tests may be necessary!
    TODO: test with rf_ringdown_time
    TODO: test with bandwith
    TODO: test with slice_thickness
    TODO: test with time_bw_product

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def __init__(self, *args, **kwargs):
        super(TestMakeArbitraryRF, self).__init__(*args, **kwargs)
        self.N = 1024
        self.signal = np.random.rand(self.N)

    def test_basic_freq_and_phase_provided(self):
        system = opts()
        t = np.arange(1, self.N+1)*system['rf_raster_time']
        freq_offset = abs(random())
        phase_offset = abs(random())
        flip_angle = abs(random())
        signal = self.signal / np.sum(self.signal * system['rf_raster_time'])\
            * flip_angle / (2*np.pi)
        gz = None
        rfp = RFPulse('rf', signal, t, freq_offset, phase_offset,
                      system['rf_dead_time'], system['rf_ringdown_time'], gz)
        rf = make_arbitrary_rf(self.signal, flip_angle,
                               freq_offset=freq_offset,
                               phase_offset=phase_offset)
        self.assertTrue((rf.signal == rfp.signal).all())
        self.assertTrue((rf.t == rfp.t).all())
        self.assertTrue(rf.freq_offset == rfp.freq_offset)
        self.assertTrue(rf.phase_offset == rfp.phase_offset)
        self.assertTrue(rf.dead_time == rfp.dead_time)
        self.assertTrue(rf.ringdown_time == rfp.ringdown_time)
        self.assertTrue(rf.gz == rfp.gz)


class TestMakeSincPulse(ut.TestCase):
    """Test make_sinc_pulse()

    TODO: implement

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def __init__(self, *args, **kwargs):
        super(TestMakeSincPulse, self).__init__(*args, **kwargs)


class TestMakeTrapezoid(ut.TestCase):
    """Test make_trapezoid()

    TODO: implement

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def __init__(self, *args, **kwargs):
        super(TestMakeTrapezoid, self).__init__(*args, **kwargs)


class TestMakeArbitraryGrad(ut.TestCase):
    """Test make_arbitrary_grad()

    TODO: implement

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def __init__(self, *args, **kwargs):
        super(TestMakeArbitraryGrad, self).__init__(*args, **kwargs)


class TestMakeADC(ut.TestCase):
    """Test make_adc()

    TODO: implement

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def __init__(self, *args, **kwargs):
        super(TestMakeADC, self).__init__(*args, **kwargs)


class TestMakeBlockPulse(ut.TestCase):
    """Test make_block_pulse()

    TODO: implement

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def __init__(self, *args, **kwargs):
        super(TestMakeBlockPulse, self).__init__(*args, **kwargs)


class TestCalcDuration(ut.TestCase):
    """Test calc_duration()

    TODO: implement

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def __init__(self, *args, **kwargs):
        super(TestCalcDuration, self).__init__(*args, **kwargs)


class TestCompressShape(ut.TestCase):
    """Test compress_shape()

    TODO: implement

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def __init__(self, *args, **kwargs):
        super(TestCompressShape, self).__init__(*args, **kwargs)


class TestMakeDelay(ut.TestCase):
    """Test make_delay()

    Author: Stefan Kroboth <stefan.kroboth@uniklinik-freiburg.de>
    """
    def __init__(self, *args, **kwargs):
        super(TestMakeDelay, self).__init__(*args, **kwargs)
        self.num = abs(random())

    def test_type(self):
        self.assertIsInstance(make_delay(self.num),  Delay)

    def test_basic(self):
        self.assertEqual(make_delay(self.num).delay,  self.num)

    def test_negative_delay(self):
        self.assertRaises(ValueError, make_delay, -self.num)

    def test_infinite_delay(self):
        self.assertRaises(ValueError, make_delay, np.inf)

    def test_nan_delay(self):
        self.assertRaises(ValueError, make_delay, np.nan)

    def test_zero_delay(self):
        self.assertRaises(ValueError, make_delay, 0)


if __name__ == '__main__':
    ut.main()

# -*- coding: utf-8 -*-
"""

@author: Stefan Kroboth
"""

import unittest as ut
from mr import convert
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


if __name__ == '__main__':
    ut.main()

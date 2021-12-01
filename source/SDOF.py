# -*- coding: utf-8 -*-
"""
Created on Mon Nov 29 14:19:59 2021

@author: stone
"""

import numpy as np

class SDOF():
    def __init__(self, mass, stiffness, damping_ratio):
        self.m = mass
        self.k = stiffness
        self.xi = damping_ratio
        self.set_damping()
        
    def set_damping(self):
        self.omega = np.sqrt(self.k / self.m)
        self.c_crit = self.m * self.omega
        self.c = self.xi * self.c_crit
        
    def solve(self, solver='CDM'):
        pass
# -*- coding: utf-8 -*-
"""
Created on Mon May 11 15:58:36 2020

@author: Benjamin
"""

import numpy as np

def mkGrid_XI_YI(coord_x,coord_y,nbe=500):
    """ grid for interpolation """
    Xm = np.linspace(min(coord_x), max(coord_x), nbe)
    Ym = np.linspace(min(coord_y), max(coord_y), nbe)
    (XI,YI) = np.meshgrid(Xm, Ym)
    
    return XI, YI
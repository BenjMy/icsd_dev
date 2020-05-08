# -*- coding: utf-8 -*-
"""
Created on Fri May  8 10:37:46 2020

@author: Benjamin
"""

import os 

Maindir='E:/Padova/Software/SourceInversion/icsd_dev/testAPI'
os.chdir(Maindir)

path2files= './run_nail/'
from icsd3d_class import iCSD3d_Class as i3d

icsd=i3d(dirName=path2files)   
icsd.showResultsFini(method='Pearson')

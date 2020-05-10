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
icsd.estimateM0(method='F1',show=False)

icsd.inix0='None' # None or cst
icsd.x0_prior=False
icsd.x0_ini_guess=False # initial guess
icsd.wr=1
icsd.plotElecs=False
# icsd3d_Sens.clim=[0,0.1]
# icsd3d_Sens.title=SimName

icsd.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))

icsd.alphaSxy=False
icsd.alphax0=1 # weight of relative smallness
icsd.alphaSx=1 # x-smooth
icsd.alphaSy=0.2 # y-smooth
# icsd3d_Sens.mesh='Figures/mesh3d_rhomap.vtk'      

icsd.icsd_init() 
icsd.run_single()

icsd.showResults()

# %%
# -----------------------------------#
# Example 3d block sensitivity analysis
# -----------------------------------#

maindir='E:/Padova/Simulation/MALM_SensitivityAnalysis/Sens3d/'
os.chdir(maindir)


depth = 10
HAno=2.5
depthAno = -(depth-HAno/2)  ##-depth/1.2 -HAno/2
width = 130
widthAno=10
SoilR=10
AnoR=1

SimName='M' + 'SoilR' + str(SoilR) + 'AnoR' + str(AnoR) + 'Z' + str(depthAno) + 'L' + str(widthAno) +  'h' + str(HAno)

path2files= './' + SimName + '/Data/'

from icsd3d_class import iCSD3d_Class as i3d
import pyvista as pv

icsd3d_SensApriori=i3d(dirName=path2files)   
icsd3d_SensApriori.type='3d'
icsd3d_SensApriori.sim='SIM_sens3d.txt'
icsd3d_SensApriori.obs='OBS_sens3d.txt'
icsd3d_SensApriori.estimateM0(method='Pearson',show=True)


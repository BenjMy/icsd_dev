"""
Inversion of current source density apply to a salt tracer
----------------------------------------------------------
"""
import os
import numpy as np

maindir='E:/Padova/Software/SourceInversion/icsd_dev/'
os.chdir(maindir)


# -----------------------------------#
# Example time-lapse salt tracer
# show the importance of an anisotropic spatial regularisation
# -----------------------------------#
from icsd3d_class import iCSD3d_Class as i3d
import pyvista as pv

tstep=1
path2files="examples/Salt_tracer/t"+str(tstep) +'/'

icsd3d_Salt=i3d(dirName=path2files)   
icsd3d_Salt.type='2d'
icsd3d_Salt.sim="SMALMtimeReg"+str(tstep)+".txt"
icsd3d_Salt.obs="OMALMtimeReg"+str(tstep)+"_synt.txt"
icsd3d_Salt.coord_file="VRTeCoord.txt"
icsd3d_Salt.regMesh='strc'
icsd3d_Salt.x0_prior=False
icsd3d_Salt.x0_ini_guess=False # initial guess
# icsd3d_TL_RWU.wr=0
icsd3d_Salt.plotElecs=False
# icsd3d_Salt.clim=[0,0.1]

icsd3d_Salt.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))

icsd3d_Salt.alphaSxy=False
icsd3d_Salt.alphax0=1 # weight of relative smallness
icsd3d_Salt.alphaSx=1 # x-smooth
icsd3d_Salt.alphaSy=1 # y-smooth
# icsd3d_TL_RWU.mesh='Figures/ERTmodel_reg_'+str(tstep)+'.vtk'      

icsd3d_Salt.icsd_init() 
# icsd3d_Salt.run_single()

icsd3d_Salt.pareto_MinErr=0.1
icsd3d_Salt.pareto_MaxErr=200
icsd3d_Salt.knee=True
icsd3d_Salt.run_pareto()

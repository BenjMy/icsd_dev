"""
Inversion of current source density apply to a salt tracer
----------------------------------------------------------
"""
import os
import numpy as np
import matplotlib.pyplot as plt

# maindir='E:/Padova/Software/SourceInversion/branch_icsd_dev/'
# os.chdir(maindir)


# -----------------------------------#
# Example time-lapse salt tracer
# show the importance of an anisotropic spatial regularisation
# -----------------------------------#
from icsd3d_class import iCSD3d as i3d
import pyvista as pv

tstep=1
path2files="./Salt_tracer/t"+str(tstep) +'/'

icsd_Salt=i3d(dirName=path2files)   
icsd_Salt.type='2d'
# icsd3d_Salt.sim="SMALMtimeReg"+str(tstep)+".txt"
# icsd3d_Salt.obs="OMALMtimeReg"+str(tstep)+"_synt.txt"
icsd_Salt.createSurvey(fname_obs="OMALMtimeReg"+str(tstep)+"_synt.txt",fname_sim="SMALMtimeReg"+str(tstep)+".txt")

icsd_Salt.coord_file="VRTeCoord.txt"
icsd_Salt.regMesh='strc'
icsd_Salt.x0_prior=False
icsd_Salt.x0_ini_guess=False # initial guess
# icsd3d_TL_RWU.wr=0
icsd_Salt.plotElecs=True
# icsd3d_Salt.clim=[0,0.1]

icsd_Salt.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))

icsd_Salt.alphaSxy=False
icsd_Salt.alphax0=1 # weight of relative smallness
icsd_Salt.alphaSx=1 # x-smooth
icsd_Salt.alphaSy=1 # y-smooth
# icsd3d_TL_RWU.mesh='Figures/ERTmodel_reg_'+str(tstep)+'.vtk'      

#icsd3d_Salt.icsd_init() 
icsd_Salt.invert(x0_prior=False)
fig, ax = plt.subplots()
icsd_Salt.showResults(ax=ax)
plt.show()


# icsd3d_Salt.pareto_MinErr=0.1
# icsd3d_Salt.pareto_MaxErr=200
# icsd3d_Salt.knee=True
# icsd3d_Salt.run_pareto()

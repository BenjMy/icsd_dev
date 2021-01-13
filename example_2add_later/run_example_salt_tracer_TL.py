"""
Inversion of current source density apply to a salt tracer
----------------------------------------------------------
"""
import os
import numpy as np

maindir='E:/Padova/Software/SourceInversion/branch_icsd_dev/'
os.chdir(maindir)


# -----------------------------------#
# Example time-lapse salt tracer
# show the importance of an anisotropic spatial regularisation
# -----------------------------------#
from icsd3d_class import iCSD3d as i3d
import pyvista as pv

path2files='example_2add_later/Salt_tracer/all/'

 # tstepmax = 5
# obs_fnames = [for i in range(tstepmax)]
icsd3d_Salt = i3d(dirName=path2files)

obs_fnames = list(['OMALMtimeReg0_synt.txt',
                   'OMALMtimeReg1_synt.txt',
                   'OMALMtimeReg2_synt.txt'])
sim_fnames = list(['SMALMtimeReg0.txt',
                   'SMALMtimeReg1.txt',
                   'SMALMtimeReg2.txt'])

surveys = icsd3d_Salt.createTimeLapseSurvey(obs_fnames,sim_fnames)

surveys[0].obs
surveys[1].obs

# m0 = icsd3d_Salt.estimateM0(method_m0='F1',show=True)
icsd3d_Salt.TL = True
icsd3d_Salt.invert(x0_prior=False,wr=1e3)

# icsd3d_Salt.invert(x0_prior=False,wr=1e3,TL=True)



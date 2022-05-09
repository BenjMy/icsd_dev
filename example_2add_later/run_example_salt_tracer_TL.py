"""
Inversion of current source density apply to a salt tracer
----------------------------------------------------------
"""
import os
import numpy as np

# maindir='E:/Padova/Software/SourceInversion/branch_icsd_dev/'
# os.chdir(maindir)


# -----------------------------------#
# Example time-lapse salt tracer
# show the importance of an anisotropic spatial regularisation
# -----------------------------------#
from icsd3d_class import iCSD3d as i3d
import pyvista as pv

path2files='./Salt_tracer/all/'

 # tstepmax = 5
# obs_fnames = [for i in range(tstepmax)]
icsdTL_Salt = i3d(dirName=path2files)
icsdTL_Salt.logTrans = False

obs_fnames = list(['OMALMtimeReg0_synt.txt',
                   'OMALMtimeReg1_synt.txt',
                   'OMALMtimeReg2_synt.txt'])
sim_fnames = list(['SMALMtimeReg0.txt',
                   'SMALMtimeReg1.txt',
                   'SMALMtimeReg2.txt'])

surveys = icsdTL_Salt.createTimeLapseSurvey(obs_fnames,sim_fnames)


# m0 = icsd3d_Salt.estimateM0(method_m0='F1',show=True)
icsdTL_Salt.TL = True
icsdTL_Salt.invert(wr=1)

fig, axs = plt.subplots(3, 1, sharex='all', sharey='all',figsize=(20,5))
for i, j in  enumerate(range(3)):
    icsdTL_Salt.showResults(ax=axs[i],index=j)
# plt.savefig(icsdPath+'icsd_Vs.png', dpi=400)




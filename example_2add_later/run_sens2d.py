import os
import numpy as np

maindir='E:/Padova/Simulation/MALM_SensitivityAnalysis/Sens2d/'
os.chdir(maindir)

# ------------------------------------------------#
# Example 2d simple synthetic sensitivity analysis
# ------------------------------------------------#
from icsd3d_class import iCSD3d_Class as i3d
from plotters import mpl_plot

SoilR=1
AnoR=1
AnoRstr = str(AnoR).replace(".", "")
CorERTBack=1
ertSchemeFull=0
SimName='M' + 'SoilR' + str(SoilR) + 'AnoR' + str(AnoRstr) +  'Cert' + str(CorERTBack==1)  + 'ERTscheme' + str(ertSchemeFull==1) + '/'


path2files= './' + SimName
icsd3d_Sens=i3d(dirName=path2files)   
icsd3d_Sens.type='2d'
icsd3d_Sens.sim='MSoilR'+ str(SoilR) + 'AnoR' + AnoRstr + 'Cert' + str(CorERTBack==1) + 'ERTscheme' + str(ertSchemeFull==1) + ".txt"
icsd3d_Sens.obs='OMSoilR'+ str(SoilR) +'AnoR'+AnoRstr + 'Cert' + str(CorERTBack==1) + 'ERTscheme' + str(ertSchemeFull==1) + "_synt.txt"

# show observation data
# mpl_plot.showObs2d(path2files,filename=icsd3d_Sens.obs)

# # Unconstrainsted current source densities inversion
# icsd3d_Sens.invert(regMesh='unstrc',wr=10,x0_prior=True)
 
# # Estimate initial model
# icsd3d_Sens.estimateM0(methodM0='F1', show=True)

# # Constrainsted current source densities inversion
# icsd3d_Sens.invert(wr=0.01,x0_prior=True,alphax0=1)

# # Change smoothing to anisotropic
# icsd3d_Sens.invert(alphaSxy=True,alphaSx=10, alphaSy=0.2,
#                     wr=1)

# Pareto curve analysis
icsd3d_Sens.invert(regMesh='unstrc',x0_prior=True, 
                   pareto=False)
# icsd3d_Sens.ModelResolution(jacMi=305)

# icsd3d_Sens.showResults()
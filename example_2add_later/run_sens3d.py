import os
import numpy as np

maindir='E:/Padova/Simulation/MALM_SensitivityAnalysis/Sens3d/'
os.chdir(maindir)

ZZ = -18.75
Rsoil = 1
SimName = 'MSoilR' + str(Rsoil) + 'AnoR1Z'+ str(ZZ) +'L10h2.5'

path2files= './' + SimName + '/Data/'

# -----------------------------------#
# Example block sensitivity analysis
# -----------------------------------#
from icsd3d_class import iCSD3d_Class as i3d

# import mpl_plot
import pyvista as pv

icsd3d_SensApriori=i3d(dirName=path2files)   
icsd3d_SensApriori.type='3d'
icsd3d_SensApriori.sim='SIM_sens3d.txt'
icsd3d_SensApriori.obs='OBS_sens3d.txt'

icsd3d_SensApriori.estimateM0(methodM0='F1',show=True)
icsd3d_SensApriori.estimateM0(methodM0='Pearson',show=True)

icsd3d_Sens=i3d(dirName=path2files)   
icsd3d_Sens.type='3d'
icsd3d_Sens.sim='SIM_sens3d.txt'
icsd3d_Sens.obs='OBS_sens3d.txt'
icsd3d_Sens.coord_file='VRTeCoord.txt'
icsd3d_Sens.regMesh='unstrc'
icsd3d_Sens.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))

# Unconstrainsted current source densities inversion
icsd3d_Sens.invert(pareto=False,regMesh='strc',wr=1,x0_prior=False)

# Unconstrainsted current source densities inversion
icsd3d_Sens.invert(pareto=False, regMesh='unstrc',wr=1,x0_prior=True)

# constrainsted current source densities inversion
icsd3d_Sens.invert(pareto=True, regMesh='strc',x0_prior=False, pareto_MinErr=1, pareto_MaxErr=100)

# sol = icsd3d_Sens.invert(x0_prior=True, pareto=False)
# # J = sol.jac[64 + 1:,100]
# # plt.plot(J)

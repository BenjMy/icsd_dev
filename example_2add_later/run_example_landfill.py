"""
Inversion of current source density apply to a landfill leakage
---------------------------------------------------------------
"""

import os
import numpy as np

maindir='E:/Padova/Software/SourceInversion/icsd_dev/'
os.chdir(maindir)


# -----------------------------------#
# Exemple with a 3d landfill geometry --> ARTICLE comp. Geosciences ?
# -----------------------------------#
# a delineaer 
from icsd3d_class import iCSD3d_Class as i3d
from plotters import mpl_plot

# No_hole Ano1_bh_lambda_1
path2files="example_2add_later/Landfill_3d/Ano_0_E13/" # Test elec outside landfill
# path2files="example_2add_later/Landfill_3d/Ano_0_E89/" # Test 2 elec outside landfill
# path2files="example_2add_later/Landfill_3d/Ano_0_EA/"  # Real A position without anomaly
# path2files="example_2add_later/Landfill_3d/Ano_1_BH_EA/"  # Real A position with big hole anomaly
# path2files="example_2add_later/Landfill_3d/Ano_1_SH_EA/"  # Real A position  with small hole anomaly
# path2files="examples/Landfill_3d/RealData_Ano_0/"  # Real data / without anomaly
# path2files="examples/Landfill_3d/RealData_Ano_1_SH/"  # Real data / with (small) anomaly

# icsd3d_landfill.pareto

icsd3d_landfill=i3d(dirName=path2files)   
icsd3d_landfill.type='3d'
icsd3d_landfill.sim="SNoAno.txt" # or SNoAno / SAno
icsd3d_landfill.obs="ONoAno_synt.txt" # or ONoAno_synt / OAno_synt / Obs_real
icsd3d_landfill.coord_file="VRTeCoord.txt"

# icsd3d_landfill.estimateM0(methodM0='Pearson',show=True)
# icsd3d_landfill.invert(plotElecs=True,mesh='mesh3d_rhomap.vtk')


# Unconstrainsted current source densities inversion
icsd3d_landfill.invert(wr=1,x0_prior=False)
 
# Estimate initial model
icsd3d_landfill.estimateM0(methodM0='F1', show=True)

# Constrainsted current source densities inversion
icsd3d_landfill.invert(regMesh='unstrc',wr=100,x0_prior=True)

icsd3d_landfill.invert(x0_prior=False,)

# Pareto curve analysis
icsd3d_landfill.invert(x0_prior=False, 
                   pareto=True, pareto_MinErr=1, pareto_MaxErr=200, pareto_nSteps=10)

icsd3d_landfill.showResults()

# xx1,yy1= [284146,	5036300]
# xx2,yy2= [284281,	5036210]

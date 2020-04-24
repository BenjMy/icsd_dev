import os
import numpy as np

maindir='E:/Padova/Software/SourceInversion/icsd_dev/'
os.chdir(maindir)


# -----------------------------------#
# Exemple with a 3d landfill geometry --> ARTICLE comp. Geosciences ?
# -----------------------------------#
# a delineaer 
from icsd3d_class import iCSD3d_Class as i3d
# No_hole Ano1_bh_lambda_1
# path2files="examples/Landfill_3d/Ano_0_E13/" # Test elec outside landfill
# path2files="examples/Landfill_3d/Ano_0_E89/" # Test 2 elec outside landfill
# path2files="examples/Landfill_3d/Ano_0_EA/"  # Real A position without anomaly
# path2files="examples/Landfill_3d/Ano_1_BH_EA/"  # Real A position with big hole anomaly
path2files="examples/Landfill_3d/Ano_1_SH_EA/"  # Real A position  with small hole anomaly
# path2files="examples/Landfill_3d/RealData_Ano_0/"  # Real data / without anomaly
# path2files="examples/Landfill_3d/RealData_Ano_1_SH/"  # Real data / with (small) anomaly


icsd3d_landfill=i3d(dirName=path2files)   
icsd3d_landfill.type='3d'
icsd3d_landfill.sim="SAno.txt" # or SNoAno / SAno
icsd3d_landfill.obs="OAno_synt.txt" # or ONoAno_synt / OAno_synt / Obs_real
icsd3d_landfill.coord_file="VRTeCoord.txt"
icsd3d_landfill.regMesh='unstrc'
icsd3d_landfill.sc=1

# icsd3d_landfill.alphaSxy=True
# icsd3d_landfill.alphax0=1 # weight of relative smallness
# icsd3d_landfill.alphaSx=1 # x-smooth
# icsd3d_landfill.alphaSy=1 # y-smooth

icsd3d_landfill.x0_prior=True
icsd3d_landfill.x0_ini_guess=True # initial guess
# icsd3d_landfill.wr=1
icsd3d_landfill.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
icsd3d_landfill.mesh='mesh3d_rhomap.vtk'
icsd3d_landfill.plotElecs=True  
icsd3d_landfill.gif3d=False     
icsd3d_landfill.icsd_init()  
icsd3d_landfill.run_single()


icsd3d_landfill.invert(pareto=True,showfig=False)
icsd3d_landfill.showResults(clim=[0,0.0001],gif3d=True)
icsd3d_landfill.DataImport(clim=[0,0.0001],gif3d=True)


# icsd3d_landfill.pareto_MinErr=0.01
# icsd3d_landfill.pareto_MaxErr=100
# icsd3d_landfill.knee=True
# icsd3d_landfill.run_pareto()


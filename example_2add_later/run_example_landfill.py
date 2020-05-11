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
# No_hole Ano1_bh_lambda_1
# path2files="example_2add_later/Landfill_3d/Ano_0_E13/" # Test elec outside landfill
# path2files="example_2add_later/Landfill_3d/Ano_0_E89/" # Test 2 elec outside landfill
path2files="example_2add_later/Landfill_3d/Ano_0_EA/"  # Real A position without anomaly
# path2files="example_2add_later/Landfill_3d/Ano_1_BH_EA/"  # Real A position with big hole anomaly
# path2files="example_2add_later/Landfill_3d/Ano_1_SH_EA/"  # Real A position  with small hole anomaly
# path2files="examples/Landfill_3d/RealData_Ano_0/"  # Real data / without anomaly
# path2files="examples/Landfill_3d/RealData_Ano_1_SH/"  # Real data / with (small) anomaly

icsd3d_landfill=i3d(dirName=path2files)   
icsd3d_landfill.type='3d'
icsd3d_landfill.sim="SNoAno.txt" # or SNoAno / SAno
icsd3d_landfill.obs="ONoAno_synt.txt" # or ONoAno_synt / OAno_synt / Obs_real
icsd3d_landfill.coord_file="VRTeCoord.txt"

icsd3d_landfill.estimateM0(method='Pearson',show=True)


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

icsd3d_landfill.showResults()

xx1,yy1= [284146,	5036300]
xx2,yy2= [284281,	5036210]

def path(y):
    """Equation: x = a(y-h)^2 + k"""
    a = 1
    x = a * y + 0.0
    return x, y

cmap = plt.cm.get_cmap("viridis", 4)

xx1,yy1= [284146,	5036300]
xx2,yy2= [284281,	5036210]

def path(y):
    """Equation: x = a(y-h)^2 + k"""
    a = 1
    x = a * y + yy1
    return x, y

x, y = path(np.arange(xx2, xx1, 15.0))
zo = np.linspace(9, 11, num=len(y))

x, y = path(np.arange(model.bounds[2], model.bounds[3], 15.0))
zo = np.linspace(9.0, 11.0, num=len(y))
points = np.c_[x, y, zo]
spline = pv.Spline(points, 15)
spline
        
slc = model.slice_along_line(spline)
slc

p = pv.Plotter()
p.add_mesh(slc, cmap=cmap)
p.add_mesh(model.outline())
p.show_grid()
p.show(cpos=[1, -1, 1])


# import pyvista as pv
# from pyvista import examples
# import matplotlib.pyplot as plt
# import numpy as np
# model = examples.load_channels()


# def path(y):
#     """Equation: x = a(y-h)^2 + k"""
#     a = 1
#     x = a * y + 0.0
#     return x, y

# x1,y1= [284146,	5036300]
# x2,y2= [284281,	5036210]

# line=pv.Line(pointa=[x1, y1, 0.0], pointb=[x2, y2, 0.0], resolution=1)

# zo = np.linspace(9.0, 11.0, num=len(y))
# points = np.c_[x, y, zo]
# spline = pv.Spline(points, 15)


# slc = model.slice_along_line(line)


x, y = path(np.arange(y1, y2, 15.0))
zo = np.linspace(9.0, 11.0, num=len(y))
points = np.c_[x, y, zo]
spline = pv.Spline(points, 15)
slc = model.slice_along_line(spline)


# icsd3d_landfill.showResults(cut=True)

# icsd3d_landfill.invert(pareto=True,showfig=False)
# icsd3d_landfill.showResults(clim=[0,0.0001],gif3d=True)
# icsd3d_landfill.DataImport(clim=[0,0.0001],gif3d=True)


# icsd3d_landfill.pareto_MinErr=0.01
# icsd3d_landfill.pareto_MaxErr=100
# icsd3d_landfill.knee=True
# icsd3d_landfill.run_pareto()


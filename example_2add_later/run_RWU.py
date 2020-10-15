# # -----------------------------------#
# # Example time-lapse data SOIL article
# #-> Importance of the smallness minimisation
# # -----------------------------------#
# from icsd3d_class import iCSD3d_Class as i3d
# import pyvista as pv
# import matplotlib.pyplot as plt

# path2files="examples/RWU_TimeLapse/ADAM/"

# icsd3d_TL_RWU=i3d(dirName=path2files)   
# icsd3d_TL_RWU.type='3d'
# icsd3d_TL_RWU.sim="VRTeSimADAMMS0.txt"
# icsd3d_TL_RWU.obs="ObsData_ADAMMS0.txt"
# icsd3d_TL_RWU.coord_file="VRTeCoord.txt"
# icsd3d_TL_RWU.regMesh='unstrc'
# icsd3d_TL_RWU.x0_prior=True  
# icsd3d_TL_RWU.alphax0=1 # weight of relative smallness
# icsd3d_TL_RWU.x0_ini_guess=True # initial guess
# # icsd3d_TL_RWU.wr=1
# icsd3d_landfill.obs_err='sqrt' # choose between constant weight and w = 1/sqrt(abs(obs))
# icsd3d_TL_RWU.plotElecs=True
# icsd3d_TL_RWU.icsd_init() 
# # icsd3d_TL_RWU.mesh='mesh3d_rhomap.vtk'      
# # icsd3d_TL_RWU.run_single()


# icsd3d_TL_RWU.pareto_MinErr=0.1
# icsd3d_TL_RWU.pareto_MaxErr=200
# icsd3d_TL_RWU.knee=True
# icsd3d_TL_RWU.run_pareto()
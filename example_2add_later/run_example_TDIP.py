"""
Inversion of current source density: application on TDIP data
------------------------------------------------------------
"""
import os
import matplotlib.pyplot as plt
# from TDIP.TDIP_MALM import  TDIP_Green_simulation

# -----------------------------------#
# Exemple TDIP
# -----------------------------------#

path2files="./TDIP/"

# -----------------------------------#
# apply here a time lapse contrainst to the regularisation
from icsd3d_class import iCSD3d as i3d
from plotters import mpl_plot

# sim = TDIP_Green_simulation()
fname_obs = 'MALMIP1217.bin'
fname_sim = 'VRTeSim.txt'

from importers.read import *
# tdip_obs = loadTDIPSurvey(path2files + fname_obs) # *.data (pygimli format)
# tdip_obs.MA
# tdip_sim = load_sim(self.dirName + fname_sim) # *.data (pygimli format)
        
icsd_TDIP = i3d(dirName=path2files)   
testTDIP = icsd_TDIP.createTDIPSurvey(fname_obs, # Observation TDIP file
                           fname_sim) #Simulation primary voltages only
m0 = icsd_TDIP.estimateM0(method_m0='F1',show=True)


testTDIP[0].TDIP_flag
testTDIP[0].dirName
testTDIP[0].regMesh
testTDIP[0].A

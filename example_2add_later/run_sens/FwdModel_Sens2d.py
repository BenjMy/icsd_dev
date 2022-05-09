#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
    Extended minimal pygimli example to simulate Darcy velocity,
    mass transport and time-lapse ERT measurements
"""
import os 

MainPath='E:/Padova/Experiments/8_2019_Rhizotron_DAFNAE_Lancaster_ERT_MALM/'
os.chdir(MainPath)

from Z_Repo_Codes_Rhizotron.MALM4RootsRhizo_class import MALM4RootsRhizo_class as MR
from Z_Repo_Codes_Rhizotron.MALM4RootsRhizo_class import Read_class as RC


os.chdir(r'E:\Padova\Simulation\icsd_SalineTracer')

import numpy as np
# import pybert as pb

from pygimli.mplviewer import drawStreams
from pygimli.mplviewer import drawSensors
import pygimli as pg
import pygimli.meshtools as mt
from pygimli.physics.ert import ERTManager
import pygimli.physics.ert as ert
import pygimli.physics.petro as petro

import matplotlib.pyplot as plt

import itertools


PathSimu_salt=r'E:\Padova\Simulation\MALM_SensitivityAnalysis\Sens2d'
os.chdir(PathSimu_salt)

SoilR=1000
AnoR=1
AnoRstr = str(AnoR).replace(".", "")
CorERTBack=1
ertSchemeFull=0
# SimName='M' + 'SoilR' + str(SoilR) + 'AnoR' + str(AnoRstr) +  'Cert' + str(CorERTBack==1)  + 'ERTscheme' + str(ertSchemeFull==1) 
SimName='M' + 'SoilR' + str(SoilR) + 'AnoR' + str(AnoRstr) +  'Cert' + str(CorERTBack==1)  + 'ERTscheme' + str(ertSchemeFull==1)

path_Figures = PathSimu_salt  + '/Sens/' + SimName + '/Figures/'
if not os.path.exists(path_Figures):
        os.makedirs(path_Figures)

path_Data = PathSimu_salt + '/Sens/'  + SimName + '/Data/'
if not os.path.exists(path_Data):
        os.makedirs(path_Data)
    
    
# Create geometry definition for the modelling domain
world = mt.createWorld(start=[-20, 0], end=[20, -16],worldMarker=False)
# Create a heterogeneous block
block = mt.createRectangle(start=[-5, -2.5], end=[5, -5.0],
                            marker=4,  boundaryMarker=10, area=0.1)
# block = mt.createRectangle(start=[-5, -5], end=[5, -7.5],
#                             marker=4,  boundaryMarker=10, area=0.1)
# Merge geometrical entities
# geom = mt.mergePLC([world, block])
geom = mt.mergePLC([world, block])

pg.show(geom, boundaryMarker=True, savefig='geometry.pdf')

# Create a mesh from the geometry definition
mesh = mt.createMesh(geom, quality=32, area=0.2, smooth=[1, 10])
# pg.show(mesh, savefig='mesh.pdf')

rMap = [[1, SoilR], [4, AnoR]]

# Map conductivity value per region to each cell in the given mesh
R = pg.solver.parseMapToCellArray(rMap, mesh)
# Dirichlet conditions for hydraulic potential


# Create survey measurement scheme
# ertScheme = ert.createERTData(pg.utils.grange(-20, 20, dx=1.0),
#                               schemeName='dd')


if ertSchemeFull==1:
    a=np.arange(-9.3,10,1)
    b=np.arange(-6.3,0,1)
    c = list(itertools.product(a, b))
    elecs= np.c_[np.array(c),np.zeros(len(c))]
    # orig_scheme = pg.DataContainerERT()
    # orig_scheme.setSensorPositions(elecs)
    len(elecs)      
    ertScheme = ert.createERTData(elecs=elecs,
                              schemeName='dd')
else:
    ertScheme = ert.createERTData(pg.utils.grange(-20, 20, dx=1.0),
                              schemeName='dd')
ertScheme.save(path_Data + 'SchemeERT_dd_2d_SENS.data')

# Create suitable mesh for ert forward calculation
# meshERT = mt.createParaMesh(ertScheme, quality=33, paraMaxCellSize=0.2,
#                             boundaryMaxCellSize=50, smooth=[1, 2])

## --------- make sequence MALM --------- ##

sensors = ertScheme.sensorPositions()
EA = len(sensors) -1
EB = len(sensors) -2
EN = len(sensors) -3

sensors_n= np.array(sensors)
# Combine two injection - 2 current sources
# sensors_n[EA,0:2]=[-3.5, -1]
# sensors_n[EA,0:2]=[3, -6.25]

# injection in the conductive boday
sensors_n[EA,0:2]=[0, -3.25]

sensors_n[EB,0:2]=[-150, 0]
sensors_n[EN,0:2]=[150, 0]

fig, ax = plt.subplots()
pg.show(geom,ax=ax)
drawSensors(ax, sensors_n, diam=0.6)
# ax.figure.savefig(directoryFig +'/'
#                   +'MeshERTFwd.png',dpi=350)
# mesh_root.save(path_Scheme +'/'+ SimuName + '/'  '/MeshERTFwd.bms')
# drawSensors(ax, sensors=sensors_n, diam=0.01)

# Create suitable mesh for ert forward calculation

World_VRTs = mt.createWorld(start=[-160, -160], end=[160, 0]) #Extend world if pole pole acquisition in the field ; Take care of boundary conditions !
geom_VRTs = mt.mergePLC([World_VRTs,geom])
# geom = mt.mergePLC([world, block])
pg.show(geom_VRTs)

## --------- create vRTE positions in regular grid 2d --------- ##
dx = (block.xmax()- block.xmin())/15
xreg = np.arange(min(sensors_n[:-3,0]), max(sensors_n[:-3,0])+dx, dx, 'float')
yreg = np.arange(-10, 0, dx, 'float')
print('Create quadrilateral mesh suitable for iscd -- Virtual sources must no coincide with electrodes !')
#     Create quadrilateral mesh suitable for iscd
grid_vrte = pg.Mesh(2)
grid_vrte.create2DGrid(xreg, yreg, 0)
        
# add virtual sources to the geom        
for pts in sensors_n[:-3]:
    geom_VRTs.createNode(pts,marker=-99)
    geom_VRTs.createNode(pts + pg.RVector3(0, -0.1))

# add remotes electrode to the geom        
for pts in sensors_n[-3:-1,:]:
    geom_VRTs.createNode(pts,marker=-100)
    geom_VRTs.createNode(pts + pg.RVector3(0, -0.1))
    
# add malm A electrode to the geom        
geom_VRTs.createNode(sensors_n[-1],marker=-101)
geom_VRTs.createNode(sensors_n[-1] + pg.RVector3(0, -0.1))
    
# add virtual sources to the geom        
Pos_vrtS=np.array([])
for node in grid_vrte.nodes():
    print(node.pos())
    Pos_vrtS = np.append(Pos_vrtS,np.array(node.pos()))
    geom_VRTs.createNode(node.pos(),marker=-991)
    geom_VRTs.createNode(node.pos() + pg.RVector3(0, -0.1))
Pos_vrtS = np.reshape(Pos_vrtS, [len(grid_vrte.nodes()),3])

# mesh = mt.createMesh(geom, quality=32, area=0.2, smooth=[1, 10])

meshMALM_ERT = mt.createMesh(geom_VRTs, quality=33, smooth=[1, 10]) #refine close to the electrodes
# mesh_VRTs.node(RemoteNPos).setMarker(-99)
# pg.show(meshMALM_ERT)

malmScheme = pg.DataContainerERT()
meas = []
dataABMN= []
for j in range(len(sensors)-3):
    M=j
    meas.append([EA, EB,M, EN])

dataABMN = np.array(meas)
malmScheme.resize(len(dataABMN))        

for i, j in enumerate("abmn"):
    malmScheme.set(j, dataABMN[:, i])
    
malmScheme.setSensorPositions(sensors_n)
malmScheme.set("valid", np.ones(len(meas)))
malmScheme.save(path_Data+'SchemeMALM_2d_SENS.data')

#See MALM4ROOTSmain file to insert VRTE into the geom and not into the mesh

        
## ---------  return VRTE file format icsd MarkerVRTE=991 --------- ##
VRTEpos= MR.VRTEpos(mesh=meshMALM_ERT,dim=2)
print('Nb of virtual sources=' + str(len(VRTEpos)))
print('Total nb of sensors=' + str(len(sensors)-3 + len(VRTEpos)))


## ---------  collect and plot all electrode config --------- ##
Elecs, Elec99, Elec991 = MR.SensorsPos(meshMALM_ERT,pltRem=False,pltVRTE=True, dim2=True,savefig=True)
# Elecs = MR.SensorsPos(mesh3d,pltRem=True,pltVRTE=True)

# ax, _ = pg.show(mesh, data=K, label='Hydraulic conductivity $K$ in m$/$s',
#                 cMin=1e-5, cMax=1e-2, nLevs=4, cMap='viridis')
ax, _ = pg.show(meshMALM_ERT)
drawSensors(ax, sensors_n,diam=1)
ax.scatter(Pos_vrtS[:,0],Pos_vrtS[:,1], s=1, c='r', marker='v')
ax.figure.savefig(path_Figures + 'MeshVRTE'+'.png',dpi=350)


## --------- prepare sequence for green functions --------- ##
SeqFullMALM= MR.mk_full_malm(dataABMN, VRTe = range(len(Elec99)+3, len(VRTEpos)+len(Elec99)+3),
                             mesh=meshMALM_ERT, R3=False) # output .shm with sensors
                  

## ---------  Calculate bulk resistivity based on Archie's Law --------- ##

ax, _ = pg.show(meshMALM_ERT, data=rMap, label='ER $\Omega$.m',
                nLevs=2, cMap='viridis')
offset =  1
ax.set_xlim([min(sensors_n[:-3,0]), max(sensors_n[:-3,0])+offset])
ax.set_ylim(-10, 0)
ax.figure.savefig(path_Figures+ 'ERTmodel' + SimName,dpi=350)

# meshMALM_ERT.save('meshMALM.vtk')
meshMALM_ERT.addData("rhomap", R)
meshMALM_ERT.exportVTK(path_Figures+ 'ERTmodel' + SimName+ '.vtk')
# ax, _ = pg.show(meshMALM_ERT)
# offset =  1
# ax.set_xlim([min(sensors_n[:-3,0]), max(sensors_n[:-3,0])+offset])
# ax.set_ylim(-10, 0)

## ---------  Forward simulate ERT and MALM --------- ##
## --------- simulate solution (synthetic case) --------- ##

MALM_SynthSol = MR.SimulateSynthMALMSol(meshMALM_ERT,dataABMN,rMap,Name=SimName)

# resBulk_new?

# plot potential saw by surface electrodes

plt.figure()
plt.plot(sensors_n[:-3,0],np.array(MALM_SynthSol('r')), label=SimName)
plt.xlabel('#Elec')
plt.ylabel('Voltage (V)')
plt.legend(loc='upper right')
plt.tight_layout(.5)
plt.savefig(path_Figures+ SimName + "PotGraph_reg.png",dpi=350)

len(np.array(MALM_SynthSol('r')))
len(sensors_n[:-3,:])



from scipy.interpolate import interp1d
xnew = np.arange(min(sensors_n[:-3,0]), max(sensors_n[:-3,0]),step=0.1)
f = interp1d(sensors_n[:-3,0], np.array(MALM_SynthSol('r')), kind='cubic')
# f2 = interp1d(x, y, kind='cubic')

plt.plot(sensors_n[:-3,0], np.array(MALM_SynthSol('r')), 'o', 
         xnew, f(xnew), '-')

data2write=np.zeros([len(xnew),4])
data2write[:,0]=xnew
data2write[:,3]=f(xnew)

f = open(path_Data + SimName +'_wlt.dat','w')
np.savetxt(f, data2write, fmt='%1.2f\t %d\t %d\t %1.6f\t', delimiter='\t')   # X is an array
f.close()

f = open(path_Data + SimName +'_wlt.txt','w')
np.savetxt(f, data2write[:,[0,1,3]], fmt='%1.3f\t %1.3f\t %1.6f\t', delimiter='\t')   # X is an array
f.close()

            
fig, ax = plt.subplots()

u = MR.SimulateSynthMALMSol(meshMALM_ERT,dataABMN-1,rMap,Name=SimName,returnFields_opt=True)

# len(u)
# u[0]
uT=u[EA-1]-u[EB-1]
# np.array(u[EA-1])
# np.array(u[EB-1])

print('Plot equipotential distribution and quivers')

# min(uT)
cmin=np.percentile(uT, 10)
cmax=np.percentile(uT, 90)
#drawquiver 
gridCoarse = pg.createGrid(x=np.linspace(meshMALM_ERT.xmin(), meshMALM_ERT.xmax(), 20),
                           y=np.linspace(meshMALM_ERT.ymin(), meshMALM_ERT.ymax(), 20))


fig, ax = plt.subplots()
fig.subplots_adjust(top = 0.8)
# drawStreams(ax, mesh, uT, color='red')
pg.show(meshMALM_ERT, data=uT, fillContour=True, colorBar=True, 
          cMin=cmin,cMax=cmax,
          orientation='horizontal', label='Voltage u (V)', nLevs=10,
          logScale=False, hold=True, showMesh=True, ax=ax,linewidths=1)[0]
# pg.show(meshMALM_ERT, data=uT, fillContour=True, colorBar=True, 
#           cMin=cmin,cMax=cmax, orientation='horizontal', label='Voltage u (V)', nLevs=10,
#           logScale=True, hold=True, showMesh=True, xlabel='x(m)', ylabel='y(m)',ax=ax)
# drawStreams(ax, meshMALM_ERT, uT, coarseMesh=gridCoarse, color='red')
# pg.show(geom_VRTs)
ax.set_xlim([min(sensors_n[:-3,0])-offset, max(sensors_n[:-3,0])+offset])
ax.set_ylim(-10, 0)
ax.scatter(sensors_n[EA-1][0],sensors_n[EA-1][1], c='red')
ax.scatter(sensors_n[EB-1][0],sensors_n[EB-1][1])
ax.figure.savefig(path_Figures+ SimName+'_CurrentmapFwdMALM_reg.png',dpi=350)
        


# ## --------- simulate green functions --------- ##
if CorERTBack==1:
    MR.SimulateGreenFcts(mesh_VRTs=meshMALM_ERT,rhomapVRTs=rMap,schemeVrts=SeqFullMALM, 
                          Name=SimName)
else:
    MR.SimulateGreenFcts(mesh_VRTs=meshMALM_ERT,rhomapVRTs=1,schemeVrts=SeqFullMALM, 
                          Name=SimName)







# # apply background resistivity model
# rho0 = np.zeros(meshMALM_ERT.cellCount()) + 1000.
# for c in meshMALM_ERT.cells():
#     if c.center()[1] < -8:
#         rho0[c.id()] = 150.
#     elif c.center()[1] < -2:
#         rho0[c.id()] = 500.
# resis = pg.RMatrix(resBulk)
# # for i, rbI in enumerate(resBulk):
# #     resis[i] = 1. / ((1./rbI) + 1./rho0)
# # Initialize ert method manager
# # ERT = ERTManager(verbose=False)
# # Run  simulation for  the apparent resistivities
# # rhoa = ERT.simulate(meshERT, rho0, ertScheme)
# # rhoa = ERT.simulate(meshERT, resis, ertScheme, verbose=0, returnArray=True)

# resBulk_new.save('test_resis.dat')
# ert = pb.ERTManager()
# u= ert.simulate(mesh=meshMALM_ERT, res=resBulk_new, scheme=ertScheme)
# u.save('ertsimu.data')
# pb.showData(ertScheme, vals=u)
# plt.show()


# malm = pb.ERTManager()



# pg.show()
# #         Scheme_Sol.setSensorPositions(sensors)
# #         Scheme_Sol.set("valid", np.ones(len(SeqSol)))
# #         Scheme_Sol.save('seq_sol_tmp.data')

# #         # Simulate
# # #        schemeMALMSol = pb.DataContainerERT('Scheme_Sol.shm')
# #         MALM_sol = pb.ERTManager()
# #         Sol = MALM_sol.simulate(mesh=mesh, res=rhomap, scheme=Scheme_Sol)
        
# # pg.show(meshERT)
# ertScheme.save('schemeERT.data')
# meshERT.addExportData("rhomap", resis)
# meshERT.exportVTK('mesh2d_rhomap.vtk')

# MR.PyMesh3d(path_Figures_ano + 'mesh3d_rhomap.vtk',scalars='rhomap')


# # Solve the electrical forward problem using the ERT method manager
# axs = pg.plt.subplots(3, 2, sharex=True, sharey=True)[1].flatten()
# for i in range(6):
#     # ERT.showData(ertScheme, vals=rhoa[i*10]/rhoa[0], ax=axs[i])
#     pg.show(ertScheme, vals=rhoa[i]/rhoa[0], ax=axs[i])

# # just hold figure windows open if run outside from spyder, ipython or similar
# pg.wait()

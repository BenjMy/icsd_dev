# -*- coding: utf-8 -*-
"""
Created on Sun Jun 21 13:31:39 2020

@author: Benjamin
"""


import numpy as np
import pygimli as pg
import pygimli.meshtools as mt
# import pybert as pb
from pygimli.physics.ert import simulate as simulateERT
import matplotlib.pyplot as plt
from Z_Repo_Codes_Rhizotron.MALM4RootsRhizo_class import Mesher as MyMS
import Z_Repo_Codes_Rhizotron.MALM4RootsRhizo_class as MR
from fatiando import gridder
import plot_dEXP as pEXP
from pygimli.viewer.pv import drawSensors, drawSlice, pgMesh2pvMesh

import pickle

import os
MainPath = 'E:\Padova\Simulation\MALM_SensitivityAnalysis'
os.chdir(MainPath)

pyvista = pg.optImport("pyvista")

depth = 60  # mesh depth
width = 350 # mesh width
# -------------------
HAno = 2.5  # 2.5, 12.5, 22.5 anomaly thickness
depthAno = -13.75 # 3.75, 13.75 or 23.75
widthAno = 15 # 5, 15, 25
shift = 0
lengthAno = 5
# -------------------
noiseLevel = 0 #0.01 #1 
noiseAbs= 0 #1e-6 # 1e-6
# seed=1337
# -------------------

ERmodel=True
SoilR=1
AnoR=1
rMap = [[1, SoilR], [4, AnoR]]
add_bor=0

# comments bugs
# 3d world area parameters not working
# quality 33 not good compare to 1.3
# noise not working with returnfield option
                    
# %% create path figure and data

SimName=('M' + 'SoilR' + str(SoilR) + 'AnoR' + str(AnoR) + 'Z' + str(depthAno) 
        + 'W' + str(widthAno) +  'H' + str(HAno) + 'L' + str(lengthAno) + 'S' + str(shift)
        + 'Noise' + str(noiseLevel))
path_Figures = MainPath  + '/Sens3d/' + SimName + '/Figures/'
if not os.path.exists(path_Figures):
    os.makedirs(path_Figures)
path_Data = MainPath + '/Sens3d/'  + SimName + '/Data/'
if not os.path.exists(path_Data):
    os.makedirs(path_Data)
os.chdir(MainPath + '/Sens3d/'  + SimName)

# %% create model geometry

# world = mt.createWorld(start=[-width/2, -width/2, 0], end=[width/2, width/2, -depth],
#                         worldMarker=True) # ,area=10 markerPosition=[width,width]
world = mt.createCube(size=[width, width, depth], pos=[width/2, width/2, -depth/2],marker=1,area=75) #, boundaryMarker=0)

for i, b in enumerate(world.boundaries()):
    # if worldMarker is True:
    if b.norm()[2] == 1.0:
        b.setMarker(pg.core.MARKER_BOUND_HOMOGEN_NEUMANN)
    else:
        b.setMarker(pg.core.MARKER_BOUND_MIXED)
                    
ano = mt.createCube(size=[widthAno, lengthAno, HAno], pos=[width/2+shift/2, width/2, depthAno], marker=4, area=0.01)
plc = mt.mergePLC([world, ano])

# pg.show(plc)
# pg.show(world)

# plc.exportVTK('plc.vtk')
maxA=np.array(ano.boundaryCenters()).max(axis=0)
minA=np.array(ano.boundaryCenters()).min(axis=0)

#%%
n_sensors_surf = 16*16
sensors_surf = np.zeros((n_sensors_surf, 3))
# area = [-10.1,10.1,-10.1,10.1]

area = [minA[0]-widthAno*2,maxA[0]+widthAno*2,minA[0]-widthAno*2,maxA[0]+widthAno*2]
# area = [2.5*i for i in area]

shape=(32,32)
xy= np.vstack(MyMS.regular(area,shape))
id_sensors = np.arange(0,len(xy[0,:]),4)
sensors_surf=xy[0:2,id_sensors].T

if add_bor==1:
    n_sensors_b=50
    nb_bor = 5
    sensors_b = np.zeros((n_sensors_b, 3))
    sensors_b[:, 0] = np.repeat(np.linspace(-10.1, 10.1,nb_bor),nb_bor*2)
    sensors_b[:, 1] = 0
    sensors_b[::, 2] = np.reshape([np.linspace(-10,-1,nb_bor*2),]*nb_bor,[1,n_sensors_b])
    
    n_sensors= n_sensors_surf + n_sensors_b
    sensors = np.zeros((n_sensors, 3))
    sensors= np.vstack([sensors_surf,sensors_b])
else:
    sensors = np.zeros((n_sensors_surf, 3))
    sensors[:,0:2]= sensors_surf
    
# or read elecs from file
# s99 = np.genfromtxt('elecs99.dat',delimiter=',')
# # place A,B and M (remote) first
# elecs = np.vstack([s99[-1,:],s99[-3:-1,:],s99[:-3,:]])

# %% add MALM electrodes 


EAp= (minA+maxA)/2
EAp[2] = depthAno
# EAp = [0.01,0.01,0]

# EBp= [world.xmin()/1.2, world.ymin()/1.2, 0]
# ENp= [world.xmax()/1.2, world.ymax()/1.2, 0]

ENp= [world.xmax()/1.2, world.ymax()/1.2, 0]
# EBp = [i * -1 for i in ENp]
EBp = [world.xmax()/3.5, world.ymax()/3.5, 0]


# elecs= np.vstack([EAp,EBp,ENp,sensors])
elecs= np.vstack([sensors,EAp,EBp,ENp])

# %% create malmScheme
# EA = 0
# EB = 1
# EN = 2

EA = len(elecs) -3
EB = len(elecs) -2
EN = len(elecs) -1

malmScheme = pg.DataContainerERT()
for ei in elecs:
    malmScheme.createSensor([ei[0], ei[1], ei[2]])

meas = []
for j in np.arange(0,len(elecs)-3):
    M=j
    meas.append([EA, EB,M, EN])

meas = np.array(meas)
malmScheme.resize(len(meas))    

for i, j in enumerate("abmn"):
    malmScheme.set(j, meas[:, i])
    
malmScheme.save('malmScheme.shm')
print(malmScheme)

# %% put the electrodes as node into the mesh
for pos in malmScheme.sensorPositions():
    plc.createNode(pos, marker=-99)  # electrode
    plc.createNode(pos + [0, 0, -0.01])  # refinement

# pos.array()
# %% refine again close to injection electrode (?) 

# %% create node along p1,p2 line

shape = (100,100)
xnew_min = world.xmin()
xnew_max = world.xmax()
ynew_min = xnew_min
ynew_max = xnew_max

# xnew_min = -20
# xnew_max = 20
# ynew_min = xnew_min
# ynew_max = xnew_max

xnew, ynew = gridder.regular((xnew_min, xnew_max, ynew_min, ynew_max), shape)
len(xnew)
grid_to_calc = np.zeros([len(xnew),3])
grid_to_calc[:,0] = xnew
grid_to_calc[:,1] = ynew

# xnew_min = world.xmin()
# xnew_max = world.xmax()
# xlin = np.linspace(xnew_min,xnew_max,300)

# for xy in grid_to_calc:
#     plc.createNode([xy[0],xy[1],0])  # electrode
#     plc.createNode([xy[0],xy[1],-0.01])  # refinement
#     # plc.createNode([xy[0],xy[1],depthAno/2])  # electrode
#     # plc.createNode([xy[0],xy[1],depthAno*2])  # electrode

for xyz in grid_to_calc:
    plc.createNode(pg.RVector3(xyz))
    # plc.createNode(pg.RVector3(xyz) + [0,0,-0.1])  # refinement

# pg.RVector3(xyz).array()

# for pos in malmScheme.sensorPositions():
#     plc.createNode(pos, marker=-99)  # electrode
#     plc.createNode(pos + [0, 0, -0.01])  # refinement


# %% create the mesh
mesh = mt.createMesh(plc, quality=1.3)
# mesh.createSecondaryNodes(1)

# pg.show(plc)

# print(mesh)
# pg.show(mesh)
# mesh.exportVTK('mesharea.vtk')

# mesh = mesh.createP2()
# mesh.exportVTK('meshfwd.vtk')

# mesh = pg.load('./mesh/mesh_Anocube_test.msh')

# %% compute geometric factors using homogeneous resistivity
# hom = simulateERT(mesh, res=100.0, scheme=malmScheme,useBert=True, calcOnly=True)
# malmScheme['k'] = 1. / hom['u']
#  # noiseLevel=noiseLevel,                    noiseAbs=noiseAbs, seed=1337)
# # max(np.array(hom['u']))
# # min(np.array(hom['u']))

# plt.figure()
# plt.scatter(elecs[:-3,0], elecs[:-3,1],c=np.array(hom['u']), cmap='viridis')
# plt.colorbar()
# plt.show()

# %% simluate for a given rhomap

u_elecs = simulateERT(mesh, scheme=malmScheme, res= [[1,SoilR],[4,AnoR]],
                        calcOnly=False, verbose=True, sr=False,
                        noiseLevel=noiseLevel, noiseAbs=noiseAbs, seed=1337)

u_elecs['u'] = u_elecs['rhoa']/u_elecs['k']
plt.figure()
plt.scatter(elecs[:-3,0], elecs[:-3,1],c=np.array(u_elecs['u']),
            cmap='viridis')
plt.colorbar()
plt.show()

# %% simluate for a given rhomap

# u_elecs_nonoise = simulateERT(mesh, scheme=malmScheme, res= [[1,SoilR],[4,AnoR]],
#                         sr=False, calcOnly=True, verbose=True)
# plt.figure()
# plt.scatter(elecs[:-3,0], elecs[:-3,1],c=np.array(u_elecs_nonoise['u']), 
#             cmap='viridis', vmin=5, vmax=10)
# plt.colorbar()
# plt.show()


#%% simulate for a given rhomap and return fields

# U = simulateERT(mesh, scheme=malmScheme, res=[[1,SoilR],[2,AnoR]],
#                       calcOnly=True, verbose=True, returnFields=True)

#%% 
Rmap = pg.solver.parseMapToCellArray([[1,SoilR],[4,AnoR]], mesh)
U = simulateERT(mesh=mesh, res=Rmap, scheme=malmScheme,returnFields=True,useBert=True, calcOnly=False,
                noiseLevel=noiseLevel, noiseAbs=noiseAbs, seed=1337)

#%% 
uT=U[EA]-U[EB]       
uA=U[EA]

np.array(U[EA]) - np.array(U[EB])
np.array(uT)

mesh.addData("rhomap", Rmap)
mesh.addData("uA", uA)
mesh.addData("uT", uT)
mesh.exportVTK('u.vtk')


# pg.show(mesh)


# gridpv = pgMesh2pvMesh(mesh,data=uA, label="uA")
# slices = gridpv.slice_orthogonal()
# slices.plot(notebook=False)


# p = pyvista.Plotter(notebook=False)
# ax =  drawSlice(p,mesh, normal=[0,1,0], data=Rmap, label="Rmap")
# drawSensors(p, malmScheme.sensorPositions(), diam=1, color='yellow')
# p.view_xz()
# ax.show()


#%%  compute zero level field data  on a regular grid -----
# i = EA
# ui = U[i]
# ui = U[EA]


uAz0_grid = []
uTz0_grid = []
for ni, node in enumerate(grid_to_calc):
    nearest = mesh.findNearestNode([node[0], node[1], node[2]])
    uAz0_grid.append(uA[nearest])                      
    uTz0_grid.append(uT[nearest])                      



plt.xlim([area[0],area[1]])
plt.ylim([area[0],area[1]])


p1 = [min(xnew),(area[0]+area[1])/2]
p2 = [max(xnew),(area[0]+area[1])/2]
pEXP.plot_line(xnew, ynew, uAz0_grid,p1,p2, interp=True)

plt.figure()
plt.subplot(1,2,1)
plt.scatter(xnew, ynew,c=uAz0_grid, cmap='viridis')
plt.colorbar()
plt.scatter(elecs[i,0],elecs[i,1])
plt.axis('square')
plt.subplot(1,2,2)
plt.scatter(xnew, ynew,c=uTz0_grid, cmap='viridis')
plt.colorbar()
plt.scatter(elecs[i,0],elecs[i,1],c='black')
plt.axis('square')


# Y=np.random.normal
# uTz0_grid = np.array(uTz0_grid) + max(uTz0_grid)*noiseLevel*np.random.randn(len(uTz0_grid))
# pEXP.plot_line(xnew, ynew, uTz0_grid,p1,p2, interp=True)


#%% Add noise
# uAz0_grid = np.array(uAz0_grid) + max(uAz0_grid)*noiseLevel*np.random.randn(len(uAz0_grid))
# pEXP.plot_line(xnew, ynew, uAz0_grid,p1,p2, interp=True)

uAz0_grid = np.vstack([grid_to_calc.T,uAz0_grid]).T
uTz0_grid = np.vstack([grid_to_calc.T,uTz0_grid]).T


#%%  save field data -----

# dstruct = { "SoilR" : SoilR, "AnoR" : AnoR , "HWD" : [HAno,widthAno,depthAno], "XYU" : uz0_grid,
#            "shape" : shape, "p12": [p1,p2], "Fieldu": np.array(U), "uA": uA, "uT": uT, 
#            "elecs": elecs, "EA_EB_EN": [EAp, EBp, ENp],
#            "mesh": mesh}

# dstruct = { "SoilR" : SoilR, "AnoR" : AnoR , "HWD" : [HAno,widthAno,depthAno], "XYU" : uz0_grid,
#            "shape" : shape, "p12": [p1,p2], 
#            "uTz0_grid": uTz0_grid,
#            "uA": uA, "uT": uT.array(), "u_elecs": hom['u'].array(),
#            "elecs": elecs, "EA_EB_EN": [EAp, EBp, ENp]}

dstruct = { "SoilR" : SoilR, "AnoR" : AnoR , "HWD" : [HAno,widthAno,depthAno], 
            "uAz0_grid" : uAz0_grid,
            "uTz0_grid" : uTz0_grid,
            "u_elecs": u_elecs['u'].array(),
            "elecs": elecs, 
            "EA_EB_EN": [EAp, EBp, ENp],
            "shape" : shape, 
            "p12": [p1,p2]}

# dstruct = { "SoilR" : SoilR, "AnoR" : AnoR , "HWDLs" : [HAno,widthAno,depthAno,lengthAno,shift], "XYU" : uz0_grid,
#            "shape" : shape, "p12": [p1,p2],
#            "noiseLevel": noiseLevel}

afile = open(SimName + '.pkl', 'wb')
pickle.dump(dstruct, afile)
afile.close()


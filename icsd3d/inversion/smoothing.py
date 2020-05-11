# -*- coding: utf-8 -*-
"""
Created on Mon May 11 17:29:01 2020
@author: Benjamin
"""

import numpy as np
from scipy.sparse import diags

#%%
def nx_ny(coord):
    """find number of nodes in each direction, has to be a regular grid"""
    nx = np.unique(np.round(coord[:, 0], 3)).shape[0]
    ny = np.unique(np.round(coord[:, 1], 3)).shape[0]
    
    return nx,ny
    
def nx_ny_nz(coord):
    """find number of nodes in each direction, has to be a regular grid"""
    nx = np.unique(np.round(coord[:, 0], 3)).shape[0]
    ny = np.unique(np.round(coord[:, 1], 3)).shape[0]
    nz = np.unique(np.round(coord[:, 2], 3)).shape[0]
    
    return nx,ny,nz

#%% Smoothing methods for different mesh types

#%% 2d

def regularize_A_x_y(coord,alphaSx,alphaSy):
    """create and append rows for spatial regularization to A, second derivative is applied in both direction x and y
    math:: Dx = ?? Dy=
    We used a ponderate diagonal matrix with coeffcient 1,-2, 1
    """
    nx,ny= nx_ny(coord)

    ncells = nx*ny;       
    Dx=diags([1, -2, 1], [-1, 0, 1], shape=(ncells, ncells)).todense()
    idx=np.arange(0,len(Dx),nx)
    Dx[idx,:] = 0
    idx2=np.arange(nx,len(Dx),nx)
    Dx[idx2,:] = 0
    
    Dy=diags([1, -2, 1], [-nx, 0, nx], shape=(ncells, ncells)).todense()
    idy=np.arange(0,nx)
    Dy[idy,:] = 0
    idy2=np.arange(((ny-1)*nx+1),(nx)*(ny))
    Dy[idy2,:] = 0
    
    reg_Ax = alphaSx*np.array(Dx.transpose()*Dx)
    reg_Ay = alphaSy*np.array(Dy.transpose()*Dy)
    
    return reg_Ax, reg_Ay
    
#%% 3d
def regularize_A_3d(nVRTe,coord):
    """create and append rows for spatial regularization to A"""
    nx,ny,nz= nx_ny_nz(coord)
    reg = []
    vrte = range(1, nVRTe + 1)
    vrte = np.reshape(vrte,(ny, nx, nz))
    for z in range(nz):
	        for y in range(ny):
	            for x in range(nx):
	                minus = vrte[y, x]
	                if x + 1 in range(nx):
	                    plus = vrte[y , x + 1]
	                    row = np.zeros(nVRTe, int)
	                    row[minus -1] = - 1
	                    row[plus -1] = + 1
	                    reg.append(row)
	                    if y + 1 in range(ny):
	                        plus = vrte[y + 1, x]
	                        row = np.zeros(nVRTe)
	                        row[minus -1] = - 1
	                        row[plus -1] = + 1
	                        reg.append(row)
    reg_A = np.array(reg)
    
    return reg_A
        
        
def regularize_A_UnstructuredMesh3d(coord,nVRTe,k_neighbors=4): # should also work for the 2d case

    reg = []
    for VRTEnb in range(nVRTe):
        dist =  np.linalg.norm(coord[VRTEnb]-coord, axis=1)
        closest = np.argsort(dist)
        k = k_neighbors  # For each point, find the k closest current sources
        Ind = closest[1:k+1]
        row = np.zeros(nVRTe) # add a line to the regularisation A with k non-null coefficients
        knorm = dist[closest[1:k+1]]/dist[closest[1:k+1]].sum(axis=0,keepdims=1)
        row[Ind]= -knorm
        row[VRTEnb]= 1 # = one for the actual current source
        reg.append(row)
        test=[1]
        mask = np.in1d(test, VRTEnb)
        if mask.any()==True: 
            print('nll')
            # self.fc = plt.figure('TEST regularisation')
            # ax = self.fc.add_subplot(111, projection='3d')
            # ax.scatter(self.coord[VRTEnb,0], self.coord[VRTEnb,1], self.coord[VRTEnb,2], linewidths=12,
            #            facecolor = 'green', edgecolor = 'green')
            # ax.scatter(self.coord[Ind,0], self.coord[Ind,1], self.coord[Ind,2], linewidths=12,
            #            facecolor = 'red', edgecolor = 'red')
            # ax.set_xlim([min(self.coord_x),max(self.coord_x)])
            # ax.set_ylim([min(self.coord_y),max(self.coord_y)])
            # ax.set_zlim([min(self.coord_z),max(self.coord_z)])
            # self.fc.savefig(self.path2save+ 'TEST regularisation', dpi = 600)
            # #plt.show()
            # #plt.close()
        reg_A = np.array(reg)
    return reg_A
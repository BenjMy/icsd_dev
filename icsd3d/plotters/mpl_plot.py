# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:44:26 2020

@author: Benjamin
"""

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
import numpy as np
from scipy.interpolate import griddata

from importers.read import *
from gridder.mkgrid import mkGrid_XI_YI 

       
def _fig_Interpolation_(coord, data, **kwargs):
    """ plot the interpolation of the VRTe current fractions """
    
    coord_x, coord_y = parseCoord(coord,dim='2d')
    XI, YI= mkGrid_XI_YI(coord_x,coord_y)
    points = np.column_stack((coord_x, coord_y))
    
    grid = griddata(points, data, (XI, YI), method = 'linear') # grid is a np array
    plt.imshow(grid,
               extent = (min (coord_x), max(coord_x), min(coord_y), max(coord_y)),
               aspect = 'auto', origin = 'lower', cmap= 'jet')
    cbar = plt.colorbar()
    
    if kwargs.get('clim') is not None:
        plt.clim(clim[0],clim[1])
    
    if kwargs.get('lgd_label') is not None:
        cbar.set_label(kwargs.get('lgd_label'), labelpad = 10)
    else:
        cbar.set_label('Fraction of Current Source', labelpad = 10)

def _fig_RealSources_(sc):
    """ add known point sources if present """
    if sc == None:
        return
    for s in sc:
        sx = float(s.split(',')[0])
        sy = float(s.split(',')[1])
        plt.plot(sx, sy,'ow', markersize = 10, markeredgecolor = 'k')

def _fig_ReturnElec_(retElec):
    """ plot the return electrode """
    if retElec == None:
        return
    print('reading return electrode: ', retElec)
    retElecx = float(retElec.split(',')[0])
    retElecy = float(retElec.split(',')[1])
    plt.plot(retElecx, retElecy,'sw', markersize = 10)

def _fig_VRTe_(coord,data_sol):
    """ plot the VRTe current franctions """
    coord_x, coord_y = parseCoord(coord,dim='2d')
    norm_z = (data_sol - min(data_sol)) / (max(data_sol) - min(data_sol))
    grey_cm = plt.cm.get_cmap('Greys')
    edgecolor_norm_z = grey_cm(norm_z)
    jet_cm = plt.cm.get_cmap('jet')
    facecolor_norm_z = jet_cm(norm_z)
    plt.scatter(coord_x, coord_y, facecolor = facecolor_norm_z, edgecolor = edgecolor_norm_z, cmap = 'jet')

            
def _fig_Axis_Labels_(title):
    plt.title(title)
    plt.ylabel('y [m]',fontsize=12)
    plt.xlabel('x [m]',fontsize=12)
    axes = plt.gca()
    # axes.set_xlim([0,0.53])
    # axes.set_ylim([0,0.52])
    plt.tick_params(axis='both', which='major')
    plt.tight_layout()
    axes.set_aspect('equal')
    
def parseCoord(coord,dim):

    coord_x = coord[:,0]
    coord_y = coord[:,1]
    
    if dim == '3d':
        coord_z = coord[:,2]
        return coord_x, coord_y, coord_z
    else:
        return coord_x, coord_y
    
def plotRemotes(path,dim,pltRemotes=False):
    
    if dim == '2d':
        print('not yet implemented')
        return 
    else:
        if pltRemotes == True:
            RemLineNb, Injection, coordE, pointsE= load_geom(path) # geometry file containing electrodes position includinf remotes 
            ax.scatter(coordE[RemLineNb:RemLineNb+2,1], coordE[RemLineNb:RemLineNb+2,2], coordE[RemLineNb:RemLineNb+2,3],
                        marker="v", color='black',s=60, label = 'Remotes')
            ax.scatter(coordE[Injection,1], coordE[Injection,2], coordE[Injection,3],
                        marker="*", color='black',s=60, label = 'A')
            ax.scatter(coordE[:RemLineNb,1], coordE[:RemLineNb,2], coordE[:RemLineNb,3],
                        marker="s", color='black',s=60, label = 'Velecs')
            ax.view_init(azim=-101, elev=35)
    
#%% Specific plot functions for ICSD outputs

def plotPareto(wr,pareto_list_FitRes,pareto_list_RegRes,IdPtkneew,path):
    p, ax = plt.subplots()
    # ax.annotate('Wr=' + str(int(wr)), xy=(float(np.asarray(pareto_list_FitRes)[IdPtkneew]), 
    #                                  float(np.asarray(pareto_list_RegRes)[IdPtkneew])), 
    #                               xytext=(float(np.asarray(pareto_list_FitRes)[IdPtkneew])+max(pareto_list_FitRes)/3, 
    #                                  float(np.asarray(pareto_list_RegRes)[IdPtkneew])+max(pareto_list_RegRes)/3),
    #                               arrowprops=dict(facecolor='black', shrink=0.05))
    plt.plot(float(np.asarray(pareto_list_FitRes)[IdPtkneew]), float(np.asarray(pareto_list_RegRes)[IdPtkneew]), 'og')
    plt.plot(pareto_list_FitRes, pareto_list_RegRes, 'or')
    ax = plt.gca()
    ax.tick_params(axis = 'both', which = 'both', direction = 'out')
    ax.grid()
    plt.xlabel('Residual')
    plt.ylabel('Roughness')
    plt.tight_layout()
    plt.show()
    plt.savefig(path+'ParetoFront.png', dpi = 600)
    
    return p, ax

def plot_knee_icsd(wr,kn):
    """ Plot CSD for the best regularisation parameter after L-curve automatic analysis using a knee-locator
    
    Parameters
    ------------
    self
    """
    KneeWr=wr
    kn.plot_knee_normalized()
    plt.show(kn)
    run_single()


def plotFIT(b,b_w,xfun,path):
    stopat=len(b) # 204
    plt.figure()
    plt.subplot(121)
    plt.plot(xfun[:stopat] + b_w[:stopat], 'or', label = 'Inverted CSD')
    plt.plot(b_w[:stopat], 'ob', label = 'True model')
    plt.xlabel('Measurement number')
    plt.ylabel('R [Ohm]')
    plt.legend()
    plt.subplot(122)
    plt.plot(xfun[:stopat] + b_w[:stopat], b_w[:stopat], 'or')
    plt.xlabel('Inverted CSD, R [Ohm]')
    plt.ylabel('True model, R [Ohm]')
    plt.tight_layout()
    plt.savefig(path+'Fit.png', dpi = 600)
    plt.show()

    
def plotCSD2d(coord,data_sol,b,b_w,xfun,path,pareto,retElec=None, sc=None, **kwargs):
    """ Plot CSD in 2d, using matplotlib and scipy interpolation
    
    Parameters
    ------------
    self
    """
    f = plt.figure('surface')
    _fig_Interpolation_(coord,data_sol)
    _fig_VRTe_(coord,data_sol)
    _fig_RealSources_(sc)
    _fig_ReturnElec_(retElec)
    
    if kwargs.get('title_wr') is not None:
        title=r'$\lambda$=' + str(kwargs.get('title_wr') )
        _fig_Axis_Labels_(title)

    return f

    if not pareto:
        plotFIT(b,b_w,xfun,path)

            
def plotCSD3d(wr,coord,data,path,filename,knee,KneeWr,title=None,pltRemotes=False,**kwargs):
    """ plot scattered 3d current sources density for a given regularisation weight wr 
    (can be the knee location if pareto-curve mode is run)
    
    Parameters
    ----------
    sc: sources coordinates
    
    kwargs (to add) : 'sc' (plot source position (for a synthetic case experiment)

    """
    coord_x, coord_y, coord_z = parseCoord(coord,dim='3d')

    f = plt.figure('volume')

    step=(max(coord_x)-min(coord_x))/10

    xlin=np.arange(min(coord_x),max(coord_x),step)
    ylin=np.arange(min(coord_y),max(coord_y),step)
    zlin=np.arange(min(coord_z),max(coord_z),step)
    #generate new grid
    X,Y,Z=np.meshgrid(xlin,ylin,zlin)
    
    data_2_plot= data  
    ax=f.gca(projection='3d')
    sc=ax.scatter(coord_x, coord_y, coord_z, c=data_2_plot, cmap ='coolwarm', s=data_2_plot*1e4,
              )
    # if self.clim is None:
    #     print('none clim')
    # sc.set_clim(self.clim)
    cbar = plt.colorbar(sc)
    # self.labels()
    # cbar.set_label(self.physLabel)
    # ax.set_zlim([-10,0])
    
    # label=kwargs.get('label_nor', 'normal')
            
    for key, value in kwargs.items():
        print("{0} = {1}".format(key, value))
        
        if key == 'zlim':
            ax.set_zlim(
                [kwargs.get(key)[0],
                  kwargs.get(key)[1]]
                )
    
    plotRemotes(path,dim='3d',pltRemotes=False) # plot remotes and injection position
            
    if title==None:
        title= 'Scattered current sources density, wr=' + str(wr) 
    else:
        title= title + ', wr=' + str(wr)
    plt.legend()
    plt.title(title)
    plt.savefig(path +  filename + '_icsd_scatter.png' , dpi=550,bbox_inches='tight',pad_inches = 0)

    if knee==True:
        if wr==KneeWr:
            plt.savefig(path+ filename + 'icsd_knee_scatter'+ str(KneeWr) + '.png',dpi=550,bbox_inches='tight',pad_inches = 0)
    
    plt.show()
    
    return f
       
        
#%% Generic plot functions


def scatter2d(coord, data, label, path, filename, pltRemotes=False, **kwargs):
      
    coord_x, coord_y = parseCoord(coord,dim='2d')

    f = plt.figure('volume')
    step=(max(coord_x)-min(coord_x))/10
    xlin=np.arange(min(coord_x),max(coord_x),step)
    ylin=np.arange(min(coord_y),max(coord_y),step)

    X,Y=np.meshgrid(xlin,ylin)
    ax=f.gca()
    sc=ax.scatter(coord_x, coord_y, c=data, cmap ='coolwarm')
    cbar = plt.colorbar(sc)
    cbar.set_label(label)      
    ax.set_ylabel('y [m]',fontsize=15)
    ax.set_xlabel('x [m]',fontsize=15)
    
    plotRemotes(path,dim='2d',pltRemotes=False) # plot remotes and injection position

    plt.legend()
    # plt.title(title)
    plt.savefig(path+  filename + '_icsd_scatter.png' , dpi=550,bbox_inches='tight',pad_inches = 0)
    plt.show()
    
    return f
        

def scatter3d(coord, data, label, path, filename, pltRemotes=False, **kwargs):
      
    coord_x, coord_y, coord_z = parseCoord(coord,dim='3d')

    
    f = plt.figure('volume')
    step=(max(coord_x)-min(coord_x))/10
    xlin=np.arange(min(coord_x),max(coord_x),step)
    ylin=np.arange(min(coord_y),max(coord_y),step)
    zlin=np.arange(min(coord_z),max(coord_z),step)
    X,Y,Z=np.meshgrid(xlin,ylin,zlin)
    ax=f.gca(projection='3d')
    sc=ax.scatter(coord_x, coord_y, coord_z, c=data, cmap ='coolwarm')
    cbar = plt.colorbar(sc)
    cbar.set_label(label)      
    ax.set_ylabel('y [m]',fontsize=15)
    ax.set_xlabel('x [m]',fontsize=15)
    
    plotRemotes(path,dim='3d',pltRemotes=False) # plot remotes and injection position

    plt.legend()
    # plt.title(title)
    plt.savefig(path+  filename + '_icsd_scatter.png' , dpi=550,bbox_inches='tight',pad_inches = 0)
    plt.show()
    
    return f


def labels(method):
    """ Parse graphical labels to plot
    """
    if method=='F1':      
        physLabel= 'normed misfit F1'
    if method=='Pearson':      
        physLabel= 'Pearson r coeff' 
    
    return physLabel
            
            
def plotContour2d(coord,data_sol,physLabel,path,retElec=None, sc=None, **kwargs):
    """ Plot contour in 2d, using matplotlib and scipy interpolation
    
    Parameters
    ------------
    self
    """
    f = plt.figure('surface')
    _fig_Interpolation_(coord,data_sol,lgd_label=physLabel)
    _fig_VRTe_(coord,data_sol)
    _fig_RealSources_(sc)
    _fig_ReturnElec_(retElec)

    if kwargs.get('title_wr') is not None:
        title=r'$\lambda$=' + str(title_wr)
        _fig_Axis_Labels_(title)

    return f


def showObs2d(path, **kwargs):
    """ Plot contour in 2d, using matplotlib and scipy interpolation. Required surface and borehole electrode to make the 2d interpolation possible
    
    Parameters
    ------------
    self
    """
    filename='ObsData.txt'
    
    if kwargs.get('filename') is not None:
        filename = kwargs.get('filename')
        
    RemLineNb, Injection, coordE, pointsE= load_geom(path) # geometry file containing electrodes position includinf remotes 
    # print(path)
    # print(path)

    data_obs = load_obs(path,filename)
    f = plt.figure('surface')
    _fig_Interpolation_(pointsE,data_obs,lgd_label='U/I')
    _fig_VRTe_(pointsE,data_obs)
    _fig_RealSources_(sc=None)
    _fig_ReturnElec_(retElec=None)
    _fig_Axis_Labels_(title='Observations')
    plt.show(f)
    plt.close(f)

    return f

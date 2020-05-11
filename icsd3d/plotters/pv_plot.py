# -*- coding: utf-8 -*-
"""
Created on Mon May 11 14:44:26 2020
@author: Benjamin
3D plots using pyvista
"""

import pyvista as pv
import os
import matplotlib.pyplot as plt
import numpy as np


def plotCSD3d_pyvista(coords, path=None, filename_root='Solution.dat', **kwargs): # add kwards or *args for pyvista plot

    if path is None:
        cwd = os.getcwd()
        path= cwd
        
    filename = path + filename_root
    
    data_2_plot =np.genfromtxt(filename)
    coord_x= data_2_plot[:,0]
    coord_y= data_2_plot[:,1]
    coord_z= data_2_plot[:,2]
    coord= data_2_plot[:,:-1]

    opacity = [0, 0, 0.1, 0.3, 0.6, 0.9, 1]
    grid = pv.UniformGrid()
    spc=(max(coord_x)-min(coord_x))/10
    xdim = int(round((max(coord_x)-min(coord_x))/spc))
    ydim = int(round((max(coord_y)-min(coord_y))/spc))
    zdim = int(round((max(coord_z)-min(coord_z))/spc))
    grid.dimensions = (xdim, ydim, zdim)
    grid.dimensions = np.array(grid.dimensions) +1
    grid.origin = (min(coord_x), min(coord_y), min(coord_z)) # The bottom left corner of the data set
    grid.spacing = (spc, spc,spc) # These are the cell sizes along each axis
    

    pv.set_plot_theme('document')
    poly = pv.PolyData(coord)
    pvfig = pv.Plotter(notebook=False,window_size=[600, 600])
    
    # if self.mesh!=None:
    #     print('plot mesh.vtk')
    #     ModelVtk = pv.read(self.path2load + self.mesh)
    #     cmap = plt.cm.get_cmap('viridis', 2)
    #     pvfig.add_bounding_box()
    #     pvfig.add_mesh(cmap=cmap,mesh=ModelVtk,scalars='rhomap', opacity=0.2)    # add a dataset to the scene
    # 
    pvfig.add_mesh(poly, point_size=15.0, scalars=data_2_plot[:,3], opacity=opacity, render_points_as_spheres=True,cmap='jet')
    print('interpolation spacing='+ str(spc))
    interpolated = grid.interpolate(poly, radius=spc*2)
    cmap = plt.cm.get_cmap('jet',10)
    contours = interpolated.contour()
    # pvfig.add_mesh(interpolated, show_scalar_bar=False, cmap=cmap,opacity=0.3)
    pvfig.add_mesh(contours, show_scalar_bar=False, opacity= opacity,cmap='jet')
    pvfig.show_bounds(bounds=[min(coord_x), max(coord_x), 
                                    min(coord_y),max(coord_y),
                                    min(coord_z), 0],font_size=16)
    pvfig.add_axes()
    pvfig.show_axes()
    pvfig.add_scalar_bar('Normalized Current density',width=0.25,vertical=False,position_x=0.3)

    # if self.plotElecs==True:
    #     pvfig.add_points(self.pointsE)
    #     # p.add_point_labels(pointsE,coordE[:RemLineNb,0].astype(int), point_size=15, font_size=35,
    #     #     shape_opacity=0.01, margin=4.)
    #     pvfig.add_point_labels(self.pointsE,self.coordE[:self.RemLineNb,0].astype(int), point_size=int(spc)-5, font_size=int(spc)-5,
    #         shape_opacity=0.01, margin=4.)
    # # set_focus(self.pointsE)
    # # pvfig.set_focus(point=self.pointsE[0])
    # # pvfig.set_scale(xscale=1, yscale=1, zscale=0.2, reset_camera=True)
    pvfig.show(auto_close=True)  


    # if self.knee==True:
    #    if self.wr==self.KneeWr:
    #        # pvfig.screenshot('Pts_iCSD_knee'+ str(self.ObsName) + '.png')
    #         pvfig.screenshot(self.path2save+ 'Pts_iCSD_knee_wr'+ self.obs +  str(self.KneeWr) + '.png')
    #         if self.gif3d==True:
    #             viewup = [0.5, 0.5, 1]
    #             path = pvfig.generate_orbital_path(factor=2.0, shift=poly.length, viewup=viewup, n_points=36)
    #             # p.open_gif(path2file + simName + "orbit.gif")
    #             pvfig.open_movie(self.path2save +  self.obs + "orbit.gif", framerate=4)
    #             pvfig.orbit_on_path(path, write_frames=True, viewup=[0, 0, 1])
    # else:
    #         pvfig.screenshot(self.path2save+  self.obs + 'Pts_iCSD_wr'+ str(self.wr) + '.png')
    #         if self.gif3d==True:
    #             viewup = [0.5, 0.5, 1]
    #             path = pvfig.generate_orbital_path(factor=2.0, shift=poly.length, viewup=viewup, n_points=36)
    #             # p.open_gif(path2file + simName + "orbit.gif")
    #             pvfig.open_movie(self.path2save + self.obs +"orbit.gif", framerate=4)
    #             pvfig.orbit_on_path(path, write_frames=True, viewup=[0, 0, 1])
                
    # pvfig.close()
    
 
    
    


    # xyB4321 = np.loadtxt("digitized_B4321.dat",skiprows=0,delimiter=',')
    # x=xyB4321[:,0]
    # y=xyB4321[:,1]
    
    # zo = np.linspace(-5, 0, num=len(y))
    # points = np.c_[x, y, zo]
    # spline = pv.Spline(points, 15)
    # spline
    
    # slc = interpolated.slice_along_line(spline)       
    # # p = pv.Plotter()
    # pvfig.add_mesh(slc, cmap=cmap)
    # pvfig.add_mesh(interpolated.outline())
    # # pvfig.show(cpos=[1, -1, 1])
    # pv.set_plot_theme('document')
    # single_slice = interpolated.slice(normal=[0, 1, 0])
    # p = pv.Plotter(notebook=False)
    # p.add_mesh(interpolated.outline(), color="k")
    # p.add_mesh(single_slice, cmap=cmap)
    # p.view_xz()
    # p.add_axes()
    # p.show_axes()
    # p.show_grid()
    # p.show_bounds(font_size=16)      
    # p.show()
    
    


def plotCSD3d_pyvistaSlice(coords, path=None, filename_root='Solution.dat', **kwargs): # add kwards or *args for pyvista plot

    if path is None:
        cwd = os.getcwd()
        path= cwd
        
    filename = path + filename_root
    
    data_2_plot =np.genfromtxt(filename)
    coord_x= data_2_plot[:,0]
    coord_y= data_2_plot[:,1]
    coord_z= data_2_plot[:,2]
    coord= data_2_plot[:,:-1]

    opacity = [0, 0, 0.1, 0.3, 0.6, 0.9, 1]
    grid = pv.UniformGrid()
    spc=(max(coord_x)-min(coord_x))/10
    xdim = int(round((max(coord_x)-min(coord_x))/spc))
    ydim = int(round((max(coord_y)-min(coord_y))/spc))
    zdim = int(round((max(coord_z)-min(coord_z))/spc))
    grid.dimensions = (xdim, ydim, zdim)
    grid.dimensions = np.array(grid.dimensions) +1
    grid.origin = (min(coord_x), min(coord_y), min(coord_z)) # The bottom left corner of the data set
    grid.spacing = (spc, spc,spc) # These are the cell sizes along each axis
    

    pv.set_plot_theme('document')
    poly = pv.PolyData(coord)
    print('interpolation spacing='+ str(spc))
    interpolated = grid.interpolate(poly, radius=spc*2)
    cmap = plt.cm.get_cmap('jet',10)
    contours = interpolated.contour()

    single_slice = interpolated.slice(normal=[0, 1, 0])
    p = pv.Plotter(notebook=False)
    p.add_mesh(interpolated.outline(), color="k")
    p.add_mesh(single_slice, cmap=cmap)
    p.view_xz()
    p.add_axes()
    p.show_axes()
    p.show_grid()
    p.show_bounds(font_size=16)      
    p.show(auto_close=True)  

    
    
    
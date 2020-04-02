import re
import argparse
import numpy as np
import matplotlib.pyplot as plt
import matplotlib
from scipy.interpolate import griddata
from scipy.linalg import lu
from scipy.optimize import lsq_linear, curve_fit, least_squares, leastsq
from matplotlib.backends.backend_pdf import PdfPages
import matplotlib.font_manager
from mpl_toolkits.mplot3d import Axes3D

#font = {'family' : 'normal',
#        'weight' : 'normal',
#        'size'   : 12}
#matplotlib.rc('font', **font)
import pyvista as pv
import vtk
from numpy import sin, cos, pi
from skimage import measure
from scipy.interpolate import griddata as gd
from mpl_toolkits.mplot3d.art3d import Poly3DCollection
from kneed import KneeLocator
import os



class iCSD_class():
    def get_args(self):
        ap = argparse.ArgumentParser(description = ' current source density inversion')
        ap.add_argument('-sim', metavar = 'BERT_VRTe_calculation', type = str,
            help = 'file containing the forward calculation for each virtual source, going into the A matrix',
            default = 'VRTeSim.txt')
        ap.add_argument('-obs', metavar = 'Observed_MALM_data', type = str,
            help = 'measured MALM data, b column of the linear system',
            default = 'ObsData.txt')
        ap.add_argument('-coord_file', metavar = 'VRTe_coordinates', type = str,
            help = 'coordinates of the VRTe for plotting',
            default = 'VRTeCoord.txt')
        ap.add_argument('-wr', metavar = 'weight regularization', type = float,
            help = 'first derivative regularization weight', default = 25, nargs = '*')
        ap.add_argument('-wc', metavar = 'weight constrain', type = float,
            help = 'current conservation constrain, sum current fractions = 1', default = 10000)
        ap.add_argument('-sc', metavar = 'source coords', type = str,
            help = 'coordinates of the sources, format = x1,y1 x2,y2', nargs = '*')
        ap.add_argument('-retElec', metavar = 'Return Elec', type = str,
            help = 'coordinates of the return electrode, format = x1,y1') # rhi elec 64 = 0.44,0.8
        ap.add_argument('-pareto', action = 'store_true',
            help = 'if True run many icsd to explore the Pareto front', default = False)
        ap.add_argument('-errRmin', metavar = 'min R in err', type = float,
            help = 'min R to which the err is applied before passing to constant error', default = 1)
        ap.add_argument('-pareto_MinErr', metavar = 'min Pareto error', type = float,
            help = 'starting error for pareto mode', default = 0.0001)
        ap.add_argument('-pareto_MaxErr', metavar = 'max Pareto error', type = float,
            help = 'maximum error for pareto mode', default = 50)
        ap.add_argument('-pareto_nSteps', metavar = 'number steps Pareto', type = int,
            help = 'number of error steps for the pareto mode', default = 21)
        ap.add_argument('-obs_err', metavar = 'obs weight mode', type = str,
            help = 'choose between constant weight and w = 1/sqrt(abs(obs))', default = 'const',
            choices = ['const', 'sqrt'])
        ap.add_argument('-k', metavar = 'nb of neighbours', type = int,
            help = 'choose the number of neighbours sources for the regularisation', default = 2)
        ap.add_argument('-TL', action = 'store_true',
            help = 'if True run Time-Lapse icsd', default = False)
        ap.add_argument('-F1w_file', metavar = 'F1w_file', type = str,
            help = 'guided after F1 inversion',
            default = None)
        #ap.add_argument('-Thresh', metavar = 'min current threshold', type = tuple,
        #    help = 'min current threshold for volume visualisation', default = [0.03,0.1])
        self.args = ap.parse_args()


# miscellaneous functions 
    # def search_string_in_file(self):
    #     """Search for the given string in file and return lines containing that string,
    #     along with line numbers"""
    #     line_number = 0
    #     line_of_results = []
    #     # Open the file in read only mode
    #     with open(file_name, 'r') as read_obj:
    #         # Read all lines in the file one by one
    #         for line in read_obj:
    #             # For each line, check if line contains the string
    #             line_number += 1
    #             if string_to_search in line:
    #                 # If yes, then add the line number & line as a tuple in the list
    #                 line_of_results.append((line_number))
                            
    #     # Return list of tuples containing line numbers and lines where string is found
    #     return line_of_results

    ### PLOT 

    def _fig_VRTeColored_(self):
        norm_z = (self.x.x - min(self.x.x)) / (max(self.x.x) - min(self.x.x))
        #print('Norm z:' + str(norm_z))
        edgecolor_norm_z = plt.cm.jet(norm_z)
        facecolor_norm_z = plt.cm.jet(norm_z)
        #plt.scatter(self.coord_x, self.coord_y, facecolor = facecolor_norm_z, edgecolor = edgecolor_norm_z, cmap = 'jet')
        #ax = plt.gca()
        #print(self.coord_z)
        self.f = plt.figure()
        ax = self.f.add_subplot(111, projection='3d')
        xs,ys,zs = np.random.random(50),np.random.random(50),np.random.random(50)
        values = np.random.random(50)*10
        p = ax.scatter3D(self.coord_x,self.coord_y, self.coord_z, c=norm_z, 
            facecolor = facecolor_norm_z, edgecolor = edgecolor_norm_z, cmap = 'jet')
        self.f.colorbar(p, ax=ax)

    def _fig_Axis_Labels_(self):
        plt.ylabel('y [m]',fontsize=12)
        plt.xlabel('x [m]',fontsize=12)
        axes = plt.gca()
        axes.set_xlim([min(self.coord_x),max(self.coord_x)])
        axes.set_ylim([min(self.coord_y),max(self.coord_y)])
        axes.set_zlim([min(self.coord_z),max(self.coord_z)])
        plt.tick_params(axis='both', which='major')
        plt.tight_layout()

    def _plotCSD_(self):
        self._fig_VRTeColored_()
        self._fig_Axis_Labels_()
        if not self.args.pareto:
            self._plotFIT_()

    def _Interp3d_(self):
        print('interp with matplotlib')

        # def fun(x, y, z):
        #     return cos(x) + cos(y) + cos(z)

        # x, y, z = pi*np.mgrid[-1:1:31j, -1:1:31j, -1:1:31j]
        # vol = fun(x, y, z)
        # verts, faces, _, _ = measure.marching_cubes(vol, 0, spacing=(0.1, 0.1, 0.1))

        # fig = plt.figure()
        # ax = fig.add_subplot(111, projection='3d')
        # ax.plot_trisurf(verts[:, 0], verts[:,1], faces, verts[:, 2],
        #                 cmap='Spectral', lw=1)
        #plt.show()
        print(min(self.coord_z))
        print(max(self.coord_z))
        step=(max(self.coord_x)-min(self.coord_x))/10
        xlin=np.arange(min(self.coord_x),max(self.coord_x),step)
        ylin=np.arange(min(self.coord_y),max(self.coord_y),step)
        zlin=np.arange(min(self.coord_z),max(self.coord_z),step)

        #generate new grid
        X,Y,Z=np.meshgrid(xlin,ylin,zlin)

        #interpolate "data.v" on new grid "inter_mesh"
        V = gd((self.coord_x,self.coord_y,self.coord_z), self.x.x, (X,Y,Z), method='linear')

        fileName= 'ExportSolint.dat'
        ExpSOLint = np.vstack([X,Y,Z,V])
        ExpSOLint= ExpSOLint.T
        f = open(fileName,'w')
        np.savetxt(f, ExpSOLint, fmt='%1.2f %1.2f %1.2f %1.5f', delimiter='\t',header='X Y Z i')   # X is an array
        f.close()

        #Plot values
        fig = plt.figure()
        ax=fig.gca(projection='3d')
        sc=ax.scatter(X, Y, Z, V, cmap=plt.hot())
        #plt.colorbar(sc)
        plt.title('This is a test')
        plt.show()

        verts, faces, _, _ = measure.marching_cubes_lewiner(V, level=0.1, spacing=(0.1, 0.1, 0.1))
        fig = plt.figure(figsize=(10, 10))
        ax = fig.add_subplot(111, projection='3d')

        mesh = Poly3DCollection(verts[faces], alpha=0.1)
        face_color = np.array((255.0/255.0, 54.0/255.0, 57/255.0, 1.0))
        # Plot
        ax.add_collection3d(mesh)
        ax.set_xlim([min(self.coord_x),max(self.coord_x)])
        ax.set_ylim([min(self.coord_y),max(self.coord_y)])
        ax.set_zlim([min(self.coord_z),max(self.coord_z)])



    def _Export_sol_(self):
        fileName= 'ExportSol.dat'
        ExpSOL = np.vstack([self.coord_x,self.coord_y,self.coord_z,self.x.x])
        ExpSOL= ExpSOL.T
        print(ExpSOL)
        f = open(fileName,'w')
        np.savetxt(f, ExpSOL, fmt='%1.2f %1.2f %1.2f %1.5f', delimiter='\t',header='X Y Z i')   # X is an array
        f.close()



    def _Interp3dPyvista_(self):

        #pvfig = plt.figure('pareto', figsize=(10,4))
        pv.set_plot_theme('document')
        poly = pv.PolyData(self.coord)
        pvfig = pv.Plotter(notebook=False,window_size=[1000, 1000])
        # pvfig = pv.BackgroundPlotter()

        grid = pv.UniformGrid()
        #min(self.coord_x)
        spc=(max(self.coord_x)-min(self.coord_x))/10
        xdim = int(round((max(self.coord_x)-min(self.coord_x))/spc))
        ydim = int(round((max(self.coord_y)-min(self.coord_y))/spc))
        zdim = int(round((max(self.coord_z)-min(self.coord_z))/spc))

        grid.dimensions = (xdim, ydim, zdim)
        grid.dimensions = np.array(grid.dimensions) +1
        grid.origin = (min(self.coord_x), min(self.coord_y), min(self.coord_z)) # The bottom left corner of the data set
        grid.spacing = (spc, spc,spc) # These are the cell sizes along each axis
        pvfig.add_mesh(grid,show_edges=True, opacity=0.05)
        pvfig.add_mesh(poly, point_size=15.0, scalars=self.x.x, render_points_as_spheres=True,cmap='jet')
        #pvfig.add_bounding_box()
        pvfig.show_bounds(bounds=[min(self.coord_x), max(self.coord_x), 
                                        min(self.coord_y),max(self.coord_y),
                                        min(self.coord_z), 0],font_size=16)
        # pvfig.add_scalar_bar('Normalized Current density',width=0.1,vertical=False)
        #pvfig.view_isometric()
        #pvfig.camera_position = [2,1,0]
        pvfig.camera_position = [(3.926506482167545, -2.7565045773414787, 0.2967034053253811), (0.6, 0.5, -0.65), (-0.1420999573447329, 0.13972559831202394, 0.9799409978661837)]
        pvfig.add_points(self.pointsE)
        # pvfig.add_point_labels(self.pointsE,self.coordE[:69,0].astype(int), point_size=16, font_size=16)
        pvfig.add_point_labels(self.pointsE,self.coordE[:self.RemLineNb,0].astype(int), point_size=16, font_size=16)

        #pvfig.screenshot('SourcesColored_WR_' + str(self.wr) + '.png')
        # pvfig.screenshot('Pts_iCSD_knee'+ str(self.ObsName) + '.png')
        # pvfig.show(auto_close=True)  
        # pvfig.close()          
        #%matplotlib inline
        #plt.figure('test', figsize=(10,4))
        #plt.imshow(pvfig.image)
        #plt.show()




        interpolated = grid.interpolate(poly, radius=0.2)
        cmap = plt.cm.get_cmap('jet',10)
        contours = interpolated.contour()
        #slices = interpolated.slice_orthogonal(x=0.5, y=0.5, z=-np.array([0.1,0.2]))
        slices = interpolated.slice_along_axis(n=6, axis="z")

        # glyphs = interpolated.glyph(factor=5e-3, geom=pv.Sphere(radius=3.14))
        outline = interpolated.outline()
        # #print(self.args.Thresh)

        # # Make two points to construct the line between
        # a = [interpolated.bounds[1]/2, interpolated.bounds[3]/2, interpolated.bounds[4]/2]
        # b = [interpolated.bounds[1]/2, interpolated.bounds[3]/2, interpolated.bounds[5]]

        # Preview how this line intersects this mesh
        # line = pv.Line(a, b)
        # p = pv.Plotter()
        # p.add_mesh(interpolated, style="wireframe", color="w")
        # p.add_mesh(line, color="b")
        # p.screenshot('Test'+ str(self.ObsName) + '.png')
        # p.show()

        # dd=[]
        # LineInt =interpolated.plot_over_line(a, b, resolution=1000)
        # dd.append(LineInt[1])
        # dd.append(LineInt[2])
        # dd= np.array(dd)
        # f = open('testoverlineoutputT0PlantB.dat','w')
        # np.savetxt(f, dd.T, fmt='%1.5f %1.5f', delimiter='\t')   # X is an array
        # f.close()

        minTs= 0.3
        maxTs = 1
        cbarmax= 0.01
        # #minTs= int(self.args.Thresh[0])
        # #maxTs = int(self.args.Thresh[1])
        threshed = interpolated.threshold([minTs, 1])

        # threshed = interpolated.threshold_percent([0.15, 0.50], invert=False)
        # threshed = interpolated.threshold([0.01, 0.1])
        #bodies = threshed.split_bodies()
        #bodies.plot(show_grid=True, multi_colors=True, cpos=[-2, 5, 3],screenshot='Testbody')

        # print(str(minTs) + '<Threshold<' + str(1))
        # #print(threshed.cells)
        # #while not threshed.cells:
        
        while not threshed.cells.size:
           #print('change threshold')
            threshed = interpolated.threshold([minTs, 1])
            minTs= minTs - 1e-3
            #print('New threshold=' + str(minTs))
            #print(threshed.cells)
            continue

        print('New threshold=' + str(minTs))
        minTs=minTs-0.1  

        # p = pv.Plotter(window_size=[1000, 1000])
        # p.add_mesh(outline, color="k")
        # p.add_mesh(slices, show_scalar_bar=False)
        # #p.camera_position = [-2, 5, 3]
        # p.add_bounding_box()
        # p.show_bounds(grid=True,bounds=[min(self.coord_x), max(self.coord_x), 
        #                                 min(self.coord_y),max(self.coord_y),
        #                                 min(self.coord_z), 0],font_size=26,
        #                                 xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)')
        # p.update_scalar_bar_range([minTs, cbarmax])  
        # p.add_scalar_bar('Normalized Current density')
        # p.camera_position = [2,1,0]
        # #p.view_isometric()
        # p.add_points(self.pointsE)
        # p.add_point_labels(self.pointsE,self.coordE[:self.RemLineNb,0].astype(int), point_size=16, font_size=16)
        # p.screenshot('Current_iCSD_knee_glyphs'+ str(self.ObsName) + '.png')

        #########################################
        # Show the threshold
        #########################################
        p = pv.Plotter(shape=(1,2),window_size=[10000, 10000])
        # Show the theshold
        p.add_mesh(outline, color='k')
        opacity = [0.0, 0.1, 0.3, 0.45, 0.6, 0.9, 1]
        # opacity = [0.9, 0.9, 0.9, 0.9, 0.9, 0.9, 9]
        # #opacity = [0.0, 0.0, 0.1, 0.3, 0.6, 0.9, 1]

        p.add_mesh(threshed, show_scalar_bar=False,cmap='jet') #, opacity=opacity) # or add volume
        #p.add_mesh(threshed, show_scalar_bar=False, opacity = [0, 0, 0, 0.1, 0.3, 0.6, 1])
        #p.add_bounding_box()
        p.show_bounds(grid=True,bounds=[min(self.coord_x), max(self.coord_x), 
                                        min(self.coord_y),max(self.coord_y),
                                        min(self.coord_z), 0],font_size=26,
                                        xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)')
        #p.camera_position = [-2,5,3]
        p.add_scalar_bar('Normalized Current density')
        p.update_scalar_bar_range([minTs, cbarmax])  
        #p.view_xz()
        #p.camera_position = [2,1,0]
        p.view_isometric()
        p.add_points(self.pointsE)
        p.add_point_labels(self.pointsE,self.coordE[:self.RemLineNb,0].astype(int), point_size=16, font_size=16)

        # p.subplot(0,1)
        # p.add_mesh(threshed, show_scalar_bar=False,cmap='jet') #, opacity=opacity) # or add volume
        # #p.add_mesh(threshed, show_scalar_bar=False, opacity = [0, 0, 0, 0.1, 0.3, 0.6, 1])
        # #p.add_bounding_box()
        # p.show_bounds(grid=True,bounds=[min(self.coord_x), max(self.coord_x), 
        #                                 min(self.coord_y),max(self.coord_y),
        #                                 min(self.coord_z), 0],font_size=26,
        #                                 xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)')
        # #p.camera_position = [-2,5,3]
        # p.add_scalar_bar('Normalized Current density')
        # p.update_scalar_bar_range([minTs, cbarmax])  
        # p.view_xz()
        # #p.camera_position = [2,1,0]
        # p.add_points(self.pointsE)
        # p.add_point_labels(self.pointsE,self.coordE[:69,0].astype(int), point_size=16, font_size=16)

        # p.screenshot('Current_iCSD_knee_thresh'+ str(self.ObsName) + '.png')

        #########################################
        # Show the contour
        #########################################



        p = pv.Plotter(notebook=False,window_size=[1000, 1000])
        # # pv.BackgroundPlotter()
        # # p = pv.BackgroundPlotter(notebook=False,window_size=[10000, 10000])

        p.add_mesh(outline, color='k')
        p.add_mesh(contours, show_scalar_bar=False, opacity= opacity,cmap='jet')
        # p.camera_position = [(3.926506482167545, -2.7565045773414787, 0.2967034053253811), (0.6, 0.5, -0.65), (-0.1420999573447329, 0.13972559831202394, 0.9799409978661837)]
        p.add_bounding_box()
        p.show_bounds(grid=True,bounds=[min(self.coord_x), max(self.coord_x), 
                                        min(self.coord_y),max(self.coord_y),
                                        min(self.coord_z), 0],
                                        font_size=60,
                                        xlabel='X (m)', ylabel='Y (m)', zlabel='Z (m)')
        #p.view_xz()
        #p.view_isometric()
        p.add_axes()
        p.show_axes()
        p.add_scalar_bar('Normalized Current density',width=0.25,vertical=False,position_x=0.3)
        p.update_scalar_bar_range([minTs, cbarmax])  

        p.add_points(self.pointsE)
        # p.add_point_labels(self.pointsE,self.coordE[:self.RemLineNb,0].astype(int), point_size=15, font_size=35,
        #     shape_opacity=0.01, margin=4.)
        p.add_point_labels(self.pointsE,self.coordE[:self.RemLineNb,0].astype(int), point_size=150, font_size=35,
            shape_opacity=0.01, margin=4.)
        #poly['My Labels'] = ['E{}'.format(i+1) for i in range(poly.n_points)]

        # #p.add_point_labels(poly,   'My Labels', point_size=10, font_size=6)
        # # hsize=256000
        # # p.ren_win.OffScreenRenderingOn()
        # #p.enable_anti_aliasing()
        # # p.screenshot("try.png", transparent_background=True, window_size=[hsize, int(hsize*p.window_size[1]/p.window_size[0])])
        # # p.screenshot('Current_iCSD_knee_cont'+ str(self.ObsName) + '.png', window_size=[hsize, int(hsize*p.window_size[1]/p.window_size[0])])
        # # p.ren_win.SetSize([1000, 800])
        # # p.ren_win.OffScreenRenderingOff()
        # # p.ren_win.Render()
        # # p.show(auto_close=True)        
        p.screenshot('Current_iCSD_knee_cont'+ str(self.ObsName) + '.png')
        # p.save('Current_iCSD_knee_cont.vtkjs')

        # dd=[]
        # contInt =contours.plot_over_line(a, b, resolution=1000)
        # dd.append(contInt[1])
        # dd.append(contInt[2])
        # dd= np.array(dd)
        # f = open('testoverlineoutputT0PlantB.dat','w')
        # np.savetxt(f, dd.T, fmt='%1.5f %1.5f', delimiter='\t')   # X is an array
        # f.close()


        # pv_test = pv.BackgroundPlotter()
        # pv_test.add_mesh(contours, show_scalar_bar=False, opacity= opacity,cmap='jet')
        #p.show()

        # contours.save('contours.vtk')
        # interpolated.save('interpolated.vtk')
        # poly.save('poly.vtk')
        # grid.save('grid.vtk')

   
        # Grab the largest connected volume present
        p = pv.Plotter(shape=(1,2),window_size=[1000, 1000])
        threshed = interpolated.threshold([minTs, 1])
        largest = threshed.connectivity(largest=True)
        # # or: largest = threshed.extract_largest()
        # # Get volume as numeric value
        # large_volume = largest.volume
        # p.add_mesh(large_volume)
        # p.show(show_grid=True, clim= [minTs, 1], cpos=[-2, 5, 3],
        #    screenshot='Current_iCSD_knee_Ctvity'+ str(self.ObsName) + '.png')
        largest.plot(show_grid=True, clim= [minTs, cbarmax], cpos=[-2, 5, 3], 
            screenshot='Current_iCSD_knee_Ctvity'+ str(self.ObsName) + '.png',window_size=[1000, 1000])

        # p = pv.Plotter()
        # p.add_mesh(largest, color='#965434')
        # p.add_mesh(slices, show_scalar_bar=False)
        # p.add_mesh(outline)
        # #p.camera_position = cpos
        # p.show()
        # p.screenshot('Test__'+ str(self.ObsName) + '.png')

        # Linelargest =largest.plot_over_line(a, b, resolution=1000)

        #vol = interpolated.threshold_percent(70, invert=0)
        #surf = vol.extract_geometry()
        #smooth = surf.smooth(n_iter=1000)
        #cpos = [-2, 5, 3]
        #smooth.plot(show_grid=True,show_edges=True, cpos=cpos,
        #    screenshot='Current_iCSD_knee_vol'+ str(self.ObsName) + '.png',window_size=[10000, 10000])
        #vol.plot(show_edges=True, cpos=cpos)

    def _plotPareto_(self):
        self.p, self.ax = plt.subplots()
        #self.p = plt.figure('pareto', figsize=(10,4))
        plt.plot(self.pareto_list_FitRes, self.pareto_list_RegRes, 'or')
        self.ax.annotate('Wr=' + str(int(self.wr)), xy=(float(np.asarray(self.pareto_list_FitRes)[self.IdPtkneew]), 
                                         float(np.asarray(self.pareto_list_RegRes)[self.IdPtkneew])), 
                                      xytext=(float(np.asarray(self.pareto_list_FitRes)[self.IdPtkneew])+max(self.pareto_list_FitRes)/3, 
                                         float(np.asarray(self.pareto_list_RegRes)[self.IdPtkneew])+max(self.pareto_list_RegRes)/3),
                                      arrowprops=dict(facecolor='black', shrink=0.05))
        plt.plot(float(np.asarray(self.pareto_list_FitRes)[self.IdPtkneew]), float(np.asarray(self.pareto_list_RegRes)[self.IdPtkneew]), 'og')

        ax = plt.gca()
        ax.tick_params(axis = 'both', which = 'both', direction = 'out')
        ax.grid()
        plt.xlabel('Residual')
        plt.ylabel('Roughness')
        plt.tight_layout()
        plt.savefig('ParetoFront.png', dpi = 600)

    def _plotFIT_(self):
        plt.figure()
        plt.subplot(121)
        plt.plot(self.x.fun[:204] + self.b_w[:204], 'or')
        plt.plot(self.b_w[:204], 'ob')
        plt.subplot(122)
        plt.plot(self.x.fun[:204] + self.b_w[:204], self.b_w[:204], 'or')
        plt.tight_layout()
        plt.show()

    ### LOAD 
    def load_names(self):
        self.ObsName, file_extension = os.path.splitext(str(self.args.obs))

    def _load_coord3d_(self):
        self.coord = np.loadtxt(self.args.coord_file)
        self.coord_x, self.coord_y, self.coord_z = self.coord[:, 0], self.coord[:, 1], self.coord[:, 2]

    def _load_geom_(self):
        geom_files = [f for f in os.listdir('.') if f.endswith('.geom')]
        if len(geom_files) != 1:
            raise ValueError('should be only one geom file in the current directory')

        fileNameElec = geom_files[0]
        #fileNameElec= '*.geom'
        # RemLineNb= search_string_in_file(fileNameElec, '#Remote')

        line_number = 0
        line_of_results = []
        # Open the file in read only mode
        with open(fileNameElec, 'r') as read_obj:
            # Read all lines in the file one by one
            for line in read_obj:
                # For each line, check if line contains the string
                line_number += 1
                if ('#Remote') in line:
                    # If yes, then add the line number & line as a tuple in the list
                    line_of_results.append((line_number))
        self.RemLineNb= int(line_of_results[0])-1

        self.coordE = np.loadtxt(fileNameElec)
        # print(self.coordE)
        self.coord_xE, self.coord_yE, self.coord_zE = self.coordE[:, 0], self.coordE[:, 1], self.coordE[:, 2]
        # self.pointsE= np.vstack(self.coordE[:69,1:4])
        self.pointsE= np.vstack(self.coordE[:self.RemLineNb,1:4])

    def _load_F1w_(self):
        self.F1w = np.loadtxt(self.args.F1w_file)

    def _load_coord_(self):
        self.coord = np.loadtxt(self.args.coord_file)
        self.coord_x, self.coord_y = self.coord[:, 0], self.coord[:, 1]

    def _load_obs_(self):
        self.b = np.loadtxt(self.args.obs)

    def _load_sim_(self):
        self.A = np.loadtxt(self.args.sim)

    ### MAKE LINEAR SYSTEM 

    def _check_nVRTe_(self):
        if self.A.shape[0] / self.b.shape[0] ==  self.coord.shape[0]:
            self.nVRTe = self.coord.shape[0]
        else:
            raise ValueError('### dimensions of the files do not agree')

    def _reshape_A_(self):
        self.A = self.A.reshape((-1, self.nVRTe), order = 'F')

    def _obs_w_(self):
        """weight the observations, can also ignore observations by setting w = 0"""
        if self.args.obs_err == 'const':
            print('const weight')
            self.obs_w = np.ones(self.b.shape[0])
        elif self.args.obs_err == 'sqrt':
            print('sqrt obs weight')
            self.obs_w = 1 / np.sqrt(np.abs(self.b))
            self.obs_w[self.obs_w >= self.args.errRmin] = 1


    ### Individual Misfit  
    def normF1(self):
        self.F1=[]
        for i in range(np.shape(self.A)[1]):
            F1i = LA.norm((self.b-self.A[:,i]))
            self.F1.append(F1i)
        self.norm_F1 = (self.F1 - min(self.F1)) / (max(self.F1) - min(self.F1))

    def misfit_2_initialX0(self):
        alpha=1
        self.x0F1=1./((self.norm_F1+1)*(self.norm_F1+1))
        self.x0F1_sum= alpha*self.x0F1/sum(self.x0F1)
        #print('max x0F1_sum =' + str(max(self.x0F1_sum)))
        #print('sum x0F1_sum =' + str(sum(self.x0F1_sum)))

    ### CONSTRAIN 

    def _con_A_(self):
        self.con_A = np.ones(self.A.shape[1])

    def _con_b_(self):
        self.con_b = np.ones(1)

    def _con_w_(self):
        self.con_w = np.ones(1) * self.args.wc

    ### REGULARIZATION 

    def _nx_ny_(self):
        """find number of nodes in each direction, has to be a regular grid"""
        self.nx = np.unique(np.round(self.coord[:, 0], 3)).shape[0]
        self.ny = np.unique(np.round(self.coord[:, 1], 3)).shape[0]
        print(self.nx, self.ny)

    def _nx_ny_nz(self):
        """find number of nodes in each direction, has to be a regular grid"""
        self.nx = np.unique(np.round(self.coord[:, 0], 3)).shape[0]
        self.ny = np.unique(np.round(self.coord[:, 1], 3)).shape[0]
        self.nz = np.unique(np.round(self.coord[:, 2], 3)).shape[0]
        # print(self.nx, self.ny)
        
    def regularize_A_UnstructuredMesh3d(self): # should also work for the 2d case
        print('regularize_A_UnstructuredMesh3d')
        reg = []
        for VRTEnb in range(self.nVRTe):
            dist =  np.linalg.norm(self.coord[VRTEnb]-self.coord, axis=1)
            closest = np.argsort(dist)
            k = self.args.k  # For each point, find the k closest current sources
            #print('k=' + str(k))
            Ind = closest[1:k+1]
            row = np.zeros(self.nVRTe) # add a line to the regularisation A with k non-null coefficients
            #row[Ind]= -1/k # ponderate coeff for k closest current sources
            knorm = dist[closest[1:k+1]]/dist[closest[1:k+1]].sum(axis=0,keepdims=1)
            row[Ind]= knorm
            row[VRTEnb]= 1 # = one for the actual current source
            reg.append(row)
            test=[1]
            mask = np.in1d(test, VRTEnb)
            #print(mask)
            if mask.any()==True: 
                self.fc = plt.figure('TEST regularisation')
                ax = self.fc.add_subplot(111, projection='3d')
                ax.scatter(self.coord[VRTEnb,0], self.coord[VRTEnb,1], self.coord[VRTEnb,2], linewidths=12,
                           facecolor = 'green', edgecolor = 'green')
                ax.scatter(self.coord[Ind,0], self.coord[Ind,1], self.coord[Ind,2], linewidths=12,
                           facecolor = 'red', edgecolor = 'red')
                ax.set_xlim([min(self.coord_x),max(self.coord_x)])
                ax.set_ylim([min(self.coord_y),max(self.coord_y)])
                ax.set_zlim([min(self.coord_z),max(self.coord_z)])
                self.fc.savefig('TEST regularisation', dpi = 600)
                #plt.show()
                #plt.close()
        self.reg_A = np.array(reg)

    def _regularize_A_3d(self):
        """create and append rows for spacial regularization to A"""
        reg = []
        vrte = range(1, self.nVRTe + 1)
        vrte = np.reshape(vrte,(self.ny, self.nx, self.nz))
        for z in range(self.nz):
	        for y in range(self.ny):
	            for x in range(self.nx):
	                minus = vrte[y, x]
	                if x + 1 in range(self.nx):
	                    plus = vrte[y , x + 1]
	                    row = np.zeros(self.nVRTe, int)
	                    row[minus -1] = - 1
	                    row[plus -1] = + 1
	                    reg.append(row)
	                    if y + 1 in range(self.ny):
	                        plus = vrte[y + 1, x]
	                        row = np.zeros(self.nVRTe)
	                        row[minus -1] = - 1
	                        row[plus -1] = + 1
	                        reg.append(row)
        self.reg_A = np.array(reg)

    def _regularize_A_(self):
        """create and append rows for spacial regularization to A"""
        reg = []
        vrte = range(1, self.nVRTe + 1)
        vrte = np.reshape(vrte,(self.ny, self.nx))
        for y in range(self.ny):
            for x in range(self.nx):
                minus = vrte[y, x]
                if x + 1 in range(self.nx):
                    plus = vrte[y , x + 1]
                    row = np.zeros(self.nVRTe, int)
                    row[minus -1] = - 1
                    row[plus -1] = + 1
                    reg.append(row)
                    if y + 1 in range(self.ny):
                        plus = vrte[y + 1, x]
                        row = np.zeros(self.nVRTe)
                        row[minus -1] = - 1
                        row[plus -1] = + 1
                        reg.append(row)
        self.reg_A = np.array(reg)

    def _regularize_b_(self):
        """append 0's to b"""
        #self.reg_b = np.zeros(self.reg_A.shape[0])
        if self.args.F1w_file is not None:
            self.reg_b = np.ones(self.reg_A.shape[0])
        else:
            self.reg_b = np.zeros(self.reg_A.shape[0])

    def _regularize_w_(self):
        """create vector with weights, the length is determined by the number of regul rows in A"""
        #self.reg_w = np.ones(self.reg_A.shape[0]) * self.args.wr
        if icsd.args.pareto == True:
            self.reg_w = np.ones(self.reg_A.shape[0]) * self.wr
            if self.args.F1w_file is not None:
                self._load_F1w_()
                print('!! Wm*m0 !!')
                print('self.F1w= ' + str(np.shape(self.F1w)))
                print('np.ones= ' + str(np.shape(np.ones(self.reg_A.shape[0]))))
                self.reg_w_x0 = np.ones(self.reg_A.shape[0]) * self.wr * self.F1w
                #np.savetxt('reg_w_x0.txt', self.reg_w_x0)
            #np.savetxt('reg_w.txt', self.reg_w)
            #print('reg_w= ' + str(np.shape(self.reg_w)))
            #print('reg_w_x0= ' + str(np.shape(self.reg_w_x0)))

        else:
            self.reg_w = np.ones(self.reg_A.shape[0]) * self.args.wr
        #self.reg_w = np.ones(self.reg_A.shape[0]) * self.args.wr

    ### VERTICAL STACK EQUATIONS

    def _stack_A_(self):
        self.A_s = np.vstack((self.A, self.con_A, self.reg_A)) # TO DO add Relative smallness regularisation ! 

    def _stack_b_(self):
        self.b_s = np.concatenate((self.b, self.con_b, self.reg_b))

    def _stack_w_(self):
        """create vector with weights for observation, constrain, and regularization
        then use it as diagonal for the weight matrix"""
        w = np.concatenate((self.obs_w, self.con_w, self.reg_w))
        #np.savetxt('wb.txt', w)
        W = np.zeros((w.shape[0], w.shape[0]))
        np.fill_diagonal(W, w)
        self.W_s = W

        if self.args.F1w_file is not None:
            print('write wb')
            wb = np.concatenate((self.obs_w, self.con_w, self.reg_w_x0))
            #np.savetxt('wb_x0.txt', wb)
            Wb = np.zeros((wb.shape[0], wb.shape[0]))
            np.fill_diagonal(Wb, wb)
            self.W_s_x0 = Wb
            #np.savetxt('self.W_s_x0.txt', self.W_s_x0)

        #np.savetxt('self.W_s.txt', self.W_s)

    ### APPLY WEIGHTS 

    def _weight_A_(self):
        """Apply the weights to A"""
        self.A_w = np.matmul(self.W_s, self.A_s)

    def _weight_b_(self):
        """Apply the weights to b"""

        if self.args.F1w_file is not None:
            #self.b_w = np.matmul(self.b_s, self.W_s_x0)
            self.b_w = np.dot(self.b_s, self.W_s_x0)
            #np.savetxt('b_w_x0.txt',  self.b_w)
            #np.savetxt('b_s_x0.txt',  self.b_s)
            #np.savetxt('W_s_x0.txt',  self.W_s_x0)
            #print('W_s_x0= ' + str(np.shape(self.W_s_x0)))
            #print('W_s= ' + str(np.shape(self.W_s)))

        else:
            self.b_w = np.matmul(self.b_s,  self.W_s)
            #np.savetxt('b_w.txt', self.W_s)
            #np.savetxt('b_w.txt',  self.b_w)
            #np.savetxt('b_s.txt',  self.b_s)
            #np.savetxt('W_s.txt',  self.W_s)

    ### RUN INITILIZATION METHODS

    def icsd_init(self):
        """ these functions are called only once, even for pareto,
        as they do not depend on the regularization weight"""
        # load
        icsd.load_names()
        self._load_coord3d_() # CHANGE HERE
        self._load_geom_()
        self._load_obs_()
        self._load_sim_()
        # check vector sizes
        self._check_nVRTe_()
        # reshape VRTe vector into matrix A 
        self._reshape_A_()
        # define mode to weights the data        
        self._obs_w_()
        # constrain (curent conservation)
        self._con_A_()
        self._con_b_()
        self._con_w_()
        # get VRTe grid size and make grid
        # self._nx_ny_()
        self._nx_ny_nz()
        #self._mkGrid_XI_YI_() # COMMENTED
        # append spatial regularization
        # self.regularize_A_UnstructuredMesh3d() # CHANGE
        self._regularize_A_3d() # CHANGE
        self._regularize_b_()
        # stack data, constrain, and regularization 
        self._stack_A_()
        self._stack_b_()


    def _prepare4iCSD_(self):
        """ this fucntion are called for each weight, keep them separated for pareto"""
        # create regularization part of the weight matrix
        self._regularize_w_()
        self._stack_w_()
        # apply weights with matrix multiplication
        self._weight_A_()
        self._weight_b_()

    ### LSQ 

    def _iCSD_(self):
        """solve linear system, given A matrix (VRTe, constrain, regul) and b (observations)"""
        #self.x = lsq_linear(self.A_w, self.b_w, bounds = (0, 1))
        #print(type(self.x.x))
        # with open('X' + str(self.args.wr)+ '.txt', 'w') as f:
        #     for item in self.x.x:
        #         f.write("%f\n" % item)
        """solve linear system, given A matrix (VRTe, constrain, regul) and b (observations)"""
        print(self.args.F1w_file)
        if self.args.F1w_file is not None:
            print('Initial guess')
            a = self.A_w 
            b = self.b_w
            def func(x, a, b):
                return (b - np.dot(a, x))
            self.x = least_squares(func, x0=self.F1w, bounds = (0, 1), args=(a, b)) # Add initial guess
            print('Sum=' + str(np.sum(self.x.x)))

        else:
            print('No initial guess')
            self.x = lsq_linear(self.A_w, self.b_w, bounds = (0, 1))
            print('*' * 20)
        icsd._Export_sol_()


    def _ResidualAnalysis_(self):
        fitting_res = self.x.fun[0 : self.b.shape[0]]
        # constrain_res = self.x.fun[self.b.shape[0] + 1] / self.args.wc
        if icsd.args.pareto == True:
            regularization_res = self.x.fun[self.b.shape[0] + 1 :] / self.wr # this constrain is not included in the reg function
        else:
            regularization_res = self.x.fun[self.b.shape[0] + 1 :] / self.args.wr # this constrain is not included in the reg function
        self.reg_sum = np.sum(np.square(regularization_res)) # ||Ax - b||2
        self.fit_sum = np.sum(np.square(fitting_res)) 
        self.pareto_list_FitRes.append(self.fit_sum)
        self.pareto_list_RegRes.append(self.reg_sum)

    def DetectKneePt(self):
        self.kn = KneeLocator(self.pareto_list_FitRes,self.pareto_list_RegRes, 
                         curve='convex', direction='decreasing')
#        print('knee xloc=' + str(self.kn.knee))
        self.IdPtkneew= np.where(self.kn.knee==self.pareto_list_FitRes)[0]
#        print('Id pt knee=' + str(self.IdPtkneew))
        self.pareto_weights[self.IdPtkneew]
#        print('Wr knee=' + str(self.pareto_weights[self.IdPtkneew]))
        if len(self.IdPtkneew)<1:
             self.IdPtkneew=1
             print('No knee detection possible, put 1 as default')
#            raise ValueError('No knee detection possible')
    ### MAINS

    def run_pareto(self):
        """
        run iCSD multiple times while changing the weights to explore the L-curve
        """
        self.pareto_weights = np.linspace(self.args.pareto_MinErr, self.args.pareto_MaxErr, self.args.pareto_nSteps)
        print('pareto weights are\n', self.pareto_weights)

        self.pareto_list_FitRes = []
        self.pareto_list_RegRes = []
        self.ObsName, file_extension = os.path.splitext(str(self.args.obs))
        with PdfPages('iCSD_pareto.pdf') as pdf:
            for w in self.pareto_weights:
                self.wr = w
                print('-' * 80)
                print('Pareto mode, running new iCSDwith regularizaiton weight: ', self.wr)
                self._prepare4iCSD_()
                self._iCSD_()
                print(min(self.x.x), max(self.x.x))
                self._plotCSD_()
                self._ResidualAnalysis_()
                pdf.savefig(self.f)
                plt.close(self.f)

            icsd.DetectKneePt()
            print(self.IdPtkneew)
            self.wr=float(self.pareto_weights[self.IdPtkneew])
            icsd._plotPareto_()
            self.p.savefig('Pareto_Lcurve'+ str(self.ObsName), dpi = 450)
            pdf.savefig(self.p)
            plt.close(self.p)
            #self._Interp3d_()
            self._Interp3dPyvista_()

            icsd._prepare4iCSD_()
            icsd._iCSD_()
            icsd._plotCSD_()
            self.f.savefig('iCSD_knee'+ str(self.ObsName), dpi = 450)
            plt.show()

    def run_single(self): 
        self.ObsName, file_extension = os.path.splitext(str(self.args.obs))
        self._prepare4iCSD_()
        self._iCSD_()
        self._plotCSD_()
        self.f.savefig('iCSD', dpi = 450)
        plt.show()
        # print('je suis la')
        # self._Interp3d_()
        self._Interp3dPyvista_()




if __name__ == "__main__":

    icsd = iCSD_class()
    icsd.get_args()
    icsd.icsd_init()
    # run iCSD, either pareto or single or TL (time lapse)
    if icsd.args.pareto == True:
        print('run pareto')
        icsd.run_pareto()

    elif icsd.args.pareto == False:
        print('run single inversion')
        icsd.run_single()

    elif icsd.args.TL == True:
        print('run single inversion')
        icsd.run_TL
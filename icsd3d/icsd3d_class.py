import os


import numpy as np
from scipy.optimize import lsq_linear, least_squares

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm

from kneed import KneeLocator

import matplotlib.pyplot as plt
import pyvista as pv

from plotters.mpl_plot import plotCSD3d, plotCSD2d, scatter3d, scatter2d, plotContour2d
from plotters.pv_plot import plotCSD3d_pyvista
from inversion.smoothing import *
from inversion.priorM0 import *
from importers.read import *
from exporters.save import * 
from gridder.mkgrid import mkGrid_XI_YI 

class iCSD3d_Class():
    
    
    """ Create a icsd inversion object.
    
    Parameters
    ----------
    coord_file : str, mandatory
        coordinates of the VRTe for plotting
    wr : float, optional
        Weight regularization
    wc : float, optional
    """
    def __init__(self,dirName):
        self.dirName = dirName
        self.clim = []
       # load
        self.type='2d'
        self.sim='VRTeSim.txt'
        self.obs='ObsData.txt'
        self.coord_file='VRTeCoord.txt'
        self.wr=25 #weight regularization
        self.wc=10000 #current conservation constrain, sum current fractions = 1
        self.sc=[] #coordinates of the sources, format = x1,y1 x2,y2'
        self.retElec=None #coordinates of the return electrode, format = x1,y1')
        self.pareto=False #if True run many icsd to explore the Pareto front
        self.errRmin=1 #min R to which the err is applied before passing to constant error
        self.pareto_MinErr=0.001
        self.pareto_MaxErr=1
        self.pareto_nSteps=10
        self.obs_err='const' #const or sqrt - choose between constant weight and w = 1/sqrt(abs(obs))
        # IMPLEMENT obs_err based on reciprocal analysis i.e. estimated standard deviation of the data errors;                         % estimated standard deviation of the traveltime data errors
        self.k=4  # For each point, find the k closest current sources
        self.TL=False # Time lapse inversion (see slides VIU Venice to implement)
        self.x0_prior=False #  relative smallness regularization as a prior criterion for the inversion; i.ethe algorithm minimizes ||m−m0||2
        self.x0_ini_guess=False # initial guess
        self.knee=False # L-curve knee automatic detection
        self.KneeWr=[]
        self.regMesh='strc' # strc or unstrc
        self.plotElecs=False 
        self.alphax0=1 # weight on model smallness relative to m0
        self.inix0=None # or 'cst' if constant vector *0.1
        self.alphaSxy=False # weight on model smoothness in z-direction 
        self.alphaSx=1 # weight on model smoothness in x-direction
        self.alphaSy=1 # weight on model smoothness in y-direction
        self.alphaSz=1 # weight on model smoothness in z-direction [TO IMPLEMENT]
        self.mesh=None # mesh3d .vtk to plot with the results of icsd
        self.gif3d=False # create gif orbit
        self.title=None #


    def icsd_init(self):
        """ these functions are called only once, even for pareto,
        as they do not depend on the regularization weight"""
        # load
        self.createdirs()
        # self.load_coord()
        # self.load_coord(path, filename, dim)
        if self.type == '2d':
            self.coord_x, self.coord_y, self.coord = load_coord(self.path2load, self.coord_file, self.type)
        else:
            self.coord_x, self.coord_y, self.coord_z, self.coord = load_coord(self.path2load, self.coord_file, self.type)
        self.b = load_obs(self.path2load, self.obs)
        self.A = load_sim(self.path2load, self.sim)
        if self.plotElecs==True:
            self.RemLineNb, self.Injection, self.coordE, self.pointsE= load_geom(self.path2load) # geometry file containing electrodes position includinf remotes 
        # check vector sizes
        self.check_nVRTe()
        # # reshape VRTe vector into matrix A 
        self.reshape_A()
        # # define mode to weights the data        
        self.obs_w_f()
        # # constrain (curent conservation)
        self._con_A_f()
        self._con_b_f()
        self._con_w_f()
        self.XI, self.YI= mkGrid_XI_YI(self.coord_x,self.coord_y)

        # # append spatial regularization
        self._parseModelReg()
        # self.parseDataReg()
        self._regularize_b()
        # # stack data, constrain, and regularization 
        self._stack_A()
        self._stack_b()
        
        
      
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
             


  





    ### mkdirs 
    def createdirs(self):
        self.path2load = self.dirName
        print(self.path2load)
        self.path2save= self.path2load + 'figs'
        # self.path2save= self.dirName + 'fig/'
        # print(self.path2save)
        # cwd = os.getcwd()
        # print(cwd+self.path2save)
        # try:
        #     # Create target Directory
        #     os.mkdir(self.path2save)
        #     print("Directory " , self.path2save ,  " Created ") 
        # except FileExistsError:
        #     print("Directory " , self.path2save ,  " already exists")
            
    ### MAKE LINEAR SYSTEM 

    def check_nVRTe(self):
        if self.A.shape[0] / self.b.shape[0] ==  self.coord.shape[0]:
            self.nVRTe = self.coord.shape[0]
        else:
            raise ValueError('### dimensions of the files do not agree')

    def reshape_A(self):
        self.A = self.A.reshape((-1, self.nVRTe), order = 'F')

    def obs_w_f(self):
        """weight the observations, can also ignore observations by setting w = 0"""
        if self.obs_err == 'const':
            self.obs_w = np.ones(self.b.shape[0])
        elif self.obs_err == 'sqrt':
            self.obs_w = 1 / np.sqrt(np.abs(self.b))
            print('Selfb = 0 could be a problem, check presence of 0 and filter if needed')
            self.obs_w[self.obs_w >= self.errRmin] = 1    
        elif self.obs_err == 'reciprocals': #[TO IMPLEMENT]
            self.obs_w = 1 / np.sqrt(np.abs(self.sd_rec))

    ### CONSTRAIN 

    def _con_A_f(self):
        """Set current conservation constrainst on A (rows of ones)"""
        self.con_A = np.ones(self.A.shape[1])

    def _con_b_f(self):
        """Set current conservation constrainst on b"""
        self.con_b = np.ones(1)

    def _con_w_f(self):
        """Set current conservation constrainst weight; default is wc=1e6 """
        self.con_w = np.ones(1) * self.wc
        
    ### REGULARIZATION 



    def regularize_A(self):
        """create and append rows for spatial regularization to A"""
        print('Reg A (Luca"s implementation)')
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
        


    def regularize_A_x_y_z(self):
        self.estimateM0()
        self.reg_Axz = []
        self.reg_Ayz = []
        for z in range(self.nz):
            ncells = self.nx*self.ny;       
            Dx=diags([1, -2, 1], [-1, 0, 1], shape=(ncells, ncells)).todense()
            idx=np.arange(0,len(Dx),self.nx)
            Dx[idx,:] = 0
            idx2=np.arange(self.nx,len(Dx),self.nx)
            Dx[idx2,:] = 0
            
            Dy=diags([1, -2, 1], [-self.nx, 0, self.nx], shape=(ncells, ncells)).todense()
            idy=np.arange(0,self.nx)
            Dy[idy,:] = 0
            idy2=np.arange(((self.ny-1)*self.nx+1),(self.nx)*(self.ny))
            Dy[idy2,:] = 0
            
            self.reg_Ax = self.alphaSx*np.array(Dx.transpose()*Dx)
            self.reg_Ay = self.alphaSy*np.array(Dy.transpose()*Dy)
            print('reg_Ax_shape ='+ str(np.shape(self.reg_Ax)))

            # self.reg_Axz=np.append([[self.reg_Axz],[self.reg_Ax]], axis=0)
            self.reg_Axz.append(self.reg_Ax)
            self.reg_Ayz.append(self.reg_Ay)
            print('reg_Axz ='+ str(np.shape(self.reg_Axz)))

        self.reg_Ax= np.reshape(self.reg_Axz,[160,20])
        self.reg_Ay= np.reshape(self.reg_Ayz,[160,20])
        print('reshape ='+ str(np.shape(self.reg_Ax)))
        print('length reg_Ax =' + str(len(self.reg_Ax)))
        print('length reg_Ay =' + str(len(self.reg_Ay)))
        print('length reg_Axz =' + str(len(self.reg_Axz)))
        print('length reg_Ayz =' + str(len(self.reg_Ayz)))
        
    ### RELATIVE SMALLNESS conditions (m-m0)
        
    def regularize_smallnessX0(self):
        """ Create relative smallness instance 
        
        .. math :: X_{0} = A*\alpha_{x_{0}}

        Parameters
        ------------
        self
        """
        if self.alphaSxy==True:
            self.reg_smallx0 = np.ones(self.reg_Ax.shape)*self.alphax0
        else:
            self.reg_smallx0 = np.ones(self.reg_A.shape)*self.alphax0
        # self.reg_smallx0 = np.ones(self.x0.shape)*self.alphax0
        
    def regularize_sum_AX0(self):
        """sum smallness and spatial regularisation

        .. math:: W_{m}=\alpha_{s}I+{D_{x}}^{T}D_{x} + D_{z}}^{T}D_{z}

        Parameters
        ------------
        self
        """
        if (self.alphaSxy==True and self.x0_prior==True):
            print("""sum small x0, reg Ax, reg Ay""") 
            self.reg_A= self.reg_smallx0 + self.reg_Ax + self.reg_Ay
        elif (self.alphaSxy==True and self.x0_prior==False):
            print("""sum reg Ax, reg Ay""") 
            self.reg_A= self.reg_Ax + self.reg_Ay
        elif (self.alphaSxy==False and self.x0_prior==True):
            print("""reg_A= reg_A + small x0""") 
            self.reg_A= self.reg_A + self.reg_smallx0
            
    def _regularize_b(self):
        self.reg_b = np.zeros(self.reg_A.shape[0])

    def regularize_w(self):
        """create vector with weights, the length is determined by the number of regul rows in A"""
        self.reg_w = np.ones(self.reg_A.shape[0]) * self.wr
        if self.x0_prior==True:
            print('reg Wm (smallness + spatial reg) * lambda=' + str(self.wr))
            self.reg_w_0_b = np.ones(self.reg_A.shape[0]) * self.x0 * self.wr
            self.reg_w_0_A = np.ones(self.reg_A.shape[0])* self.wr
            
    ### VERTICAL STACK EQUATIONS
    def _stack_A(self):
        """Stack A (green fcts), constrainsts and regularisation"""
        print('shape A=' + str(np.shape(self.A)))
        print('shape con_A=' + str(np.shape(self.con_A)))
        print('shape reg_A=' + str(np.shape(self.reg_A)))
        self.A_s = np.vstack((self.A, self.con_A, self.reg_A))

    def _stack_b(self):
        """Stack b, constrainsts and regularisation"""
        self.b_s = np.concatenate((self.b, self.con_b, self.reg_b))

    def _stack_w(self):
        """create vector with weights for observation, constrain, and regularization
        then use it as diagonal for the weight matrix"""     
        w = np.concatenate((self.obs_w, self.con_w, self.reg_w))
        W = np.zeros((w.shape[0], w.shape[0]))
        np.fill_diagonal(W, w)
        self.W_s = W
        
        if self.x0_prior==True: # if relative smallness 
            wa = np.concatenate((self.obs_w, self.con_w, self.reg_w_0_A))
            wb = np.concatenate((self.obs_w, self.con_w, self.reg_w_0_b))
            W = np.zeros((wa.shape[0], wa.shape[0]))
            np.fill_diagonal(W, wa)
            self.W_s_A = W
            np.fill_diagonal(W, wb)
            self.W_s_b = W

    ### APPLY WEIGHTS 

    def _weight_A(self):
        """Apply the weights to A"""
        if self.x0_prior==True:
            self.A_w = np.matmul(self.W_s_A, self.A_s)
        else:
            self.A_w = np.matmul(self.W_s, self.A_s)
            
    def _weight_b(self):
        """Apply the weights to b"""
        if self.x0_prior==True:
            self.b_w = np.matmul(self.b_s,  self.W_s_b)
        else:
            self.b_w = np.matmul(self.b_s,  self.W_s)
    
    ### PREPARE FOR ICSD

    def prepare4iCSD(self):
        """ this function is called for each weight, keep them separated for pareto"""
        # create regularization part of the weight matrix
        if self.x0_ini_guess==True:
            self.estimateM0()
        elif self.x0_prior==True:
            self.estimateM0()
        self.regularize_w()
        self._stack_w()
        # apply weights with matrix multiplication
        self._weight_A()
        self._weight_b()

    ### LSQ 

    def iCSD(self):
        """solve linear system, given A matrix (VRTe, constrain, regul) and b (observations)"""
        if self.x0_ini_guess==False:
            print('No initial guess')
            self.x = lsq_linear(self.A_w, self.b_w, bounds = (0, 1))
            print('*' * 20)
            print('CURRENT Sum=' + str(np.sum(self.x.x)))
            # TO IMPLEMENT RETURN JAC Matrice to evaluate MALM sensitivity
        else:
            print('Initial guess x0')
            a = self.A_w 
            b = self.b_w
            def func(x, a, b):
                return (b - np.dot(a, x))
            self.x = least_squares(func, x0=self.x0, bounds = (0, 1), args=(a, b)) # Add initial guess
            print('CURRENT Sum=' + str(np.sum(self.x.x)))
        Export_sol(self.coord, self.x.x, self.type,path=self.path2load,filename_root='Solution.dat')
        
        

        
    def run_pareto(self):
        """
        run iCSD multiple times while changing the weights to explore the L-curve
        """           
        self.pareto_weights = np.linspace(self.pareto_MinErr, self.pareto_MaxErr, self.pareto_nSteps)
        print('pareto weights are\n', self.pareto_weights)

        self.pareto_list_FitRes = []
        self.pareto_list_RegRes = []
        with PdfPages(self.path2save+ self.obs +'iCSD_pareto.pdf') as pdf:
            for self.wr in self.pareto_weights:              
                if self.wr== self.pareto_weights[0]:
                    if self.x0_prior==True:
                        self.estimateM0()
                    if self.x0_ini_guess==True:
                        print('initial x0 guess')
                        self.estimateM0()
            
                print('-' * 80)
                print('Pareto mode, running new iCSDwith regularizaiton weight: ', self.wr)
                self.prepare4iCSD()
                self.iCSD()
                if self.type=="2d":
                    self.plotCSD()
                else:
                    self.f = plotCSD3d(self.wr,self.coord,self.x.x,self.path2load,self.obs,self.knee,self.KneeWr,title=None,**kwargs)
                    # self.plotCSD3d_pyvista()

                pdf.savefig(self.f)
                plt.close(self.f)

                self.ResidualAnalysis()

            
            self.DetectKneePt()
            self.wr=float(self.pareto_weights[self.IdPtkneew])
            print('Knee detected for wr=' + str(self.wr))
            self._plotPareto_()
            pdf.savefig(self.p)
            plt.close(self.p)

            if self.knee==True:
               self._plot_knee_icsd()
               
    ### PLOT
    def _plot_knee_icsd(self):
        """ Plot CSD for the best regularisation parameter after L-curve automatic analysis using a knee-locator
        
        Parameters
        ------------
        self
        """
        self.KneeWr=self.wr
        # self.wr=float(self.pareto_weights[self.IdPtkneew])
        self.kn.plot_knee_normalized()
        # self.kn.savefig('Lc', dpi = 450)
        # pdf.savefig(self.p)
        # plt.close(self.p)
        self.run_single()
        #self.f.savefig('iK'+ str("{:02d}".format(int(ObsName.translate({ord(i): None for i in 'OW'})))), dpi = 450,
        #    bbox_inches='tight',pad_inches = 0)
        plt.show()
            


    # def _mkGrid_XI_YI_(self):
    #     """ grid for interpolation """
    #     Xm = np.linspace(min(self.coord_x), max(self.coord_x), 500)
    #     Ym = np.linspace(min(self.coord_y), max(self.coord_y), 500)
    #     self.XI, self.YI = np.meshgrid(Xm, Ym)


        
    def _plotPareto_(self):
        self.p, self.ax = plt.subplots()
        # self.p = plt.figure('pareto', figsize=(10,4))
        self.ax.annotate('Wr=' + str(int(self.wr)), xy=(float(np.asarray(self.pareto_list_FitRes)[self.IdPtkneew]), 
                                         float(np.asarray(self.pareto_list_RegRes)[self.IdPtkneew])), 
                                      xytext=(float(np.asarray(self.pareto_list_FitRes)[self.IdPtkneew])+max(self.pareto_list_FitRes)/3, 
                                         float(np.asarray(self.pareto_list_RegRes)[self.IdPtkneew])+max(self.pareto_list_RegRes)/3),
                                      arrowprops=dict(facecolor='black', shrink=0.05))
        plt.plot(float(np.asarray(self.pareto_list_FitRes)[self.IdPtkneew]), float(np.asarray(self.pareto_list_RegRes)[self.IdPtkneew]), 'og')
        plt.plot(self.pareto_list_FitRes, self.pareto_list_RegRes, 'or')
        ax = plt.gca()
        ax.tick_params(axis = 'both', which = 'both', direction = 'out')
        ax.grid()
        plt.xlabel('Residual')
        plt.ylabel('Roughness')
        plt.tight_layout()
        plt.savefig(self.path2save+'ParetoFront.png', dpi = 600)



    def writeFIT(self):
        print('write FIT')
        np.savetxt(self.path2load+'invertedR.txt', self.x.fun[:] + self.b_w[:])
        np.savetxt(self.path2load+'truemodelR.txt', self.b_w[:])
        # np.savetxt(self.path2load+'Ap.txt', self.A_w[:])
        np.savetxt(self.path2load+'b_w.txt', self.b_w[:])
        np.savetxt(self.path2load+'b_s.txt', self.b_s[:])
        # np.savetxt(self.path2load+'reg_b.txt', self.reg_b[:])
        # np.savetxt(self.path2load+'b_s.txt', self.b_s[:])
        if self.x0_prior==True:
            np.savetxt(self.path2load+'x0F1_sum.wtxt', self.x0)
        np.savetxt(self.path2load+'reg_A.txt',self.reg_A)


        
    def run_single(self):
        """Run a single inversion (unique regularisation weight)
        Equivalent to several steps::
            self.prepare4iCSD()
            self.plotCSD()
            self.RMSAnalysis()
            
        """
        
        self.prepare4iCSD()            
        self.iCSD()
        self.showResults()
        self.RMSAnalysis()
        self.writeFIT()
        self.f.savefig(self.path2save+'iCSD', dpi = 600)
        plt.show()
        

        
    def Invert(self,pareto=False):
        """Invert the voltage to current densities.
        
        Parameters
        ----------

        """
        if self.pareto==False:
            self.run_single()
        else:
             self.run_pareto()
       
        
    def _parseModelReg(self):
        """ Parse regularisation parameters before inversion
        """    
        if self.type=='2d': # 2D CASE -----------------------------------------
            if self.regMesh=='strc': # structure mesh of virtual sources
                if self.alphaSxy==True:
                    self.reg_Ax, self.reg_Ay = regularize_A_x_y(self.coord,self.alphaSx,self.alphaSy)    
                    if self.x0_prior==True:
                        self.regularize_smallnessX0() #add smallness regularisation
                    self.regularize_sum_AX0()
                else:
                    if self.x0_prior==True:
                        # print('not possible')
                        raise ValueError('#### dimensions of matrices do not agree - change regularisation types')
                    else:
                        self.regularize_A()
            else:
                self.reg_A = regularize_A_UnstructuredMesh3d(self.coord,self.nVRTe,self.k)
        else:              # 3D CASE -----------------------------------------
            if self.regMesh=='strc':
                if self.x0_prior==True:
                        self.regularize_A_x_y_z()
                        self.regularize_smallnessX0() #add smallness regularisation
                        self.regularize_sum_AX0()
                        # raise ValueError('### dimensions of matrices do not agree - change regularisation type')
                else:
                     self.reg_A = regularize_A_3d(self.nVRTe,self.coord) # working for structured mesh
            elif self.regMesh=='unstrc':
                self.reg_A = regularize_A_UnstructuredMesh3d(self.coord,self.nVRTe,self.k)
                if self.x0_prior==True:
                    self.regularize_smallnessX0() #add smallness regularisation
            elif self.alphaSxy==True:
                if self.x0_prior==True:
                    self.regularize_smallnessX0() #add smallness regularisation
                self.regularize_sum_AX0()
                
    def parseDataReg(self):
        """ Parse regularisation parameters before inversion
        """    

#%% DEFINE SURVEY container for observations and simulations     

    def CreateSurvey(self):
        """Data container for survey paramaters such as geometry file
        
        Parameters
        ----------

        """


    def CreateTimeLapseSurvey(self):
        """Description here
        
        Parameters
        ----------

        """
        

#%% INITIAL MODEL          

    def estimateM0(self,method='F1',show=True):
        self.icsd_init()  
        self.parseM0(method) # define the method to estimate M0
        self.labelsM0(method) # labeling for the plot
        if show == True:
            self.showEstimateM0()
        
    def parseM0(self,method):
        """ Parse M0 parameters
        """
        if method=='F1': 
            self.norm_F1, self.x0 = misfit_2_initialX0(self.A,self.b,self.inix0)
        elif method=='Pearson': 
            self.x0 = productmoment(self.A,self.b)
            
    def labelsM0(self,method):
        """ Parse graphical labels to plot estimate M0 model
        """
        if method=='F1':      
            self.physLabel= 'normed misfit F1'
        if method=='Pearson':      
            self.physLabel= 'Pearson r coeff' 
            
    def showEstimateM0(self):
        fig, ax = plt.subplots()
        if self.type=='2d':
            # f = Scatter2d(self.coord, self.x0, self.physLabel, self.path2load,self.obs)
            f = plotContour2d(self.coord,self.x0,path=self.path2load,retElec=None, sc=None)
        else:
            f = Scatter3d(self.coord, self.x0, self.physLabel, self.path2load,self.obs)
        fig.tight_layout()
        
    # def run_misfitF1(self):
    #     self.icsd_init()            
    #     self.norm_F1, self.x0 = misfit_2_initialX0(self.A,self.b,self.inix0)


    

#%% PLOT fcts     


    def showResults(self, ax=None, clim=None ,cmap='viridis_r',
                    plotElecs=False,sc=False,retElec=None,
                    mesh=None, gif3d=False, title=None, cut=False,**kwargs):
        """Show inverted model.
        
        Parameters
        ----------
        ax : Matplotlib.Axes, optional
            If specified, the graph will be plotted against this axis.
        clim : array, optional
            Minimum and Maximum value of the colorbar.
        cmap : str, optional
            Name of the Matplotlib colormap to use.
        plotElecs : bool, optional
            If `True` add to the ICSD plot measuring electrodes as points
        sc : bool, optional
            If `True` add to the ICSD plot remotes and injection electrodes as points
        retElec : array, optional
            Coordinates of the return electrode, format = x1,y1')
            If Not None add to the ICSD plot the return B electrode
        mesh :  str, optional
            Specify name of the vtk file
            If Not None add mesh3d.vtk to plot with the results of icsd (for 3d using pyvista)
        gif3d :  bool, optional
            If `True` record a gif using orbital function of pyvista
          title :  str, optional
            Specify inversion titlename to be add to the plot         
        """
        self.clim=clim
        self.plotElecs=plotElecs
        self.sc=sc
        self.retElec=retElec
        self.mesh_over=mesh
        self.gif3d=gif3d
        self.title=title
       
        if self.type=='2d':
            self.f = plotCSD2d(self.coord,self.x.x,self.b,self.b_w,self.x.fun,self.path2load,self.pareto,retElec=None, sc=None)
        else:
            self.f = plotCSD3d(self.wr,self.coord,self.x.x,self.path2load,self.obs,self.knee,self.KneeWr,title=None,**kwargs)
            if cut ==True:
                plotCSD3d_pyvistaSlice()
            else:
                plotCSD3d_pyvista(self.coord, path=self.path2load, filename_root='Solution.dat', **kwargs)
        return
        
#%% POST inversion analysis        

    def ResidualAnalysis(self):
        fitting_res = self.x.fun[0 : self.b.shape[0]]
        # constrain_res = self.x.fun[self.b.shape[0] + 1] / self.wc
        regularization_res = self.x.fun[self.b.shape[0] + 2 :] / self.wr # constrain not included in the reg function
        self.reg_sum = np.sum(np.square(regularization_res))
        self.fit_sum = np.sum(np.square(fitting_res))
        self.pareto_list_FitRes.append(self.fit_sum)
        self.pareto_list_RegRes.append(self.reg_sum)

    def RMSAnalysis(self):
        self.rms = np.sum(np.square(self.x.fun))/len(self.x.fun)
        print('RMS error='+ str(self.rms))


#%% TO DO         

# NEW FUNCTIONS to introduce 
# function penalized some virtual sources (split rhizotron)
# load_names for tl analysis
# Introduce error weigting using reciprocal error instead of constant or w = 1/sqrt(abs(obs))
# use the ERT mesh as mesh for virtual position of current sources

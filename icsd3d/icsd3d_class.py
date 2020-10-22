import os


import numpy as np
from scipy.optimize import lsq_linear, least_squares

from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.colors import LogNorm

from kneed import KneeLocator

import matplotlib.pyplot as plt
import pyvista as pv

from plotters.mpl_plot import (plotCSD3d, plotCSD2d, 
                                scatter3d, scatter2d, 
                                plotContour2d, plotPareto, 
                                plot_knee_icsd,labels)

from plotters.pv_plot import plotCSD3d_pv
from inversion.regularize import *
from inversion.priorM0 import *
from inversion.solve import *
from importers.read import *
from exporters.save import * 
from gridder.mkgrid import mkGrid_XI_YI 

class iCSD3d_Class():
    
    
    """ 
    Create a icsd inversion object.
    """
    
    def __init__(self,dirName):
        
        self.df = None # main dataframe
        self.surveys = []


        # input files directory
        self.dirName = dirName
        
        # survey files to import (default names)
        self.type='2d' # or '3d'
        self.sim='VRTeSim.txt'
        self.obs='ObsData.txt'
        self.coord_file='VRTeCoord.txt' # coordinates of the virtual sources
        
        # inversion default parameters
        self.wr=25 #weight regularization
        self.wc=10000 #current conservation constrain, sum current fractions = 1
        self.pareto=False #if True run many icsd to explore the Pareto front
        self.pareto_MinErr=0.001
        self.pareto_MaxErr=1
        self.pareto_nSteps=10
        self.k=4  # For each point, find the k closest current sources
        self.TL=False # Time lapse inversion (see slides VIU Venice to implement)
        self.x0_prior=False #  relative smallness regularization as a prior criterion for the inversion; i.e the algorithm minimizes ||m−m0||2
        self.x0_ini_guess=False # initial guess
        self.knee=False # L-curve knee automatic detection
        self.KneeWr=[]
        
        # mesh type
        self.regMesh='strc' # strc or unstrc
        
        # Initial model
        self.alphax0=1 # weight on model smallness relative to m0
        self.inix0='cst' #  intial physical model: None or 'cst' if constant vector *0.1 
        self.method_m0=None # method to estimate intial physical model (F1= misfit from a single current source OR ProductMoment method)

        # spatial regularization
        self.alphaSxy=False # weight on model smoothness in z-direction 
        self.alphaSx=1 # weight on model smoothness in x-direction
        self.alphaSy=1 # weight on model smoothness in y-direction
        self.alphaSz=1 # weight on model smoothness in z-direction [TO IMPLEMENT]

        # data regularisation
        self.obs_err='sqrt' #const or sqrt or reciprocal [not yet implemented] - choose between constant weight and w = 1/sqrt(abs(obs))
        self.errRmin=1 #min R to which the err is applied before passing to constant error  
    
        # graphic options
        self.sc=[] #coordinates of the sources, format = x1,y1 x2,y2'
        self.retElec=None #coordinates of the return electrode, format = x1,y1')
        self.clim = []
        self.plotElecs=False 
        self.gif3d=False # create gif orbit
        self.title=None #
        self.mesh=None # mesh3d .vtk to plot with the results of icsd

        # IMPLEMENT obs_err based on reciprocal analysis i.e. estimated standard deviation of the data errors;                         % estimated standard deviation of the traveltime data errors

        # processing outputs [TO write]
        # self.models = []
        # self.rmses = []



    def icsd_init(self):
        """ these functions are called only once, even for pareto,
        as they do not depend on the regularization weight"""
        # create directories
        self.createdirs()
        
        # load virtual sources coordinates
        if self.type == '2d':
            self.coord_x, self.coord_y, self.coord = load_coord(self.path2load, self.coord_file, self.type)
        else:
            self.coord_x, self.coord_y, self.coord_z, self.coord = load_coord(self.path2load, self.coord_file, self.type)


        # load observations resistances b
        self.b = load_obs(self.path2load, self.obs)
        
        # load simulated resistances A (i.e. Green function)
        self.A = load_sim(self.path2load, self.sim)
        
        # load observations electrode coordinates
        if self.plotElecs==True:
            self.RemLineNb, self.Injection, self.coordE, self.pointsE= load_geom(self.path2load) # geometry file containing electrodes position includinf remotes 
        
        # check vector sizes
        self.nVRTe = check_nVRTe(self.A,self.b,self.coord)
        
        
        # reshape A vector into matrix A of nVRTE collumns
        self.A = reshape_A(self.A,self.nVRTe)
        
        # define mode to weights the data (const, sqrt or reciprocal)     
        self.obs_w= obs_w_f(self.obs_err,self.b,self.errRmin,sd_rec=None)
        
        # constrain (curent conservation)
        self.con_A = con_A_f(self.A)
        self.con_b = con_b_f(self.b)
        self.con_w = con_w_f(self.wc)
        
        # make grid mesh (??)
        # self.XI, self.YI= mkGrid_XI_YI(self.coord_x,self.coord_y)
        
        # append spatial regularization (add lines to the matrice)
        self._parseModelReg()
        # self.parseDataReg()
        self.reg_b= regularize_b(self.reg_A)
        
        # stack data, constrain, and regularization 
        self.A_s = stack_A(self.A, self.con_A, self.reg_A)
        self.b_s = stack_b(self.b, self.con_b, self.reg_b)
        
        
    ### mkdirs 
    def createdirs(self):
        self.path2load = self.dirName
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

    ### PREPARE FOR ICSD

    def prepare4iCSD(self):
        """ this function is called for each weight, keep them separated for pareto
        # 1 - Estimate M0 initial if smallness guided inversion or initial guess 
        # 2-  Create vector with weights, the length is determined by the number of regul rows in A

        """
        # create regularization part of the weight matrix
               
        if (self.x0_ini_guess==True or self.x0_prior==True):
            self._estimateM0_()
            
        
        if self.x0_prior==True: # if relative smallness 
            self.reg_w_0_b, self.reg_w_0_A = regularize_w(self.reg_A,self.wr,self.x0_prior,x0=self.x0)
            
            # stack weight matrix
            self.W_s_A, self.W_s_b = stack_w(self.obs_w, self.con_w, self.x0_prior, reg_w_0_A=self.reg_w_0_A, reg_w_0_b=self.reg_w_0_b)


            # apply weight to A and b (data weigth, constrainsts weight and regularisation weigth)
            self.A_w = weight_A(self.x0_prior,self.A_s,W_s_A=self.W_s_A)
            self.b_w = weight_b(self.x0_prior,self.b_s,W_s_b=self.W_s_b)

        else:
            self.reg_w = regularize_w(self.reg_A,self.wr,self.x0_prior)
            
            self.W_s = stack_w(self.obs_w, self.con_w, self.x0_prior, reg_w=self.reg_w)
            self.A_w = weight_A(self.x0_prior,self.A_s,W_s=self.W_s)
            self.b_w = weight_b(self.x0_prior,self.b_s,W_s=self.W_s)  

    def _parseModelReg(self):
        """ Parse regularisation parameters before inversion
        """
        
        # 2D CASE -----------------------------------------
        if self.type=='2d': 
            if self.regMesh=='strc': # structured (grid) mesh of virtual sources
                if self.alphaSxy==True: # spatial smoothing with coeff  alphaSx and alphaSy
                    self.reg_Ax, self.reg_Ay = regularize_A_x_y(self.coord,self.alphaSx,self.alphaSy)    
                else:
                    if self.x0_prior==True:
                        raise ValueError('#### dimensions of matrices do not agree - change regularisation types')
                    else:
                        self.reg_A = regularize_A(self.coord,self.nVRTe)
            # else: # 3d mesh
            #     self.reg_A = regularize_A_UnstructuredMesh3d(self.coord,self.nVRTe,self.k)
        
        
        # 3D CASE -----------------------------------------
        else:              
            if self.regMesh=='strc':
                if self.x0_prior==True:
                    # self.regularize_A_x_y_z()
                    raise ValueError('### Not yet ready to work with m0 and 3d grid mesh, regularize_A_x_y_z need to be tested')
                else:
                     self.reg_A = regularize_A_3d(self.nVRTe,self.coord) # working for structured mesh
            
            elif self.regMesh=='unstrc':
                self.reg_A = regularize_A_UnstructuredMesh3d(self.coord,self.nVRTe,self.k)
            
                    
        if self.alphaSxy==True: # anisotropic smoothing
            self.reg_smallx0 = ponderate_smallnessX0(self.alphaSxy,self.alphax0,reg_Ax=self.reg_Ax)
            self.reg_A = sum_smallness_smoothness(self.alphaSxy,self.x0_prior,reg_smallx0=self.reg_smallx0, reg_Ax=self.reg_Ax, reg_Ay=self.reg_Ay)
        else:
            self.reg_smallx0 = ponderate_smallnessX0(self.alphaSxy,self.alphax0,reg_A=self.reg_A)
            self.reg_A = sum_smallness_smoothness(self.alphaSxy,self.x0_prior,reg_smallx0=self.reg_smallx0,reg_A=self.reg_A)
            
        return self.reg_A
                
                
    def _parseDataReg(self):
        """ Parse data regularisation parameters before inversion
        """    

    def _parseSolver(self):
        """ Parse solver used during inversion
        """    
        
    def run_single(self):
        """Run a single inversion (unique regularisation weight)
        Equivalent to several steps::
            self.prepare4iCSD()
            self.plotCSD()
            self.misfit()
            
        """
        self.icsd_init()
        self.prepare4iCSD()
        if (self.x0_ini_guess == True or self.x0_prior == True):
            self.x = iCSD(self.x0_ini_guess,self.A_w,self.b_w,self.type,self.coord,self.path2load,x0=self.x0)
        else:
            self.x = iCSD(self.x0_ini_guess,self.A_w,self.b_w,self.type,self.coord,self.path2load)

        ax, f = self.showResults()
        # self.RMSAnalysis()
        self.misfit()
        print(f)
        print(self.path2save)
        f.savefig(self.path2save+'iCSD', dpi = 600)
        plt.tight_layout()
        plt.show()
        
        # add result to container
        # _add_to_container(self.x)
        
        # plt.close(f)
        # self.ModelResolution()
        
    def run_pareto(self):
        """
        run iCSD multiple times while changing the weights to explore the L-curve
        """           
        self.icsd_init()
        self.pareto_weights = np.linspace(self.pareto_MinErr, self.pareto_MaxErr, self.pareto_nSteps)
        print('pareto weights are\n', self.pareto_weights)

        self.pareto_list_FitRes = []
        self.pareto_list_RegRes = []
        
        
        with PdfPages(self.path2save+ self.obs +'iCSD_pareto.pdf') as pdf:
            
            # LOOP ON REG WEIGHTS
            for self.wr in self.pareto_weights:              
              
                
                # RUN ICSD
                self.prepare4iCSD()
                if (self.x0_ini_guess == True or self.x0_prior == True):
                    self.x = iCSD(self.x0_ini_guess,self.A_w,self.b_w,
                                  self.type,
                                  self.coord,
                                  self.path2load,
                                  x0=self.x0)
                else:
                    self.x = iCSD(self.x0_ini_guess,
                                  self.A_w,self.b_w,
                                  self.type,self.coord,
                                  self.path2load)

                # PLOT ICSD
                if self.type=="2d":
                    self.f, _ = plotCSD2d(self.coord,self.x.x,self.b,self.b_w,
                                       self.x.fun,
                                       self.path2load,
                                       self.pareto,
                                       retElec=None, sc=None, 
                                       title_wr=self.wr)
                else:
                    self.f, _ = plotCSD3d(self.wr,self.coord,self.x.x,
                                       self.path2load,
                                       self.obs,
                                       self.knee,self.KneeWr,
                                       title=None)

                pdf.savefig(self.f)
                # print(self.f)
                plt.show()
                plt.close(self.f)

                self.residualAnalysis()
                self.misfit()
            
            # Detect knee point wr of the L-curve
            self.detectKneePt()
            self.wr=float(self.pareto_weights[self.IdPtkneew])
            print('Knee detected for wr=' + str(self.wr))
            
            # Plot the L-curve
            self.p, self.ax= plotPareto(self.wr,self.pareto_list_FitRes,self.pareto_list_RegRes,self.IdPtkneew, self.path2save)
            pdf.savefig(self.p)
            # plt.close(self.p)

            # Plot knee curve
            if self.knee==True:
               plot_knee_icsd(self.wr,self.kn)
               self.run_single()


    def detectKneePt(self):
        self.kn = KneeLocator(self.pareto_list_FitRes,self.pareto_list_RegRes, 
                         curve='convex', direction='decreasing')
        self.IdPtkneew= np.where(self.kn.knee==self.pareto_list_FitRes)[0]
        self.pareto_weights[self.IdPtkneew]
        if len(self.IdPtkneew)<1:
             self.IdPtkneew=1
             print('No knee detection possible, put 1 as default')
#            raise ValueError('No knee detection possible')
               
        
    def invert(self,**kwargs):
        """Invert voltage to current densities.
        
        Parameters
        ----------
        * x0_prior:
            If True, relative smallness x0_prior model
        * regMesh:'unstrc' or 'strc'
            Specify mesh type 
        * wr:
            Specify regularisation weight
        * pareto :
            If True, run L-curve analysis.
        * pareto_MinErr, pareto_MaxErr, pareto_nSteps: 
            Specify pareto paramters
            
        Returns:
    
        * x : 
            pd.dataframe
            Data container
        """
        for (k, v) in kwargs.items():
            setattr(self, k, v)
            print("{0} = {1}".format(k, v))
         
        if self.pareto==False:
             self.run_single()
        else:
             self.run_pareto()
             
        return self.x
       
        


#%% DEFINE SURVEY container for observations and simulations     

        
    def CreateSurvey(self):
        """Data container using panda dataframe for survey paramaters such as geometry file
        
        Parameters
        ----------

        """
        survey = iCSD3d_Class()
    
        # set attribute according to the first survey
        if len(self.surveys) == 0:
            self.icsd_init()
            self.surveys.append(survey)
        else: # check we have the same configuration than other survey
            check = [a == b for a,b, in zip(self.nVRTe, survey.nVRTe)]
            if all(check) is True:
                self.surveys.append(survey)



    def CreateTimeLapseSurvey(self):
        """Import multiple surveys.
        
        Parameters
        ----------
        fnames : list of str
            List of file to be parsed or directory where the files are.
        targetProjection : str, optional
            If specified, a conversion from NMEA string in 'Latitude' and 'Longitude'
            columns will be performed according to EPSG code: e.g. 'EPSG:27700'.
        """
        if isinstance(fnames, list): # it's a list of filename
            if len(fnames) < 2:
                raise ValueError('at least two files needed for timelapse inversion')
        else: # it's a directory and we import all the files inside
            if os.path.isdir(fnames):
                fnames = [os.path.join(fnames, f) for f in np.sort(os.listdir(fnames)) if f[0] != '.']
                # this filter out hidden file as well
            else:
                raise ValueError('dirname should be a directory path or a list of filenames')
        if self.projection is not None:
            targetProjection = self.projection
        for fname in fnames:
            self.createSurvey(fname, targetProjection=targetProjection)


    def _add_to_container(self, df):
        """Add a given DataFrame to the container"""
            
        

#%% INITIAL MODEL          
    def _estimateM0_(self,show=True):
        # self.icsd_init()
        # self.parseM0(self.method_m0) # define the method to estimate M0
        # self.physLabel=labels(self.method_m0) # lgd labeling for the plot
        # if show == True:
        #     self.showEstimateM0()
        self.estimateM0(method_m0=self.method_m0)

    def _parseM0_(self,method_m0):
        """ Parse M0 parameters
        """
        if self.method_m0 is not None:
            if method_m0=='F1': 
                self.norm_F1, self.x0 = misfitF1_2_initialX0(self.A,self.b)
            elif method_m0=='Pearson': 
                self.x0 = product_moment(self.A,self.b)
        elif self.inix0 is not None:
            if self.inix0=='cst':
                self.x0=np.ones(self.b.shape)*0.1
                
        return self.x0
            
    def estimateM0(self,method_m0='F1',show=True, ax=None):
        """Estimate initial M0 model for constrainst inversion
        
        Parameters
        ----------
        * method_m0: 'F1' OR 'Pearson'
            Specify the method to estimate the initial M0
        * show:
            Show M0 misfit 
        * ax:
            Specify plot axis           
        Returns:
    
        """
        self.icsd_init()
        m0 = self._parseM0_(method_m0) # define the method to estimate M0
        self.physLabel=labels(method_m0) # lgd labeling for the plot
        if show == True:
            self.showEstimateM0(ax=ax)
            
        return m0
        

            

    def showEstimateM0(self,ax=None):
        # fig, ax = plt.subplots()
        if self.type=='2d':
            f = plotContour2d(self.coord,self.x0,self.physLabel,path=self.path2load,retElec=None, sc=None)
        else:
            f = scatter3d(self.coord, self.x0, self.physLabel, self.path2load, self.obs, ax=ax)
        plt.tight_layout()
        plt.show()
        # plt.close(f)

   

#%% PLOT fcts     


    def showResults(self, ax=None, data=None, clim=None ,cmap='viridis_r',
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
        self.retElec=retElec  # return electrode
        self.mesh_over=mesh # plot mesh over
        self.gif3d=gif3d # animated
        self.title=title # plot title
        # print(ax)
        if data is None:
            data = self.x.x
       
        if self.type=='2d':
            f = plotCSD2d(self.coord,data,self.b,self.b_w,self.x.fun,self.path2load,self.pareto,retElec=None, sc=None, ax=ax, title_wr=self.wr)
        else:
            f = plotCSD3d(self.wr,self.coord,data,self.path2load,self.obs,self.knee,self.KneeWr,ax=ax,title=None,**kwargs)
            if cut ==True:
                plotCSD3d_pyvistaSlice()
            else:
                plotCSD3d_pv(self.coord, path=self.path2load, filename_root='Solution.dat', 
                             knee = self.knee, wr = self.wr, KneeWr = self.KneeWr, 
                             mesh=self.mesh, plotElecs=plotElecs, gif3d = self.gif3d, **kwargs)
        return ax, f
        
#%% POST inversion analysis        

    def residualAnalysis(self):
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


    def misfit(self):
        """ 
        
        Parameters
        ----------
        
        """   
        residuals =  self.x.fun[0 : self.b.shape[0]] - self.b
        val = np.linalg.norm(self.obs_w*residuals)

        # val = np.sum(self.obs_w*(residuals**2))
        
        misfit = val*self.wr
        # print(self.b_w[0 : self.b.shape[0]])
        # print('RMS error='+ str(misfit))

        return misfit


    def modelResolution(self,sol=None,coord=None, jacMi=None):

        if sol is None:
            sol = self.x

        if coord is None:
            coord = self.coord
            
        J = sol.jac[self.b.shape[0] + 1:,jacMi]  

        # writeFIT(self.path2load,self.x0_prior,self.x.fun,self.b_w,self.b_s,self.x0,self.reg_A,J,saveall=False)

        if self.type == '2d':
            f = plotContour2d(coord,J,'Jac matrice',path=self.path2load,retElec=None, sc=None,jac=jacMi)
        else:
            f = scatter3d(coord,J,'Jac matrice', self.path2load, self.obs)
        plt.tight_layout()
        plt.show()
        plt.close(f)

#%% TO DO         

# NEW FUNCTIONS to introduce 
# function penalized some virtual sources (split rhizotron)
# load_names for tl analysis
# Introduce error weigting using reciprocal error instead of constant or w = 1/sqrt(abs(obs))
# use the ERT mesh as mesh for virtual position of current sources

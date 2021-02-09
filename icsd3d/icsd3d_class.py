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

from copy import deepcopy


class iCSD3d(object):
    
    
    """ 
    Create a icsd inversion object.
    """
    
    def __init__(self,dirName):
        
        self.surveys = []
        self.df = None # main dataframe containing resistance values (TDIP data in progress)


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
        self.x0_prior=False #  relative smallness regularization as a prior criterion for the inversion; i.e the algorithm minimizes ||m−m0||2
        self.x0_ini_guess=False # initial guess
        self.knee=False # L-curve knee automatic detection
        self.KneeWr=[]
        
        self.logTrans=False
            # time lapse option
        # self.iTimeLapse = False # to enable timelapse inversion [Not implemented yet]
        self.TL=False # Time lapse inversion (see slides VIU Venice to implement)

        # tdip option
        self.TDIP_flag=False # import tdip data


        # mesh type
        self.regMesh='strc' # strc or unstrc
        
        # Initial model
        self.alphax0=1 # weight on model smallness relative to m0
        self.inix0='cst' #  intial physical model: None or 'cst' if constant vector *0.1 
        self.method_m0='F1' # method to estimate intial physical model (F1= misfit from a single current source OR ProductMoment method)

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
        # print(len(self.surveys))
        if self.surveys:
            print('Existing survey')
        # else:
        #     print('Create a new survey')
        #     print(self.logTrans)
        #     survey = deepcopy(self)
        #     self.surveys.append(survey)

            
    def icsd_init(self,survey):
        """ these functions are called only once, even for pareto,
        as they do not depend on the regularization weight"""
        # create directories          
            
        print('initiation ICSD')
        self.createdirs(survey)
        
        # load virtual sources coordinates
        if self.type == '2d':
            self.coord_x, self.coord_y, survey.coord = load_coord(survey.path2load, self.coord_file, self.type)
        else:
            self.coord_x, self.coord_y, self.coord_z, survey.coord = load_coord(self.surveys[0], self.coord_file, self.type)

            
        if self.TDIP_flag==True:
            # load observations resistances b
            #self.b = load_obs(self.path2load, self.obs)
            # load simulated resistances A (i.e. Green function)
            self.A = load_sim(survey.path2load, survey.sim)          
            survey.b = survey.b[0:37]
            print('icsd_init TDIP survey- reshape *.bin to fit with sim')
        else:
            #if len(self.surveys)==1:
            #    print('single file')
                # load observations resistances b
            survey.b = load_obs(survey.path2load, survey.obs)
            # load simulated resistances A (i.e. Green function)
            self.A = load_sim(survey.path2load, survey.sim)

        print('log transformation: ' + str(survey.logTrans))
        #Log transformation before inversion
        if survey.logTrans==True:
            print('log transformation')
            print(survey.b[0])
            survey.b = np.log(survey.b) 
            self.A  = np.log(self.A) 
            print(survey.b[0])

           
        # load observations electrode coordinates
        if self.plotElecs==True:
            survey.RemLineNb, survey.Injection, survey.coordE, survey.pointsE= load_geom(self.path2load) # geometry file containing electrodes position includinf remotes 
        
        # check vector sizes
        survey.nVRTe = check_nVRTe(self.A,survey.b,survey.coord)
        
        
        # reshape A vector into matrix A of nVRTE collumns
        survey.A = reshape_A(self.A,survey.nVRTe)
        
        # define mode to weights the data (const, sqrt or reciprocal)     
        survey.obs_w= obs_w_f(self.obs_err,survey.b,self.errRmin,sd_rec=None)
        
        # constrain (curent conservation)
        survey.con_A = con_A_f(survey.A)
        survey.con_b = con_b_f(survey.b)
        survey.con_w = con_w_f(self.wc)
        
        # make grid mesh (??)
        # self.XI, self.YI= mkGrid_XI_YI(self.coord_x,self.coord_y)
        
        # append spatial regularization (add lines to the matrice)
        survey.reg_A = self._parseModelReg(survey)
        # self.parseDataReg()
        survey.reg_b = regularize_b(self.reg_A)
        

        # if self.TDIP_flag==False:
        # stack data, constrain, and regularization 
        survey.A_s = stack_A(survey.A, survey.con_A, survey.reg_A)
        survey.b_s = stack_b(survey.b, survey.con_b, survey.reg_b)
        # else:
        #     # stack data, constrain, and regularization 
        #     survey.A_s = stack_A(survey.A, survey.con_A, survey.reg_A)
        #     survey.b_s = stack_b(survey.b[:,1], survey.con_b, survey.reg_b)
        
        
    ### mkdirs 
    def createdirs(self,survey):
        survey.path2load = self.dirName
        survey.path2save= survey.path2load + 'figs'
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
        # 2-  Create vector with weights-->
        (data weigth, constrainsts weight and regularisation weigth)

        """
        # create regularization part of the weight matrix
        # for i, survey in enumerate(self.surveys):
        #     print('prepare4iCSD')
        #     print(i, survey)
        #     print('survey.x0_ini_guess')
        #     print(survey.x0_ini_guess)
        #     print('survey.reg_A')
        #     print(survey.reg_A)

        if (self.x0_ini_guess==True or self.x0_prior==True):
            self._estimateM0_()
            
        
        if self.x0_prior==True: # if relative smallness 
            self.reg_w_0_b, self.reg_w_0_A = regularize_w(self.surveys[0].reg_A,
                                                          self.wr,
                                                          self.x0_prior,
                                                          x0=self.x0)
            
            # stack weight matrix
            self.W_s_A, self.W_s_b = stack_w(self.surveys[0].obs_w, 
                                             self.surveys[0].con_w, 
                                             self.x0_prior, 
                                             reg_w_0_A=self.reg_w_0_A, 
                                             reg_w_0_b=self.reg_w_0_b)

            # apply weight to A and b (data weigth, constrainsts weight and regularisation weigth)
            self.A_w = weight_A(self.x0_prior,self.surveys[0].A_s,W_s_A=self.W_s_A)
            self.b_w = weight_b(self.x0_prior,self.surveys[0].b_s,W_s_b=self.W_s_b)

        else:
            self.reg_w = regularize_w(self.surveys[0].reg_A,
                                      self.wr,
                                      self.x0_prior)
            
            self.W_s = stack_w(self.surveys[0].obs_w, 
                               self.surveys[0].con_w, 
                               self.x0_prior, 
                               reg_w=self.reg_w)
            self.A_w = weight_A(self.x0_prior,
                                self.surveys[0].A_s,
                                W_s=self.W_s)
            self.b_w = weight_b(self.x0_prior,
                                self.surveys[0].b_s,
                                W_s=self.W_s)  

    def _parseModelReg(self,survey):
        """ Parse regularisation smoothing and 
        prior constrainst parameters before inversion
        """
        
        # 2D CASE -----------------------------------------
        if self.type=='2d': 
            if self.regMesh=='strc': # structured (grid) mesh of virtual sources
                if self.alphaSxy==True: # spatial smoothing with coeff  alphaSx and alphaSy
                    self.reg_Ax, self.reg_Ay = regularize_A_x_y(self.surveys[0].coord,self.alphaSx,self.alphaSy)    
                else:
                    if self.x0_prior==True:
                        raise ValueError('#### dimensions of matrices do not agree - change regularisation types')
                    else:
                        self.reg_A = regularize_A(survey.coord,
                                                  survey.nVRTe)
            else: # 3d mesh
                 self.reg_A = regularize_A_UnstructuredMesh3d(self.surveys[0].coord,
                                                              self.nVRTe,self.k)
        
        
        # 3D CASE -----------------------------------------
        else:              
            if self.regMesh=='strc':
                if self.x0_prior==True:
                    # self.regularize_A_x_y_z()
                    raise ValueError('### Not yet ready to work with m0 and 3d grid mesh, regularize_A_x_y_z need to be tested')
                else:
                     self.reg_A = regularize_A_3d(self.nVRTe,
                                                  self.surveys[0].coord) # working for structured mesh
            
            elif self.regMesh=='unstrc':
                self.reg_A = regularize_A_UnstructuredMesh3d(self.surveys[0].coord,
                                                             self.nVRTe,self.k)
            
                    
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
        
    def run_single_TL(self):
        """Run a time-lapse inversion (unique regularisation weight)
        Equivalent to several steps::
            self.prepare4iCSD()
            self.plotCSD()
            self.misfit()
            
        """
        
    def run_single(self,index=0,showfig=False):
        """Run a single inversion (unique regularisation weight)
        Equivalent to several steps::
            self.prepare4iCSD()
            self.plotCSD()
            self.misfit()
            
        """
        if not hasattr(self.surveys[0], 'A'):
            self.icsd_init(self.surveys[0])
            
        self.prepare4iCSD()
        
        # constrainsted inversion
        if (self.x0_ini_guess == True or self.x0_prior == True):
            self.x = iCSD(self.x0_ini_guess,
                          self.A_w,self.b_w,
                          self.type,
                          self.surveys[0].coord,
                          self.surveys[0].path2load,
                          x0=self.x0)
            
        # UNconstrainsted inversion
        else:
            self.x = iCSD(self.x0_ini_guess,
                          self.A_w,self.b_w,
                          self.type,
                          self.surveys[0].coord,
                          self.surveys[0].path2load)
        if showfig == True:
            ax, f = self.showResults(index=index)
            # self.RMSAnalysis()
            self.misfit()
            # print(f)
            # print(self.path2save)
            # f.savefig(self.path2save+'iCSD', dpi = 600)
            plt.tight_layout()
            plt.show()
            plt.close(self.f)
        # add result to container
        # _add_to_container(self.x)
        
        # plt.close(f)
        # self.ModelResolution()
        
    def run_pareto(self):
        """
        run iCSD multiple times while changing the weights to explore the L-curve
        """ 
        #self.icsd_init(self.surveys[0])
        self.pareto_weights = np.linspace(self.pareto_MinErr, self.pareto_MaxErr, self.pareto_nSteps)
        print('pareto weights are\n', self.pareto_weights)

        self.pareto_list_FitRes = []
        self.pareto_list_RegRes = []
        
        
        with PdfPages(self.surveys[0].path2save+ self.obs +'iCSD_pareto.pdf') as pdf:
            
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
                                  self.type,self.surveys[0].coord,
                                  self.surveys[0].path2load)

                # PLOT ICSD
                if self.type=="2d":
                    self.f, _ = plotCSD2d(self.surveys[0].coord,self.x.x,
                                          self.surveys[0].b,
                                          self.b_w,
                                       self.x.fun,
                                       self.surveys[0].path2load,
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
            self.p, self.ax= plotPareto(self.wr,
                                        self.pareto_list_FitRes,
                                        self.pareto_list_RegRes,
                                        self.IdPtkneew, 
                                        self.surveys[0].path2save)
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
            if len(self.surveys)<2:
                self.run_single()
                return self.x

            else:
                if self.TL==True: # flag for time-lapse inversion
                    self.run_single_TL()
                else:
                    for i, survey in enumerate(self.surveys):
                        survey.run_single(index=i)               
        else:
             self.run_pareto()
             
       
        


#%% DEFINE SURVEY container for observations and simulations     

        
    def createSurvey(self,fname_obs,fname_sim,append=False):
        """ Create an instance of iCSD3d and return
        a survey object.
        
        Parameters
        ----------
        fname_obs : list of str
            File to be parsed.
        fname_sim : list of str
            File to be parsed.
        """
        # set attribute according to the first survey
        if len(self.surveys) == 0:
            print('no existing survey')
            # survey = iCSD3d(dirName=self.dirName)
            survey = deepcopy(self)
            survey.obs=fname_obs
            survey.sim=fname_sim
            survey.icsd_init(survey)
            self.surveys.append(survey)       
        else: # check we have the same configuration than other survey
            if append: 
                print('append new survey')
                survey = iCSD3d(dirName=self.dirName)
                survey.obs=fname_obs
                survey.sim=fname_sim
                survey.icsd_init(survey)
                # check = [a == b for a,b, in zip(self.nVRTe, self.surveys[len(surveys-1)].nVRTe)]
                # if all(check) is True:
                self.surveys.append(survey)
            else:
                print('Replace existing survey')
                survey = deepcopy(self)
                survey.obs=fname_obs
                survey.sim=fname_sim
                survey.icsd_init(survey)
                self.surveys[0]=survey       
        return self.surveys

    def createTimeLapseSurvey(self,fnames_obs, fnames_sim):
        """Import multiple surveys and return
        a survey object.
        
        Parameters
        ----------
        fnames_obs : list of str
            List of file to be parsed or directory where the files are.
        """
        self.surveys = []
        if isinstance(fnames_obs, list): # it's a list of filename
            if len(fnames_obs) < 2:
                raise ValueError('at least two files needed for timelapse inversion')
        else: # it's a directory and we import all the files inside
            if os.path.isdir(fnames_obs):
                fnames_obs = [os.path.join(fnames_obs, f) for f in np.sort(os.listdir(fnames_obs)) if f[0] != '.']
                # this filter out hidden file as well
            else:
                raise ValueError('dirname should be a directory path or a list of filenames')
        for fnames_obs,fnames_sim  in zip(fnames_obs,fnames_sim):
            self.createSurvey(fnames_obs,fnames_sim,append=True)
            print(fnames_obs)

        return self.surveys


    def createTDIPSurvey(self,fname_obs,fname_sim):
        """create TDIP survey and return a survey object.
        
        Parameters
        ----------
        fname : *.bin file containing TDIP infos
        """
        # read TDIP file
        tdip_obs = loadTDIPSurvey(self.dirName + fname_obs) # *.data (pygimli format)
        tdip_sim = load_sim(self.dirName, fname_sim) #

        self.dataTDIP = tdip_obs.data
        self.Vs = tdip_obs.MA.transpose()
        self.gates =  tdip_obs.t

        # if len(self.surveys) == 0:

        
        # For now, transform each gates to several surveys as it is done for timelapse icsd
        self.surveys = []
        for i  in range(self.Vs.shape[1]):
            survey = iCSD3d(dirName=self.dirName)
            survey.TDIP_flag = True # activate tdip flag
            print(survey)
            survey.b=self.Vs[:,i]
            self.A=tdip_sim
            survey.icsd_init(survey)
            self.surveys.append(survey)     
            
        # survey = iCSD3d(dirName=self.dirName)
        # survey.icsd_init() 
        
        return self.surveys
    

    # def _add_to_container(self, df):
    #     """Add a given DataFrame to the container"""
            
        

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
        for i, survey in enumerate(self.surveys):           
            if method_m0 is not None:
                if method_m0=='F1': 
                    self.surveys[i].norm_F1, self.surveys[i].x0 = misfitF1_2_initialX0(survey.A,survey.b)
                elif method_m0=='Pearson': 
                    self.surveys[i].x0 = product_moment(survey.A,survey.b)
            elif self.inix0 is not None:
                if survey.inix0=='cst':
                    self.surveys[i].x0=np.ones(survey.b.shape)*0.1
        self.x0 = np.copy(self.surveys[i].x0)
        return self.surveys[i].x0
            
    def estimateM0(self,method_m0='F1',show=True, ax=None):
        """Estimate initial M0 model for constrainsted inversion
        
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
        for i, survey in enumerate(self.surveys):
            m0 = self._parseM0_(method_m0) # define the method to estimate M0
            self.surveys[i].physLabel=labels(method_m0) # lgd labeling for the plot
            if show == True:
                self.showEstimateM0(index=i,ax=ax)
            
        return m0
        
    def showEstimateM0(self,index=0,ax=None):
        """Show initial model m0 estimation.
        
        Parameters
        ----------
        index : int, optional
            Index of the survey to plot.
        ax : Matplotlib.Axes, optional
            If specified, the graph will be plotted against this axis.
        """

        # fig, ax = plt.subplots()
        if self.type=='2d':
            f = plotContour2d(self.surveys[index].coord,self.surveys[index].x0,
                              self.surveys[index].physLabel,
                              path=self.surveys[index].path2load,
                              retElec=None, sc=None,index=index,
                              ax=ax)
        else:
            f = scatter3d(self.coord, self.x0, self.physLabel, self.path2load, self.obs, ax=ax)
        plt.tight_layout()
        plt.show()
        # plt.close(f)

    def showCurrentStreamlines(self,index=0,ax=None):
        """Show Current Streamlines.
        
        Parameters
        ----------
        index : int, optional
            Index of the survey to plot.
        ax : Matplotlib.Axes, optional
            If specified, the graph will be plotted against this axis.
        """

        # fig, ax = plt.subplots()
        if self.type=='2d':
            f = current_streamlines(path)
        else:
            print('not yet implemented in 3d')
   

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
        
        index=0
        for key, value in kwargs.items():
            if key == 'index':
                index = value
 
                
        # print(ax)
        if data is None:
            data = self.x.x
       
        if self.type=='2d':
            f, ax = plotCSD2d(self.surveys[0].coord,
                              data,
                              self.surveys[0].b,
                              self.b_w,
                              self.x.fun,
                              self.surveys[0].path2load,
                              self.pareto,
                              retElec=None, sc=None,
                              ax=ax, title_wr=self.wr, index=index)
        else:
            f = plotCSD3d(self.wr,self.coord,data,
                          self.path2load,self.obs,self.knee,
                          self.KneeWr,ax=ax,title=None,**kwargs)
            if cut ==True:
                plotCSD3d_pyvistaSlice()
            else:
                plotCSD3d_pv(self.coord, path=self.path2load, 
                             filename_root='Solution.dat', 
                             knee = self.knee, wr = self.wr, 
                             KneeWr = self.KneeWr, 
                             mesh=self.mesh, plotElecs=plotElecs, 
                             gif3d = self.gif3d, **kwargs)
        return ax, f
        
#%% POST inversion analysis        

    def residualAnalysis(self):
        fitting_res = self.x.fun[0 : self.surveys[0].b.shape[0]]
        # constrain_res = self.x.fun[self.b.shape[0] + 1] / self.wc
        regularization_res = self.x.fun[self.surveys[0].b.shape[0] + 2 :] / self.wr # constrain not included in the reg function
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
        residuals =  self.x.fun[0 : self.surveys[0].b.shape[0]] - self.surveys[0].b
        val = np.linalg.norm(self.surveys[0].obs_w*residuals)

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

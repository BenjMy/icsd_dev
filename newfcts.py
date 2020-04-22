# -*- coding: utf-8 -*-
"""
Created on Tue Apr 21 11:58:59 2020

@author: Benjamin
"""


        plotcontour=True
        if plotcontour==True:
            
            nx, ny = (50, 50)
            x = np.linspace(min(self.coordE[:self.RemLineNb,1]), max(self.coordE[:self.RemLineNb,1]), nx)
            y = np.linspace(min(self.coordE[:self.RemLineNb,2]), max(self.coordE[:self.RemLineNb,2]), nx)
            xv, yv = np.meshgrid(x, y)
            
            # X, Y = np.meshgrid(sensors[:,0], sensors[:,1])
            print(len(self.coordE[:,2::]))
            print(len(self.b))
            
            print(self.coordE[:self.RemLineNb,2])
            print(self.b)
            
            z= griddata(self.coordE[5:,2::], self.b, (xv, yv), method='cubic')
 
            # fig = plt.figure()
            # ax = fig.gca(projection='3d')
            # X, Y, Z = axes3d.get_test_data(0.05)
            # ax.plot_surface(X, Y, Z, rstride=8, cstride=8, alpha=0.3)
            cs = plt.contourf(xv, yv, z, zdir='z', offset=-5, cmap=cm.coolwarm) #, norm=normalize)
            # cset = ax.contourf(X, Y, Z, zdir='z', offset=-5, cmap=cm.coolwarm)
            # cset = ax.contourf(X, Y, Z, zdir='x', offset=-40, cmap=cm.coolwarm)
            # cset = ax.contourf(X, Y, Z, zdir='y', offset=40, cmap=cm.coolwarm)
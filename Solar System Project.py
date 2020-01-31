#!/usr/bin/env python
# coding: utf-8

# In[23]:


"""Imported Files"""
import numpy as np
import random as rand
import pylab as plb
import matplotlib.animation as animation
import time
from scipy.spatial.distance import pdist,squareform
import math
plb.rcParams['figure.figsize'] = 8,6


# In[46]:


"""Constants"""
____MASS____ = 5.974e24        #Earth mass
____TIME____ = 2897756704      # 91.887262 Year
____DIST____ = 149.60e9        # 1 AU
____G____ = 1.     # in reduced units
dt = 0.0003
#with this choose the appropriate h or step is 0.0003 which is about 10 days.


# In[66]:


Grav_Potential = lambda m,r: -____G____*m/r   # All in reduced units
Grav_Force = lambda m1,m2,r: ____G____*m1*m2/(r**2)   # All in reduced units in the r^ direction


# In[68]:


class Planet:
    """Planet calss makes a planet. Notice that all units are reduced."""
    def __init__(self,name,color,mass,dist_from_sun,radius,speed,marker="o"):
        
        self.name = name
        self.color = color
        self.marker = marker
        self.mass = mass
        self.dist_from_sun = dist_from_sun
        self.radius = radius
        self.speed = speed
        self.X = self.dist_from_sun
        self.Y = 0.
        self.Vx = 0.
        self.Vy = self.speed
        self.F_X_0 = 0.
        self.F_Y_0 = 0.

class System:
    
    def __init__(self):
        self.Planets = []
        self.NoP = 0  # Refers to Number of Planets
        
    def add_planet(self,PLANET):
        self.Planets.append(PLANET)
        self.NoP += 1
    
    def set_coordinates(self):
        """Just Run this method once after all planets are included!"""
        
        self.Coord = np.zeros((self.NoP,4),'float32')
        
        for j in range(self.NoP):
            self.Coord[j,0] = self.Planets[j].X
            self.Coord[j,1] = self.Planets[j].Y
            self.Coord[j,2] = self.Planets[j].Vx
            self.Coord[j,3] = self.Planets[j].Vy

        
    def plot(self):
        for pl in SoSy.Planets:
            plb.scatter(pl.X,pl.Y,label=pl.name,color=pl.color,marker=pl.marker)
        plb.legend()
        
    def Action(self,C=1):
    
        for i in range(C):
            for pl in self.Planets:
                #global F_X_0
                #global F_Y_0

                #dots[:,0] += dots[:,2]*dt + 0.5*F_X_0*(dt**2)
                #dots[:,1] += dots[:,3]*dt + 0.5*F_Y_0*(dt**2)

                pl.X += pl.Vx*dt 
                pl.Y += pl.Vy*dt

                """
                cross_left = dots[:,0]<0.
                cross_right = dots[:,0]>__L__
                cross_top = dots[:,1]>__L__
                cross_bot = dots[:,1]<0.


                dots[cross_left,0] += __L__
                dots[cross_right,0] -= __L__

                dots[cross_top,1] -= __L__
                dots[cross_bot,1] += __L__
                """


                #F_X_1 , F_Y_1 = InterForces(dots,cutoff)

                #dots[:,2] += 0.5*(F_X_0+F_X_1)*(dt)
                #dots[:,3] += 0.5*(F_Y_0+F_Y_1)*(dt)

                #F_X_0 = F_X_1
                #F_Y_0 = F_Y_1
                
    
    def InterplanetaryForces(self):

        D = squareform(pdist(self.Coord[:,0:2]))
        
        ind1 , ind2 = np.where(D>0)
        unq = ind1 < ind2
        ind1 = ind1[unq] 
        ind2 = ind2[unq]

        force_matrix_x = np.zeros((self.NoP,self.NoP))                
        force_matrix_y = np.zeros((self.NoP,self.NoP))             

        for i1,i2 in zip(ind1,ind2): 

            r = (self.Coord[i2,0]-self.Coord[i1,0],self.Coord[i2,1]-self.Coord[i1,1])       


            theta_x = angle_between(r,(1,0))
            theta_y = angle_between(r,(0,1))

            force = Grav_Force(self.Planets[i1].mass,self.Planets[i2].mass,D[i1,i2])


            force_matrix_x[i1,i2] += -force*np.cos(theta_x)
            force_matrix_y[i1,i2] += -force*np.cos(theta_y)

            force_matrix_x[i2,i1] += +force*np.cos(theta_x)
            force_matrix_y[i2,i1] += +force*np.cos(theta_y)


        F_X = np.sum (force_matrix_x,axis=1)
        F_Y = np.sum (force_matrix_y,axis=1)

        return F_X,F_Y             


Sun = Planet("Sun",'yellow',333110.144,0.,4.654078e-3,0.)        
Earth = Planet("Earth",'blue',1.,1.,4.26337e-5,577.226937)
Mars = Planet("Mars",'red',0.107,1.524,2.27072e-5,466.817758)

SoSy = System ()

SoSy.add_planet(Sun)
SoSy.add_planet(Earth)
SoSy.add_planet(Mars)
SoSy.set_coordinates()

#SoSy.Action(C=10)

SoSy.plot()

print(SoSy.Coord)


# In[49]:


____TIME____/____DIST____*2.41e4 


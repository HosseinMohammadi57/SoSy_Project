#!/usr/bin/env python
# coding: utf-8

# In[1]:


"""Imported Files"""
import numpy as np
import random as rand
import pylab as plb
import matplotlib.animation as animation
import time
from scipy.spatial.distance import pdist,squareform
import math
plb.rcParams['figure.figsize'] = 8,6


# In[2]:


"""Constants"""
____MASS____ = 5.974e24   #Earth mass
____TIME____ = 86400      # 1 Day
____DIST____ = 149.60e9   # 1 AU
____G____ = 8.884e-10     # in reduced units


# In[3]:


Grav_Potential = lambda m,r: -____G____*m/r   # All in reduced units
Grav_Force = lambda m,r: ____G____*m/(r**2)   # All in reduced units in the r^ direction


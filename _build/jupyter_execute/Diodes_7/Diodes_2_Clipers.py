#!/usr/bin/env python
# coding: utf-8

# In[1]:


#%%writefile Diodes_7_Codes.py
#Library import statements

from skidl.pyspice import *
#can you say cheeky 
import PySpice as pspice
#becouse it's written by a kiwi you know
import lcapy as kiwi

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt

from IPython.display import YouTubeVideo, display

import traceback


# In[2]:


#from DC_1_Codes import get_skidl_spice_ref
#notebook specific loading control statements 
get_ipython().run_line_magic('matplotlib', 'inline')
#tool to log notebook internals
#https://github.com/jrjohansson/version_information
get_ipython().run_line_magic('load_ext', 'version_information')
get_ipython().run_line_magic('version_information', 'skidl, PySpice,lcapy, sympy, numpy, matplotlib, pandas, scipy')


# In[3]:


reset()
d=D(model='1N4148')
vs=V(dc_value=1@u_V)
vs['p', 'n']+=d['p'], gnd
d['n']+=gnd
circ=generate_netlist(libs='SpiceLib')
print(circ)


# In[4]:


reset()
q=BJT(model='2n2222a')
vdc = V(dc_value=5@u_V)     # 5V power supply.
rs = R(value=5@u_kOhm)      # Source resistor in series with sine wave input voltage.
rb = R(value=25@u_kOhm)     # Bias resistor from 5V to base of transistor.
rc = R(value=1@u_kOhm)      # Load resistor connected to collector of transistor.
vs = SINEV(amplitude=0.01@u_V, frequency=1@u_kHz)  # 1 KHz sine wave input source.
q['c', 'b', 'e'] += rc[1], rb[1], gnd  # Connect transistor CBE pins to load & bias resistors and ground.
vdc['p'] += rc[2], rb[2]    # Connect other end of load and bias resistors to power supply's positive terminal.
vdc['n'] += gnd             # Connect negative terminal of power supply to ground.
rs[1,2] += vs['p'], q['b']  # Connect source resistor from input source to base of transistor.
vs['n'] += gnd              # Connect negative terminal of input source to ground.


circ=generate_netlist(libs='SpiceLib')
print(circ)


# In[ ]:





#!/usr/bin/env python
# coding: utf-8

# Just Exploring

# In[1]:


from skidl.pyspice import *
#can you say cheeky 
import PySpice as pspice
#becouse it's written by a kiwi you know
import lcapy as kiwi

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sympy as sym


from IPython.display import YouTubeVideo, display

import traceback
import warnings


# In[2]:


from DC_1_Codes import get_skidl_spice_ref
from AC_2_Codes import *

sym.init_printing()

#notebook specific loading control statements 
get_ipython().run_line_magic('matplotlib', 'inline')
#tool to log notebook internals
#https://github.com/jrjohansson/version_information
get_ipython().run_line_magic('load_ext', 'version_information')
get_ipython().run_line_magic('version_information', 'skidl, PySpice,lcapy, sympy, numpy, matplotlib, pandas, scipy')


# In[3]:


#instatate the rc_lowpass filter to 
lowpassF=rc_lowpass(C_value=.1@u_uF, R_value=1@u_kOhm)
lowpassF.lcapy_self()


# In[4]:


reset()
#create the nets
net_in=Net('In'); net_out=Net('Out'); 

#create a 1V AC test source and attache to nets
vs=SINEV( ac_magnitude=1@u_V); vs['p', 'n']+=net_in, gnd

#attaceh term_0 to net_in and term_2 to net_out per scikit-rf convention all 
#other terminals are grounded
lowpassF.SKiDl(net_in, gnd, net_out, gnd)

circ=generate_netlist()
print(circ)


# In[5]:


sim=circ.simulator()


# In[6]:


res=sim.ac_sensitivity('V(Out)', 'lin', 10, 1@u_Hz, 10@u_MHz)


# In[7]:


res.elements


# In[18]:


res=sim.dc_sensitivity('V(Out)')


# In[19]:


res.branches


# In[20]:


res.elements


# In[21]:


res.internal_parameters


# In[22]:


res.nodes


# In[ ]:





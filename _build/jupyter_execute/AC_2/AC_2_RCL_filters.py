#!/usr/bin/env python
# coding: utf-8

# IN PROGRESS
# 
# TODO:
# - add sympy to subcirucit classes
# - add ac sim tool; try to do an inhertince from last section
# - add rl complimeint filters
# - create table of low-high pass vs rc & rl
# - add rlc resonce
# - add rlc filters
# - add filter desighn tool and filter cirucit generation tool
# - discues L->C passive conversion

# In[1]:


class ac_ease():
    def __init__(self, circ_netlist_obj):
        self.circ_netlist_obj=circ_netlist_obj
        self._build_table()
    
    def _build_table(self):
        self.fsweep_DF=pd.DataFrame(columns=['Start_freq', 'Stop_Freq', 'SamplingInc', 'StepType'])
        self.fsweep_DF.at[len(self.fsweep_DF)]=[.1@u_Hz, 120@u_GHz, 10, 'decade']

    def ac_sweep_setup(self, Start_freq, Stop_Freq, SamplingInc, StepType, display=False):
        """
        TODO:
            -add assertions and translations for StepType
            -add force change of 0Hz to .01Hz and warning for sweeps other then linear
            -add display action
        """
        self.fsweep_DF.at[0]=[Start_freq, Stop_Freq, SamplingInc, StepType]
    
    def _make_sim_control(self):
        pass
    
    def do_ac_sim(self):
        pass
    
    def record_ac_nodebranch(self):
        pass
        


# In[2]:


#%%writefile AC_2_Codes.py
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


# In[3]:


#notebook specific loading control statements 
get_ipython().run_line_magic('matplotlib', 'inline')
#tool to log notebook internals
#https://github.com/jrjohansson/version_information
get_ipython().run_line_magic('load_ext', 'version_information')
get_ipython().run_line_magic('version_information', 'skidl, PySpice,lcapy, sympy, numpy, matplotlib, pandas, scipy')


# # Basic Passive Filters

# ## The low pass filter from ALL ABOUT ELECTRONICS "RC Low Pass Filter Explained"  @~ 8:35min

# In[4]:


YouTubeVideo('_2L0l-E1Wx0', width=500, height=400, start=515)


# In[5]:


class rc_lowpass():
    def __init__(self, subcirc_ref=None, C_value=1@u_F, R_value=1@u_Ohm):
        """
        TODO:
            -add assertions
        """
        #add assertions
        self.subcirc_ref=subcirc_ref
        self.C_value=C_value
        self.R_value=R_value
        

    @subcircuit
    def SKiDl(self, term_0, term_1, term_2, term_3, return_elements=False):
        """
        Termanls are defined via:
        ```
       Left_Termanals - RC_Lowpass - Right_Termanals
                          +----+  
        Postive V_i        0-|0  2|-2     Postive V_o
        Negtive V_i        1-|1  3|-3     Negtive V_o
                          +----+
        ```
        """
        
        self.c=C(value=self.C_value)
        self.r=R(value=self.R_value)
        
        self.r[1, 2]+=term_0, self.c['p']
        self.c['p', 'n']+=term_2, term_3
        
        
        
        
        if return_elements:
            return self.c, self.r
    
    def draw_me(self):
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 1 1_1; right=2')
        self.schematic.add('W 0 0_1; right')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add(f'R 0_1 2_1; right, l=R{str(self.R_value)}')
        self.schematic.add(f'C 2_1 1_1; down, l=C{str(self.C_value)}')
        
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 1_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')

        
        self.schematic.draw()
        
        
    


# In[6]:


lowpassF=rc_lowpass(C_value=.1@u_uF, R_value=1@u_kOhm)
lowpassF.draw_me()


# In[7]:


reset()
net_1=Net('N1'); net_2=Net('N2'); 

vs=SINEV(amplitude=10@u_V, frequency=2@u_kHz)
vs['p', 'n']+=net_1, gnd
lowpassF.SKiDl(net_1, gnd, net_2, gnd)




circ=generate_netlist()
print(circ)


# In[8]:


sim=circ.simulator()


# In[9]:


ac_vals=sim.ac(start_frequency=.1@u_Hz, stop_frequency=1@u_MHz, number_of_points=10, variation='dec')


# In[10]:


v_i=ac_vals[node(net_1)]
v_o=ac_vals[node(net_2)]
f=ac_vals.frequency


# In[11]:


plt.semilogx(f, 20*np.log10(np.abs(v_o)))
plt.twinx()
plt.semilogx(f, np.rad2deg(np.unwrap(np.angle(v_o))), color='g')


# In[12]:


v_o_time=np.fft.fftshift(np.fft.ifft(v_o))
plt.plot(np.abs(v_o_time))


# In[13]:


#pz_vals=sim.pole_zero(node(net_1), node(gnd), node(net_2), node(gnd), 'vol', 'pz')


# ## Equivlint Low Pass RL filter 
# 
# The eqivlint RL filter to the above RC filter may be found via the equivlint time consitnc of the RC and RL implemntaiotn such that time constance must mach ie:
# $$RC=\tau_{RC}=\tau_{RL}=L/R$$

# In[14]:


(1e3)**2 *.1e-9


# In[15]:


class rl_lowpass():
    def __init__(self, subcirc_ref=None, L_value=1@u_H, R_value=1@u_Ohm):
        """
        TODO:
            -add assertions
        """
        #add assertions
        self.subcirc_ref=subcirc_ref
        self.L_value=L_value
        self.R_value=R_value
        

    @subcircuit
    def SKiDl(self, term_0, term_1, term_2, term_3, return_elements=False):
        """
        Termanls are defined via:
        ```
       Left_Termanals - RL_Lowpass - Right_Termanals
                          +----+  
        Postive V_i        0-|0  2|-2     Postive V_o
        Negtive V_i        1-|1  3|-3     Negtive V_o
                          +----+
        ```
        """
        
        self.l=L(value=self.L_value)
        self.r=R(value=self.R_value)
        
        self.l['p', 'n']+=term_0, term_2
        self.r[1, 2]+=self.l['n'], term_1
        
        
        
        
        if return_elements:
            return self.l, self.r
    
    def draw_me(self):
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 1 1_1; right=2')
        self.schematic.add('W 0 0_1; right')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add(f'L 0_1 2_1; right, l=L{str(self.L_value)}')
        self.schematic.add(f'R 2_1 1_1; down, l=R{str(self.R_value)}')
        
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 1_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')

        
        self.schematic.draw()
        
        
    


# In[16]:


rl_l=rl_lowpass(L_value=(1e3)**2 *.1e-9 @u_H, R_value=1e3)
rl_l.draw_me()


# In[17]:


reset()
net_1=Net('N1'); net_2=Net('N2'); 

vs=SINEV(amplitude=10@u_V, frequency=2@u_kHz)
vs['p', 'n']+=net_1, gnd
rl_l.SKiDl(net_1, gnd, net_2, gnd)




circ=generate_netlist()
print(circ)


# In[18]:


sim=circ.simulator()


# In[19]:


ac_vals=sim.ac(start_frequency=1@u_Hz, stop_frequency=2@u_MHz, number_of_points=10, variation='dec')


# In[20]:


v_i=ac_vals[node(net_1)]
v_o=ac_vals[node(net_2)]
f=ac_vals.frequency


# In[21]:


plt.semilogx(f, 20*np.log10(np.abs(v_o)))
plt.twinx()
plt.semilogx(f, np.rad2deg(np.unwrap(np.angle(v_o))), color='g')


# In[22]:


v_o_time=np.fft.fftshift(np.fft.ifft(v_o))
plt.plot(np.abs(v_o_time))


# ## The high pass filter from ALL ABOUT ELECTRONICS "RC High Pass Filter Explained"  @~ 7:57min

# In[23]:


YouTubeVideo('9Dx0b0ukNAM', width=500, height=400, start=477)


# In[24]:


class rc_highpass():
    def __init__(self, subcirc_ref=None, C_value=1@u_F, R_value=1@u_Ohm):
        """
        TODO:
            -add assertions
        """
        #add assertions
        self.subcirc_ref=subcirc_ref
        self.C_value=C_value
        self.R_value=R_value
        

    @subcircuit
    def SKiDl(self, term_0, term_1, term_2, term_3, return_elements=False):
        """
        Termanls are defined via:
        ```
       Left_Termanals - RC_Highpass - Right_Termanals
                          +----+  
        Postive V_i        0-|0  2|-2     Postive V_o
        Negtive V_i        1-|1  3|-3     Negtive V_o
                          +----+
        ```
        """
        
        self.c=C(value=self.C_value)
        self.r=R(value=self.R_value)
        
        self.c['p', 'n']+=term_0, term_2
        self.r[1, 2]+=self.c['n'], term_1
        
        
        
        
        if return_elements:
            return self.c, self.r
    
    def draw_me(self):
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 1 1_1; right=2')
        self.schematic.add('W 0 0_1; right')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add(f'C 0_1 2_1; right, l=C{str(self.C_value)}')
        self.schematic.add(f'R 2_1 1_1; down, l=R{str(self.R_value)}')
        
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 1_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')

        
        self.schematic.draw()
        
        
    


# In[25]:


highpassF=rc_highpass(C_value=1.5@u_nF, R_value=10@u_kOhm)
highpassF.draw_me()


# In[26]:


reset()
net_1=Net('N1'); net_2=Net('N2'); 

vs=SINEV(amplitude=10@u_V, frequency=10@u_kHz)
vs['p', 'n']+=net_1, gnd
highpassF.SKiDl(net_1, gnd, net_2, gnd)




circ=generate_netlist()
print(circ)


# In[27]:


sim=circ.simulator()


# In[28]:


ac_vals=sim.ac(start_frequency=1@u_Hz, stop_frequency=1@u_MHz, number_of_points=10, variation='dec')


# In[29]:


v_i=ac_vals[node(net_1)]
v_o=ac_vals[node(net_2)]
f=ac_vals.frequency


# In[30]:


plt.semilogx(f, 20*np.log10(np.abs(v_o)))
plt.twinx()
plt.semilogx(f, np.rad2deg(np.unwrap(np.angle(v_o))), color='g')


# In[31]:


v_o_time=np.fft.fftshift(np.fft.ifft(v_o))
plt.plot(np.abs(v_o_time))


# ## RL Highpass

# In[32]:


class rl_highpass():
    def __init__(self, subcirc_ref=None, L_value=1@u_H, R_value=1@u_Ohm):
        """
        TODO:
            -add assertions
        """
        #add assertions
        self.subcirc_ref=subcirc_ref
        self.L_value=L_value
        self.R_value=R_value
        

    @subcircuit
    def SKiDl(self, term_0, term_1, term_2, term_3, return_elements=False):
        """
        Termanls are defined via:
        ```
       Left_Termanals - RL_Highpass - Right_Termanals
                          +----+  
        Postive V_i        0-|0  2|-2     Postive V_o
        Negtive V_i        1-|1  3|-3     Negtive V_o
                          +----+
        ```
        """
        
        self.l=L(value=self.L_value)
        self.r=R(value=self.R_value)
        
        self.r[1, 2]+=term_0, self.l['p']
        self.l['p', 'n']+=term_2, term_3
        
        
        
        
        if return_elements:
            return self.l, self.r
    
    def draw_me(self):
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 1 1_1; right=2')
        self.schematic.add('W 0 0_1; right')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add(f'R 0_1 2_1; right, l=R{str(self.R_value)}')
        self.schematic.add(f'L 2_1 1_1; down, l=L{str(self.L_value)}')
        
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 1_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')

        
        self.schematic.draw()
        
        
    


# In[33]:


(10e3)**2 *1.5e-9


# In[34]:


rl_h=rl_highpass(L_value=(10e3)**2 *1.5e-9 @u_H, R_value=10@u_kOhm)
rl_h.draw_me()


# In[35]:


reset()
net_1=Net('N1'); net_2=Net('N2'); 

vs=SINEV(amplitude=10@u_V, frequency=10@u_kHz)
vs['p', 'n']+=net_1, gnd
rl_h.SKiDl(net_1, gnd, net_2, gnd)




circ=generate_netlist()
print(circ)


# In[36]:


sim=circ.simulator()


# In[37]:


ac_vals=sim.ac(start_frequency=1@u_Hz, stop_frequency=1@u_MHz, number_of_points=10, variation='dec')


# In[38]:


v_i=ac_vals[node(net_1)]
v_o=ac_vals[node(net_2)]
f=ac_vals.frequency


# In[39]:


plt.semilogx(f, 20*np.log10(np.abs(v_o)))
plt.twinx()
plt.semilogx(f, np.rad2deg(np.unwrap(np.angle(v_o))), color='g')


# In[40]:


v_o_time=np.fft.fftshift(np.fft.ifft(v_o))
plt.plot(np.abs(v_o_time))


# # Seconed Order filters

# ## The band pass filter from ALL ABOUT ELECTRONICS "Band Pass Filter and Band Stop Filter Explained"  @~ 4:03min

# In[41]:


YouTubeVideo('dmPIydL0lyM', width=500, height=400, start=243)


# In[42]:


reset()
net_1=Net('N1'); net_2=Net('N2'); net_3=Net('N3')

vs=SINEV(amplitude=10@u_V, frequency=10@u_kHz)
vs['p', 'n']+=net_1, gnd

highpassFsection=rc_highpass(C_value=1.5@u_nF, R_value=10@u_kOhm)
highpassFsection.SKiDl(net_1, gnd, net_2, gnd)

lowpassFsection=rc_lowpass(C_value=1.5@u_nF, R_value=1@u_kOhm)
lowpassFsection.SKiDl(net_2, gnd, net_3, gnd)


circ=generate_netlist()
print(circ)


# In[43]:


sim=circ.simulator()


# In[44]:


ac_vals=sim.ac(start_frequency=1@u_Hz, stop_frequency=50@u_MHz, number_of_points=10, variation='dec')


# In[45]:


v_i=ac_vals[node(net_1)]
v_o=ac_vals[node(net_3)]
f=ac_vals.frequency


# In[46]:


plt.semilogx(f, 20*np.log10(np.abs(v_o)))
plt.twinx()
plt.semilogx(f, np.rad2deg(np.unwrap(np.angle(v_o))), color='g')


# In[47]:


v_o_time=np.fft.fftshift(np.fft.ifft(v_o))
plt.plot(np.abs(v_o_time))


# ## Citations:
# [1] ALL ABOUT ELECTRONICS. "RC Low Pass Filter Explained," YouTube, Aug 20, 2017. [Video file]. Available: https://youtu.be/_2L0l-E1Wx0. [Accessed: Nov 30, 2020].
# 
# [2] ALL ABOUT ELECTRONICS. "RC High Pass Filter Explained," YouTube, Aug 23, 2017. [Video file]. Available: https://youtu.be/9Dx0b0ukNAM. [Accessed: Nov 30, 2020].
# 
# [2] ALL ABOUT ELECTRONICS. "Band Pass Filter and Band Stop Filter Explained," YouTube, Sep 2, 2017. [Video file]. Available: https://youtu.be/dmPIydL0lyM. [Accessed: Nov 30, 2020].
# 

# In[ ]:





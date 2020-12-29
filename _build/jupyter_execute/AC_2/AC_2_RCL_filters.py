#!/usr/bin/env python
# coding: utf-8

# IN PROGRESS
# 
# TODO:
# - look into ngspice interanls and realy verfiy that there are no eqivlincys to dc internlas with .ac sims 
# - add sympy to subcirucit classes
# - add ac sim tool; try to do an inhertince from last section
# - add rl complimeint filters
# - create table of low-high pass vs rc & rl
# - add rlc resonce
# - add rlc filters
# - add filter desighn tool and filter cirucit generation tool
# - discues L->C passive conversion

# In[1]:


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
import warnings


# In[2]:


from DC_1_Codes import get_skidl_spice_ref
#from AC_2_Codes import 

#notebook specific loading control statements 
get_ipython().run_line_magic('matplotlib', 'inline')
#tool to log notebook internals
#https://github.com/jrjohansson/version_information
get_ipython().run_line_magic('load_ext', 'version_information')
get_ipython().run_line_magic('version_information', 'skidl, PySpice,lcapy, sympy, numpy, matplotlib, pandas, scipy')


# # Basic Passive Filters

# ## The low pass filter from ALL ABOUT ELECTRONICS "RC Low Pass Filter Explained"  @~ 8:35min

# In[3]:


YouTubeVideo('_2L0l-E1Wx0', width=500, height=400, start=515)


# In[4]:


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
        
        
    


# In[5]:


lowpassF=rc_lowpass(C_value=.1@u_uF, R_value=1@u_kOhm)
lowpassF.draw_me()


# In[6]:


reset()
net_1=Net('N1'); net_2=Net('N2'); 

vs=SINEV(ac_magnitude=10@u_V, dc_value=10@u_V)
vs['p', 'n']+=net_1, gnd
lowpassF.SKiDl(net_1, gnd, net_2, gnd)




circ=generate_netlist()
print(circ)


# In[7]:


class ac_ease():
    """
    TODO:
        - independ current sources can have there AC current meassued via 
        `@I<name>[acreal]` & `@I<name>[acrimag]` not shure if this is usefull
        also trying via the sensivity 
        - to some serioues testing with ngspice directly to verfiy that internal
        parmters are as limited as they apear to be with .ac
        
    """
    def __init__(self, circ_netlist_obj):
        self.circ_netlist_obj=circ_netlist_obj
        self._build_table()
        
        #dic of allowed AC sweep types
        self.allowed_steptypes_map={'linear':'lin', 'decade':'dec', 'octave': 'oct'}

    
    def _build_table(self):
        self.fsweep_DF=pd.DataFrame(columns=['Start_freq', 'Stop_Freq', 'SamplingInc', 'StepType'])
        self.fsweep_DF.at[len(self.fsweep_DF)]=[.1@u_Hz, 120@u_GHz, 10, 'decade']

    def ac_sweep_setup(self, Start_freq, Stop_Freq, SamplingInc, StepType, display=False):
        """
        TODO:
            -add display action
        """
        #check for allwed step types
        assert StepType in self.allowed_steptypes_map.keys(),  f"{StepType} is not allowed"
        #force start to non zero if sweep not linear
        if StepType != 'linear':
            if float(Start_freq)==0:
                warnings.warn('"linear" is only sweep type that can start at 0Hz,\n setting starting frequancy to 1e-1Hz')
                Start_freq=1e-1@u_Hz
                
                
        
        self.fsweep_DF.at[0]=[Start_freq, Stop_Freq, SamplingInc, StepType]
    
    def _make_sim_control(self):
        
        #check the control table struct
        assert (self.fsweep_DF.columns==['Start_freq', 'Stop_Freq', 'SamplingInc', 'StepType']).all(), 'Contorl Table Column structer has been altered'
        
        #will probily change this down the road
        assert len(self.fsweep_DF)==1, 'there should only be one entry in the control table'
        
        #check the sweep type
        self.fsweep_DF['StepType'][0] in self.allowed_steptypes_map.keys(), f"{self.fsweep_DF['StepType'][0]} is not allowed"
        
        #force start to non zero if sweep not linear
        if self.fsweep_DF['StepType'][0] != 'linear':
            if float(self.fsweep_DF['Start_freq'][0])==0:
                warnings.warn('"linear" is only sweep type that can start at 0Hz,\n setting starting frequancy to 1e-1Hz')
                self.fsweep_DF.at[0, 'Start_freq']=1e-1@u_Hz
        
        self.ac_control={
            'start_frequency':self.fsweep_DF.at[0, 'Start_freq'], 
            'stop_frequency':self.fsweep_DF.at[0, 'Stop_Freq'], 
            'number_of_points':self.fsweep_DF.at[0, 'SamplingInc'],
            'variation': self.allowed_steptypes_map[self.fsweep_DF.at[0, 'StepType']]
        }
        
        
        
    
    def do_ac_sim(self):
        self._make_sim_control()
        self.sim=self.circ_netlist_obj.simulator()
        self.ac_vals=self.sim.ac(**self.ac_control)
        
        self.record_ac_nodebranch()

    
    def record_ac_nodebranch(self):
        self.ac_resultsNB_DF=pd.DataFrame(index=self.ac_vals.frequency.as_ndarray())
        self.ac_resultsNB_DF.index.name='freq[Hz]'
        
        #get the node voltages
        for n in self.circ_netlist_obj.node_names:
            if n=='0':
                continue
            self.ac_resultsNB_DF[n+'_[V]']=self.ac_vals[n].as_ndarray()
        
        #get the current from any voltage sourcs
        for cm in self.circ_netlist_obj.element_names:
            if 'V'==cm[0]:
                self.ac_resultsNB_DF[cm+'_[A]']=-self.ac_vals[cm].as_ndarray()
        


# In[8]:


ac_sweep=ac_ease(circ)
ac_sweep.do_ac_sim()
ac_sweep.ac_resultsNB_DF


# In[9]:


data=ac_sweep.ac_resultsNB_DF
data.copy()


# In[10]:


data['power_[W]']=data['N1_[V]']*data['V1_[A]']
data


# In[11]:


data.columns[0]


# In[12]:


def dB_convert(x):
    if ['[W]', '[VAR]', '[VA]'] in x.name:
        return 10*np.log10(x)
    elif ['[V]', '[A]'] in x.name:
        return 20*np.log10(x)


# In[13]:


data.apply(lambda x: 10*np.log10(x) if x.name in ['[W]', '[VAR]', '[VA]'] else 20*np.log10(x), axis=0)


# In[14]:


data.map(dB_convert)


# In[14]:


10*np.log10(data['power_[W]'])


# In[17]:


20*np.log10(data['N1_[V]'])


# In[16]:


repr(type(data))=="<class 'pandas.core.frame.DataFrame'>"


# In[ ]:


#remove later
angle_phase_unwrap= lambda x: np.rad2deg(np.unwrap(np.angle(x)))


# In[ ]:


class ac_analysis:
    
    def __init__(self, ac_sim_raw_DF):
        #write asserts for ac_sim_DF
        assert repr(type(ac_sim_raw_DF))=="<class 'pandas.core.frame.DataFrame'>", '`ac_sim_raw_DF` must be a dataframe'
        self.ac_sim_raw_DF=ac_sim_raw_DF
    
    def make_real_imag(self):
        self.ac_sim_real_DF=self.ac_sim_real_DF.apply(np.real, axis=0)
        
        self.ac_sim_imag_DF=self.ac_sim_imag_DF.apply(np.imag, axis=0)
        
    
    def make_mag_phase(self, deg=True, phase_unwrap=True):
        self.ac_sim_mag_DF=self.ac_sim_raw_DF.apply(np.abs, axis=1)
        if phase_unwrap!=True:
            self.ac_sim_phase_DF=self.ac_sim_raw_DF.apply(np.angle, axis=1, deg=deg)
        else:
            self.ac_sim_phase_DF=self.ac_sim_raw_DF.apply(angle_phase_unwrap, axis=0)
    
    def make_time_domain(self):
        pass
    
    
        


# In[ ]:


class bode_plots:
    pass


# In[ ]:


# .ac will not record current
#ac_sweep.ac_vals['I1']


# In[ ]:


ac_sweep.fsweep_DF.at[0, 'Start_freq']=0
ac_sweep.fsweep_DF


# In[ ]:


sim=circ.simulator()


# In[ ]:


#ac_vals=sim.ac(start_frequency=1e-1@u_Hz, stop_frequency=120@u_GHz, number_of_points=10, variation='dec')
ac_vals=sim.ac(start_frequency=0@u_Hz, stop_frequency=10e3@u_Hz, number_of_points=10000, variation='lin')
#=sim.ac(start_frequency=1e-1@u_Hz, stop_frequency=1@u_THz, number_of_points=10, variation='oct')


# In[ ]:


v_i=ac_vals[node(net_1)]
v_o=ac_vals[node(net_2)]
f=ac_vals.frequency


# In[ ]:


plt.semilogx(f, 10*np.log10(np.abs(v_o)))
plt.twinx()
plt.semilogx(f, np.rad2deg(np.unwrap(np.angle(v_o))), color='g')


# In[ ]:


v_i


# In[ ]:


v_o_time=np.fft.fftshift(np.fft.ifft(v_o))
plt.plot(np.abs(v_o_time))


# In[ ]:


#pz_vals=sim.pole_zero(node(net_1), node(gnd), node(net_2), node(gnd), 'vol', 'pz')


# ## Equivlint Low Pass RL filter 
# 
# The eqivlint RL filter to the above RC filter may be found via the equivlint time consitnc of the RC and RL implemntaiotn such that time constance must mach ie:
# $$RC=\tau_{RC}=\tau_{RL}=L/R$$

# In[ ]:


(1e3)**2 *.1e-9


# In[ ]:


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
        
        
    


# In[ ]:


rl_l=rl_lowpass(L_value=(1e3)**2 *.1e-9 @u_H, R_value=1e3)
rl_l.draw_me()


# In[ ]:


reset()
net_1=Net('N1'); net_2=Net('N2'); 

vs=SINEV(amplitude=10@u_V, frequency=2@u_kHz)
vs['p', 'n']+=net_1, gnd
rl_l.SKiDl(net_1, gnd, net_2, gnd)




circ=generate_netlist()
print(circ)


# In[ ]:


sim=circ.simulator()


# In[ ]:


ac_vals=sim.ac(start_frequency=1@u_Hz, stop_frequency=2@u_MHz, number_of_points=10, variation='dec')


# In[ ]:


v_i=ac_vals[node(net_1)]
v_o=ac_vals[node(net_2)]
f=ac_vals.frequency


# In[ ]:


plt.semilogx(f, 20*np.log10(np.abs(v_o)))
plt.twinx()
plt.semilogx(f, np.rad2deg(np.unwrap(np.angle(v_o))), color='g')


# In[ ]:


v_o_time=np.fft.fftshift(np.fft.ifft(v_o))
plt.plot(np.abs(v_o_time))


# ## The high pass filter from ALL ABOUT ELECTRONICS "RC High Pass Filter Explained"  @~ 7:57min

# In[ ]:


YouTubeVideo('9Dx0b0ukNAM', width=500, height=400, start=477)


# In[ ]:


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
        
        
    


# In[ ]:


highpassF=rc_highpass(C_value=1.5@u_nF, R_value=10@u_kOhm)
highpassF.draw_me()


# In[ ]:


reset()
net_1=Net('N1'); net_2=Net('N2'); 

vs=SINEV(amplitude=10@u_V, frequency=10@u_kHz)
vs['p', 'n']+=net_1, gnd
highpassF.SKiDl(net_1, gnd, net_2, gnd)




circ=generate_netlist()
print(circ)


# In[ ]:


sim=circ.simulator()


# In[ ]:


ac_vals=sim.ac(start_frequency=1@u_Hz, stop_frequency=1@u_MHz, number_of_points=10, variation='dec')


# In[ ]:


v_i=ac_vals[node(net_1)]
v_o=ac_vals[node(net_2)]
f=ac_vals.frequency


# In[ ]:


plt.semilogx(f, 20*np.log10(np.abs(v_o)))
plt.twinx()
plt.semilogx(f, np.rad2deg(np.unwrap(np.angle(v_o))), color='g')


# In[ ]:


v_o_time=np.fft.fftshift(np.fft.ifft(v_o))
plt.plot(np.abs(v_o_time))


# ## RL Highpass

# In[ ]:


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
        
        
    


# In[ ]:


(10e3)**2 *1.5e-9


# In[ ]:


rl_h=rl_highpass(L_value=(10e3)**2 *1.5e-9 @u_H, R_value=10@u_kOhm)
rl_h.draw_me()


# In[ ]:


reset()
net_1=Net('N1'); net_2=Net('N2'); 

vs=SINEV(amplitude=10@u_V, frequency=10@u_kHz)
vs['p', 'n']+=net_1, gnd
rl_h.SKiDl(net_1, gnd, net_2, gnd)




circ=generate_netlist()
print(circ)


# In[ ]:


sim=circ.simulator()


# In[ ]:


ac_vals=sim.ac(start_frequency=1@u_Hz, stop_frequency=1@u_MHz, number_of_points=10, variation='dec')


# In[ ]:


v_i=ac_vals[node(net_1)]
v_o=ac_vals[node(net_2)]
f=ac_vals.frequency


# In[ ]:


plt.semilogx(f, 20*np.log10(np.abs(v_o)))
plt.twinx()
plt.semilogx(f, np.rad2deg(np.unwrap(np.angle(v_o))), color='g')


# In[ ]:


v_o_time=np.fft.fftshift(np.fft.ifft(v_o))
plt.plot(np.abs(v_o_time))


# # Seconed Order filters

# ## The band pass filter from ALL ABOUT ELECTRONICS "Band Pass Filter and Band Stop Filter Explained"  @~ 4:03min

# In[ ]:


YouTubeVideo('dmPIydL0lyM', width=500, height=400, start=243)


# In[ ]:


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


# In[ ]:


sim=circ.simulator()


# In[ ]:


ac_vals=sim.ac(start_frequency=1@u_Hz, stop_frequency=50@u_MHz, number_of_points=10, variation='dec')


# In[ ]:


v_i=ac_vals[node(net_1)]
v_o=ac_vals[node(net_3)]
f=ac_vals.frequency


# In[ ]:


plt.semilogx(f, 20*np.log10(np.abs(v_o)))
plt.twinx()
plt.semilogx(f, np.rad2deg(np.unwrap(np.angle(v_o))), color='g')


# In[ ]:


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





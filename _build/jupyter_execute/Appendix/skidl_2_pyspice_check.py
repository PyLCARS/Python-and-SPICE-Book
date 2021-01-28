#!/usr/bin/env python
# coding: utf-8

# In[1]:


from skidl.pyspice import *
from PySpice.Spice.Netlist import Circuit


# # Checking tool

# In[2]:


def netlist_comp_check(skidl_netlist, pyspice_netlist):
    """
    Simple dumb check tool to compare the netlist from sckidl and pyspice 
    
    Args:
        skidl_netlist (PySpice.Spice.Netlist.Circuit): resulting netlist obj from
            skidl using skidl's `generate_netlist` utlity to compare to pyspice direct
            creation
        
        pyspice_netlist (PySpice.Spice.Netlist.Circuit): circuit obj created directly in pyspice via
            `PySpice.Spice.Netlist.Circuit` to compare it's netlist to skidl produced one
    
    Returns:
        if skidl_netlist is longer then pyspice_netlist will return string statment saying: 'skidl_netlist is longer then pyspice_netlist'
        
        if skidl_netlist is shorter then pyspice_netlist will return string statment saying: 'skidl_netlist is shorter then pyspice_netlist'
        
        if skidl_netlist and pyspice_netlist are equall and but there are diffrances then will print
        message of thoes difrances(|1 indexed) and return a list of indexs where the skidl netlist is differs from the pyspice one
        
        if skidl_netlist == pyspice_netlist then will return the word: 'Match'
    
    TODO: Where should I start
    """
    #only care about the final netlist string
    skidl_netlist=skidl_netlist.str()
    pyspice_netlist=pyspice_netlist.str()
    
    #check the lengths
    if len(skidl_netlist)>len(pyspice_netlist):
        return('skidl_netlist is longer then pyspice_netlist')
    elif len(skidl_netlist)<len(pyspice_netlist):
        return('skidl_netlist is shorter then pyspice_netlist') 
    
    #compare strings char by char
    else:
        string_check=[i for i in range(len(skidl_netlist)) if skidl_netlist[i] != pyspice_netlist[i]]
        if string_check==[]:
            return 'Match'
        else:
            print('Match failed skidl_netlist:')
            print(f'{[i|1 for i in string_check]}')
            return string_check


# # Basic Elements

# ## A            | XSPICE code model (not checked)
# Will check XSPICE elements in seperate (to be done) section in the appendx 
# 
# PySpice/PySpice/Spice/BasicElement.py; (need to find):
# 
# skidl/skidl/libs/pyspice_sklib.py; name="A"

# ## B            | Behavioral (arbitrary) source (not checked)
# PySpice/PySpice/Spice/BasicElement.py; class BehavioralSource:
# 
# skidl/skidl/libs/pyspice_sklib.py; name="B"
# 
# ngspice 5.1: Bxxxx: Nonlinear dependent source (ASRC): BXXXXXXX n| n- <i=expr > <v=expr > <tc1=value > <tc2=value > <temp=value > <dtemp=value >

# ## C            | Capacitor 
# PySpice/PySpice/Spice/BasicElement.py; class Capacitor(DipoleElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="C"
# 
# ngspice 3.2.5 Capacitors: 
# 
# CXXXXXXX n| n- <value > <mname > <m=val> <scale=val> <temp=val> <dtemp=val> <tc1=val> <tc2=val> <ic=init_condition >
#     
# ### Notes

# In[3]:


reset()
net_1=Net('N1'); net_2=Net('N2')
skidl_C=C(ref='1', value=5, scale=5, temp=5, dtemp=5, ic=5, m=5)
skidl_C['p', 'n']+=net_1, net_2
skidl_circ=generate_netlist()
print(skidl_circ)


# In[4]:


pyspice_circ=Circuit('')
pyspice_circ.C('1', 'N1', 'N2', 5, scale=5, temp=5, dtemp=5, ic=5, m=5)
print(pyspice_circ)


# In[5]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## D            | Diode   
# PySpice/PySpice/Spice/BasicElement.py; class Diode(FixedPinElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="D"
# 
# ngspice 7.1 Junction Diodes: 
# 
# DXXXXXXX n| n- mname <area=val> <m=val> <pj=val> <off> <ic=vd> <temp=val> <dtemp=val>
#     
# ### Notes
# - `ic`: did not work in eather skidl or pyspice

# In[6]:


reset()
net_1=Net('N1'); net_2=Net('N2')
skidl_D=D(ref='1',model=5, area=5, m=5, pj=5, off=5, temp=5, dtemp=5)
skidl_D['p', 'n']+=net_1, net_2
skidl_circ=generate_netlist()
print(skidl_circ)


# In[7]:


pyspice_circ=Circuit('')
pyspice_circ.D('1', 'N1', 'N2', model=5, area=5, m=5, pj=5, off=5, temp=5, dtemp=5)
print(pyspice_circ)


# In[8]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## E            | Voltage-controlled voltage source (VCVS)
# PySpice/PySpice/Spice/BasicElement.py; class VoltageControlledVoltageSource(TwoPortElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="E"
# 
# ngspice 4.2.2 Exxxx: Linear Voltage-Controlled Voltage Sources (VCVS): 
# 
# EXXXXXXX N| N- NC| NC- VALUE
#     
# ### Notes

# In[9]:


reset()
net_1=Net('N1'); net_2=Net('N2'); net_3=Net('N3'); net_4=Net('N4')
skidl_E=E(ref='1', voltage_gain=5)
skidl_E['ip', 'in']+=net_1, net_2; skidl_E['op', 'on']+=net_3, net_4
skidl_circ=generate_netlist()
print(skidl_circ)


# In[10]:


pyspice_circ=Circuit('')
pyspice_circ.VoltageControlledVoltageSource('1', 'N3', 'N4', 'N1', 'N2', voltage_gain=5)
print(pyspice_circ)


# In[11]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## F            | Current-controlled current source (CCCs) 
# PySpice/PySpice/Spice/BasicElement.py; class CurrentControlledCurrentSource(DipoleElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="F"
# 
# ngspice 4.2.3 Fxxxx: Linear Current-Controlled Current Sources (CCCS): 
# 
# FXXXXXXX N| N- VNAM VALUE <m=val>
# ### Notes

# In[12]:


reset()
net_1=Net('N1'); net_2=Net('N2')
skidl_F=F(ref='1', control='V1', current_gain=5, m=5)
skidl_F['p', 'n']+=net_1, net_2;
skidl_circ=generate_netlist()
print(skidl_circ)


# In[13]:


pyspice_circ=Circuit('')
pyspice_circ.CurrentControlledCurrentSource('1', 'N1', 'N2', 'V1', current_gain=5, m=5)
print(pyspice_circ)


# In[14]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## G            | Voltage-controlled current source (VCCS)
# PySpice/PySpice/Spice/BasicElement.py; class VoltageControlledCurrentSource(TwoPortElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="G"
# 
# ngspice 4.2.1 Gxxxx: Linear Voltage-Controlled Current Sources (VCCS):
# 
# GXXXXXXX N| N- NC| NC- VALUE <m=val>
# 
# ### Notes
# - 'transconductance' did not work in skidl; but `gain` did as did `current_gain`

# In[15]:


reset()
net_1=Net('N1'); net_2=Net('N2'); net_3=Net('N3'); net_4=Net('N4')
skidl_G=G(ref='1', current_gain=5, m=5)
skidl_G['ip', 'in']+=net_1, net_2; skidl_G['op', 'on']+=net_3, net_4
skidl_circ=generate_netlist()
print(skidl_circ)


# In[16]:


pyspice_circ=Circuit('')
pyspice_circ.VoltageControlledCurrentSource('1', 'N3', 'N4', 'N1', 'N2', transconductance=5, m=5)
print(pyspice_circ)


# In[17]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## H            | Current-controlled voltage source (CCVS)
# PySpice/PySpice/Spice/BasicElement.py; class CurrentControlledVoltageSource(DipoleElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="H"
# 
# ngspice 4.2.4 Hxxxx: Linear Current-Controlled Voltage Sources (CCVS): 
# 
# HXXXXXXX n| n- vnam val
#     
# ### Notes

# In[18]:


reset()
net_1=Net('N1'); net_2=Net('N2')
skidl_H=H(ref='1', control='V1', transresistance=5)
skidl_H['p', 'n']+=net_1, net_2;
skidl_circ=generate_netlist()
print(skidl_circ)


# In[19]:


pyspice_circ=Circuit('')
pyspice_circ.CurrentControlledVoltageSource('1', 'N1', 'N2', 'V1', transresistance=5)
print(pyspice_circ)


# In[20]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## I            | Current source        
# 
# PySpice/PySpice/Spice/BasicElement.py; class CurrentSource(DipoleElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="I"
# 
# ngspice 4.1 Independent Sources for Voltage or Current: 
# 
# IYYYYYYY N| N- <<DC> 
#     
# ### Notes
# - a reduced version of ngspices IYYYYYYY only generating the arguement for <<DC> DC/TRAN VALUE >

# In[21]:


reset()
net_1=Net('N1'); net_2=Net('N2')
skidl_I=I(ref='1', dc_value=5)
skidl_I['p', 'n']+=net_1, net_2
skidl_circ=generate_netlist()
print(skidl_circ)


# In[22]:


pyspice_circ=Circuit('')
pyspice_circ.I('1', 'N1', 'N2', dc_value=5)
print(pyspice_circ)


# In[23]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## J            | Junction field effect transistor (JFET) 
# PySpice/PySpice/Spice/BasicElement.py; class JunctionFieldEffectTransistor(JfetElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="J"
# 
# ngspice 9.1 Junction Field-Effect Transistors (JFETs): 
# 
# JXXXXXXX nd ng ns mname <area > <off> <ic=vds,vgs> <temp=t>
#     
# ### Notes
# - `ic`: did not work in eather skidl or pyspice

# In[24]:


reset()
net_1=Net('N1'); net_2=Net('N2'); net_3=Net('N3')
skidl_J=J(ref='1',model=5, area=5, m=5, off=5, temp=5)
skidl_J['d', 'g', 's']+=net_1, net_2, net_3
skidl_circ=generate_netlist()
print(skidl_circ)


# In[25]:


pyspice_circ=Circuit('')
pyspice_circ.J('1', 'N1', 'N2', 'N3', model=5, area=5, m=5, off=5, temp=5)
print(pyspice_circ)


# In[26]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## K            | Coupled (Mutual) Inductors       
# 
# PySpice/PySpice/Spice/BasicElement.py; class CoupledInductor(AnyPinElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="K"
# 
# ngspice 3.2.11 Coupled (Mutual) Inductors: 
# 
# KXXXXXXX LYYYYYYY LZZZZZZZ value
# 
# ### Notes
# - the inductors must already exsist for pyspice to work
# - K can exspect any value not just in the 0-.99 range; this needs to be fixed down in pyspice

# In[27]:


reset()
net_1=Net('N1'); net_2=Net('N2')
skidl_L1=L(ref='1', value=5, m=5, temp=5, dtemp=5, ic=5)
skidl_L1['p', 'n']+=net_1, net_2

skidl_L2=L(ref='2', value=5, m=5, temp=5, dtemp=5, ic=5)
skidl_L2['p', 'n']+=net_1, net_2

skidl_K=K(ind1=skidl_L1, ind2=skidl_L2, coupling=5)
skidl_circ=generate_netlist()
print(skidl_circ)


# In[28]:


pyspice_circ=Circuit('')
#inductors need to exsist to then be coupled
pyspice_circ.L('1', 'N1', 'N2', 5, m=5, temp=5, dtemp=5, ic=5)
pyspice_circ.L('2', 'N1', 'N2', 5, m=5, temp=5, dtemp=5, ic=5)
pyspice_circ.K('1', 'L1', 'L2', coupling_factor=5)
print(pyspice_circ)


# In[29]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## L            | Inductor        
# 
# PySpice/PySpice/Spice/BasicElement.py; class Inductor(DipoleElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="L"
# 
# ngspice 3.2.9 Inductors: 
# 
# LYYYYYYY n| n- <value > <mname > <nt=val> <m=val> <scale=val> <temp=val> <dtemp=val> <tc1=val> <tc2=val> <ic=init_condition >
#     
# ### Notes

# In[30]:


reset()
net_1=Net('N1'); net_2=Net('N2')
skidl_L=L(ref='1', value=5, m=5, temp=5, dtemp=5, ic=5)
skidl_L['p', 'n']+=net_1, net_2
skidl_circ=generate_netlist()
print(skidl_circ)


# In[31]:


pyspice_circ=Circuit('')
pyspice_circ.L('1', 'N1', 'N2', 5, m=5, temp=5, dtemp=5, ic=5)
print(pyspice_circ)


# In[32]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## M            | Metal oxide field effect transistor (MOSFET)
# PySpice/PySpice/Spice/BasicElement.py; class Mosfet(FixedPinElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="M"
# 
# ngspice 11.1 MOSFET devices: 
# 
# MXXXXXXX nd ng ns nb mname <m=val> <l=val> <w=val> <ad=val> <as=val> <pd=val> <ps=val> <nrd=val> <nrs=val> <off> <ic=vds, vgs, vbs> <temp=t>
#     
# ### Notes
# - `ic`: did not work in eather skidl or pyspice

# In[33]:


reset()
net_1=Net('N1'); net_2=Net('N2'); net_3=Net('N3'); net_4=Net('N4')
skidl_M=M(ref='1', model=5, m=5, l=5, w=5, 
               drain_area=5, source_area=5, drain_perimeter=5, source_perimeter=5, 
               drain_number_square=5, source_number_square=5,
              off=5, temp=5)

skidl_M['d', 'g', 's', 'b']+=net_1, net_2, net_3, net_4
skidl_circ=generate_netlist()
print(skidl_circ)


# In[34]:


pyspice_circ=Circuit('')
pyspice_circ.M('1', 'N1', 'N2', 'N3', 'N4', model=5, m=5, l=5, w=5, 
               drain_area=5, source_area=5, drain_perimeter=5, source_perimeter=5, 
               drain_number_square=5, source_number_square=5,
              off=5, temp=5)
print(pyspice_circ)


# In[35]:


netlist_comp_check(skidl_circ, pyspice_circ)


# | N            | Numerical device for GSS                             |

# | O            | Lossy transmission line                              |

# | P            | Coupled multiconductor line (CPL)                    |

# ## Q            | Bipolar junction transistor (BJT) 
# PySpice/PySpice/Spice/BasicElement.py; class BipolarJunctionTransistor(FixedPinElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="Q"
# 
# ngspice 8.1 Bipolar Junction Transistors (BJTs): 
# 
# QXXXXXXX nc nb ne <ns> mname <area=val> <areac=val> <areab=val> <m=val> <off> <ic=vbe,vce> <temp=val> <dtemp=val>
#     
# ### Notes
# - could not get the substrate connection working in pyspice but it worked fine with skidl
# - `ic`: did not work in eather skidl or pyspice

# In[36]:


reset()
net_1=Net('N1'); net_2=Net('N2'); net_3=Net('N3'); net_4=Net('N4')
skidl_Q=Q(ref='1',model=5, 
          area=5, areab=5, areac=5,
          m=5, off=5, temp=5, dtemp=5)
skidl_Q['c', 'b', 'e']+=net_1, net_2, net_3

#skidl will make the substrate connection fine but could not get pyspice to do so
#therefore skiping for the time being
#skidl_Q['s']+=net_4

skidl_circ=generate_netlist()
print(skidl_circ)


# In[37]:


pyspice_circ=Circuit('')
pyspice_circ.Q('1', 'N1', 'N2', 'N3', model=5, area=5, areab=5, areac=5,
          m=5, off=5, temp=5, dtemp=5, 
               #could not get the substrate connection working in pyspice
              #ns='N4'
              )
               
print(pyspice_circ)


# In[38]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## R            | Resistor         
# 
# PySpice/PySpice/Spice/BasicElement.py; class Resistor(DipoleElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="R"
# 
# ngspice 3.2.1 Resistors: 
# 
# RXXXXXXX n| n- <resistance|r=>value <ac=val> <m=val> <scale=val> <temp=val> <dtemp=val> <tc1=val> <tc2=val> <noisy=0|1>
#     
# ### Notes

# In[39]:


reset()
net_1=Net('N1'); net_2=Net('N2')
skidl_R=R(ref='1', value=5, ac=5, m=5, scale=5, temp=5, dtemp=5, noisy=1)
skidl_R['p', 'n']+=net_1, net_2
skidl_circ=generate_netlist()
print(skidl_circ)


# In[40]:


pyspice_circ=Circuit('')
pyspice_circ.R('1', 'N1', 'N2', 5, ac=5, m=5, scale=5, temp=5, dtemp=5, noisy=1)
print(pyspice_circ)


# In[41]:


netlist_comp_check(skidl_circ, pyspice_circ)


# | S            | Switch (voltage-controlled)                          |

# | T            | Lossless transmission line                           |

# | U            | Uniformly distributed RC line                        |

# ## V | Voltage source      
# 
# PySpice/PySpice/Spice/BasicElement.py; class VoltageSource(DipoleElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="V"
# 
# ngspice 4.1 Independent Sources for Voltage or Current: 
# 
# VXXXXXXX N| N- <<DC> DC/TRAN VALUE >
#     
# ### Notes
# - a reduced version of ngspices VXXXXXXX only generating the arguement for <<DC> DC/TRAN VALUE >

# In[42]:


reset()
net_1=Net('N1'); net_2=Net('N2')
skidl_V=V(ref='1', dc_value=5)
skidl_V['p', 'n']+=net_1, net_2
skidl_circ=generate_netlist()
print(skidl_circ)


# In[43]:


pyspice_circ=Circuit('')
pyspice_circ.V('1', 'N1', 'N2', dc_value=5)
print(pyspice_circ)


# In[44]:


netlist_comp_check(skidl_circ, pyspice_circ)


# | W            | Switch (current-controlled)                          |

# | X            | Subcircuit                                           |

# | Y            | Single lossy transmission line (TXL)                 |

# | Z            | Metal semiconductor field effect transistor (MESFET) |

# ## Z            | Metal semiconductor field effect transistor (MESFET)
# PySpice/PySpice/Spice/BasicElement.py; class Mesfet(JfetElement)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="Z"
# 
# ngspice 10.1 MESFETs: 
# 
# ZXXXXXXX ND NG NS MNAME <AREA > <OFF> <IC=VDS, VGS>
#     
# ### Notes
# - `ic`: did not work in eather skidl or pyspice

# In[45]:


reset()
net_1=Net('N1'); net_2=Net('N2'); net_3=Net('N3')
skidl_Z=Z(ref='1',model=5, area=5, m=5, off=5)
skidl_Z['d', 'g', 's']+=net_1, net_2, net_3
skidl_circ=generate_netlist()
print(skidl_circ)


# In[46]:


pyspice_circ=Circuit('')
pyspice_circ.Z('1', 'N1', 'N2', 'N3', model=5, area=5, m=5, off=5)
print(pyspice_circ)


# In[47]:


netlist_comp_check(skidl_circ, pyspice_circ)


# # Highlevel Elements `SinusoidalMixin` Based
# 
# ## Note in Armour's fort added as_phase 
# 
# SinusoidalMixin is the base translation class for sinusoid wave waveform sources, in other words even thou ngspice compines most sinusoid source as just argument extations to exsisting DC source to create AC souces through pyspice to ngspice these elements must be used
# 
# ## `SinusoidalMixin` args:
# 
# 
# | Name | Parameter      | Default Value | Units |
# |------|----------------|---------------|-------|
# | Vo   | offset         |               | V, A  |
# |------|----------------|---------------|-------|
# | Va   | amplitude      |               | V, A  |
# |------|----------------|---------------|-------|
# | f    | frequency      | 1 / TStop     | Hz    |
# |------|----------------|---------------|-------|
# | Td   | delay          | 0.0           | sec   |
# |------|----------------|---------------|-------|
# | Df   | damping factor | 0.01          | 1/sec |
# |------|----------------|---------------|-------|
# 
# so for a AC SIN voltage sours it's output should be equilint to the following:
# 
# $$V(t) = \begin{cases}
#           V_o & \text{if}\ 0 \leq t < T_d, \\
#           V_o + V_a e^{-D_f(t-T_d)} \sin\left(2\pi f (t-T_d)\right) & \text{if}\ T_d \leq t < T_{stop}.
#         \end{cases}$$

# ## SinusoidalVoltageSource (AC)
#     
# 
# PySpice/PySpice/Spice/HighLevelElement.py; class SinusoidalVoltageSource(VoltageSource, VoltageSourceMixinAbc, SinusoidalMixin)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="SINEV"
# 
# ngspice 4.1 Independent Sources for Voltage or Current & 4.1.2 Sinusoidal: 
# 
# VXXXXXXX N+ N- <<DC> DC/TRAN VALUE > <AC <ACMAG <ACPHASE >>> <DISTOF1 <F1MAG <F1PHASE >>> <DISTOF2 <F2MAG <F2PHASE >>> 
#     
# SIN(VO VA FREQ TD THETA PHASE)
#     
# ### Notes
# - a amalgumation of ngspice's  Independent Sources for Voltage & Sinusoidal statment for transint simulations

# In[48]:


reset()
net_1=Net('N1'); net_2=Net('N2')
skidl_SINV=SINEV(ref='1', 
            #transit sim statments
            offset=5,amplitude=5, frequency=5 , delay=5, damping_factor=5,
            #ac sim statments
            ac_magnitude=5, dc_offset=5)

skidl_SINV['p', 'n']+=net_1, net_2
skidl_circ=generate_netlist()
print(skidl_circ)


# In[49]:


pyspice_circ=Circuit('')
pyspice_circ.SinusoidalVoltageSource('1', 'N1', 'N2', 
                #transit sim statments
                offset=5,amplitude=5, frequency=5 , delay=5, damping_factor=5,
                #ac sim statments
                ac_magnitude=5, dc_offset=5
                                    
                                    )
print(pyspice_circ)


# In[50]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## SinusoidalCurrentSource (AC)
#     
# 
# PySpice/PySpice/Spice/HighLevelElement.py; class class SinusoidalCurrentSource(CurrentSource, CurrentSourceMixinAbc, SinusoidalMixin):
# 
# 
# skidl/skidl/libs/pyspice_sklib.py; name="SINEI"
# 
# ngspice 4.1 Independent Sources for Voltage or Current & 4.1.2 Sinusoidal: 
# 
# IYYYYYYY N+ N- <<DC> DC/TRAN VALUE > <AC <ACMAG <ACPHASE >>> <DISTOF1 <F1MAG <F1PHASE >>> <DISTOF2 <F2MAG <F2PHASE >>> 
#     
# SIN(VO VA FREQ TD THETA PHASE)
#     
# ### Notes
# - a amalgumation of ngspice's  Independent Sources for Voltage & Sinusoidal statment for transint simulations

# In[51]:


reset()
net_1=Net('N1'); net_2=Net('N2')
skidl_SINI=SINEI(ref='1', 
            #transit sim statments
            offset=5,amplitude=5, frequency=5 , delay=5, damping_factor=5,
            #ac sim statments
            ac_magnitude=5, dc_offset=5)

skidl_SINI['p', 'n']+=net_1, net_2
skidl_circ=generate_netlist()
print(skidl_circ)


# In[52]:


pyspice_circ=Circuit('')
pyspice_circ.SinusoidalCurrentSource('1', 'N1', 'N2', 
                #transit sim statments
                offset=5,amplitude=5, frequency=5 , delay=5, damping_factor=5,
                #ac sim statments
                ac_magnitude=5, dc_offset=5
                                    
                                    )
print(pyspice_circ)


# In[53]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## AcLine(SinusoidalVoltageSource)
# PySpice/PySpice/Spice/HighLevelElement.py; class AcLine(SinusoidalVoltageSource)
# 
# skidl/skidl/libs/pyspice_sklib.py; NOT IMPLIMENTED
# 
# ngspice 4.1 Independent Sources for Voltage or Current: 
# 
# VXXXXXXX N+ N- <<DC> DC/TRAN VALUE > <AC <ACMAG <ACPHASE >>> <DISTOF1 <F1MAG <F1PHASE >>> <DISTOF2 <F2MAG <F2PHASE >>> 
#     
# 
# 
# ### Notes
# - it's a pyspice only wraper around pyspices `SinusoidalVoltageSource` that makes a pure for transisint simulation only SIN voltage source with the only arguments being `rms_voltage` and `frequency`
# - pyspice does the rms to amplitute conversion internaly
# - pyspice does not have a `offset` arg
# - pyspice does not have a `delay` arg
# - pyspice does not have a `damping_factor` arg
# - pyspice does not have a `ac_magnitude` arg
# - pyspice does not have a `dc_offset` arg
# - pspice still gives a AC output of the default 1V; this needs to be changed to be equal to `amplitude` internal value or else will give aid in producing incorect results with **ac** simulations

# In[54]:


reset()
net_1=Net('N1'); net_2=Net('N2')
# Skidle does not impliment an AcLine equivlent at this time
skidl_circ=generate_netlist()
print(skidl_circ)


# In[55]:


pyspice_circ=Circuit('')
pyspice_circ.AcLine('1', 'N1', 'N2', 
                #transit sim statments
                rms_voltage=8, frequency=5 
                                    
                                    )
print(pyspice_circ)


# In[56]:


netlist_comp_check(skidl_circ, pyspice_circ)


# # Highlevel Elements `PulseMixin` Based

# # Highlevel Elements `ExponentialMixin` Based
# 
# ExponentialMixin is the base translation class for exponential shped sources used for transisint simulations. Typicly used for simulating responce to charing and discharing events from capcitor/inductor networks. Pyspice does not include ac arguements that are technicly allowed by ngspice
# 
# ## `ExponentialMixin` args:
# 
# 
# 
# | Name | Parameter          | Default Value | Units |
# |------|--------------------|---------------|-------|
# | V1   | Initial value      |               | V, A  |
# |------|--------------------|---------------|-------|
# | V2   | pulsed value       |               | V, A  |
# |------|--------------------|---------------|-------|
# | Td1  | rise delay time    | 0.0           | sec   |
# |------|--------------------|---------------|-------|
# | tau1 | rise time constant | Tstep         | sec   |
# |------|--------------------|---------------|-------|
# | Td2  | fall delay time    | Td1|Tstep     | sec   |
# |------|--------------------|---------------|-------|
# | tau2 | fall time constant | Tstep         | sec   |
# |------|--------------------|---------------|-------|
# 
# so for a expoential based voltage source it's output should be equilint to the following:
# 
# $$V(t) = \begin{cases}
#           V_1 & \text{if}\ 0 \leq t < T_{d1}, \\
#           V_1 + V_{21} ( 1 − e^{-\frac{t-T_{d1}}{\tau_1}} )
#           & \text{if}\ T_{d1} \leq t < T_{d2}, \\
#           V_1 + V_{21} ( 1 − e^{-\frac{t-T_{d1}}{\tau_1}} ) + V_{12} ( 1 − e^{-\frac{t-T_{d2}}{\tau_2}} )
#           & \text{if}\ T_{d2} \leq t < T_{stop}
#         \end{cases}$$
# 
# where $V_{21} = V_2 - V_1$ and $V_{12} = V_1 - V_2$

# ## ExponentialVoltageSource 
#     
# 
# PySpice/PySpice/Spice/HighLevelElement.py; class ExponentialVoltageSource(VoltageSource, VoltageSourceMixinAbc, ExponentialMixin)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="EXPV"
# 
# ngspice 4.1 Independent Sources for Voltage or Current & 4.1.3 Exponential: 
# 
# VXXXXXXX N+ N- 
#     
# EXP(V1 V2 TD1 TAU1 TD2 TAU2)
#     
# ### Notes
# - should technicly also alow dc and ac values from ngspice Independent voltage source statment

# In[57]:


reset()
net_1=Net('N1'); net_2=Net('N2')
skidl_EXPV=EXPV(ref='1', 
            #transit sim statments
            initial_value=5,pulsed_value=5, rise_delay_time=5 , rise_time_constant=5, fall_delay_time=5, fall_time_constant=5,
           )

skidl_EXPV['p', 'n']+=net_1, net_2
skidl_circ=generate_netlist()
print(skidl_circ)


# In[58]:


pyspice_circ=Circuit('')
pyspice_circ.ExponentialVoltageSource('1', 'N1', 'N2', 
            #transit sim statments
            initial_value=5,pulsed_value=5, rise_delay_time=5 , rise_time_constant=5, fall_delay_time=5, fall_time_constant=5,
            
                                    )
print(pyspice_circ)


# In[59]:


netlist_comp_check(skidl_circ, pyspice_circ)


# ## ExponentialCurrentSource 
#     
# 
# PySpice/PySpice/Spice/HighLevelElement.py; class ExponentialCurrentSource(VoltageSource, VoltageSourceMixinAbc, ExponentialMixin)
# 
# skidl/skidl/libs/pyspice_sklib.py; name="EXPI"
# 
# ngspice 4.1 Independent Sources for Voltage or Current & 4.1.3 Exponential: 
# 
# IXXXXXXX N+ N- 
#     
# EXP(I1 I2 TD1 TAU1 TD2 TAU2)
#     
# ### Notes
# - should technicly also alow dc and ac values from ngspice Independent voltage source statment

# In[60]:


reset()
net_1=Net('N1'); net_2=Net('N2')
skidl_EXPI=EXPI(ref='1', 
            #transit sim statments
            initial_value=5,pulsed_value=5, rise_delay_time=5 , rise_time_constant=5, fall_delay_time=5, fall_time_constant=5,
           )

skidl_EXPI['p', 'n']+=net_1, net_2
skidl_circ=generate_netlist()
print(skidl_circ)


# In[61]:


pyspice_circ=Circuit('')
pyspice_circ.ExponentialCurrentSource('1', 'N1', 'N2', 
            #transit sim statments
            initial_value=5,pulsed_value=5, rise_delay_time=5 , rise_time_constant=5, fall_delay_time=5, fall_time_constant=5,
            
                                    )
print(pyspice_circ)


# In[62]:


netlist_comp_check(skidl_circ, pyspice_circ)


# # Highlevel Elements `PieceWiseLinearMixin` Based

# # Highlevel Elements `SingleFrequencyFMMixin` Based

# # Highlevel Elements `AmplitudeModulatedMixin` Based

# # Highlevel Elements `RandomMixin` Based

# In[ ]:





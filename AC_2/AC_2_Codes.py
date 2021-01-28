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
#chapteer 2 section1 single_arm_mod class
#used as a test load for multiphasic power circuits in this section. 

class single_arm_mod():
    """
    holding class for SKiDl Package and lcapy schematic of a 
    inductor motor amatuer model from
    https://electronics.stackexchange.com/questions/234290/using-a-single-phase-induction-motor-equivalent-circuit-in-ltspice
    used as a test load for multiphase power circuits
    """
 
    @package
    def SKiDl_pack(P=Net(), N=Net()):
        """
        SkiDl package to create a single-phase inductive motor armature
        
        Terminals:
            P: Positive terminal
            N: Negative terminal
        
        Returns:
            returns a model based on Handbook of Electric Power Calculations, Third Edition
            By H. Wayne Beaty Section 6: SINGLE-PHASE MOTORS with values sourced
            from a stackexchange quarry (see docstring right under class)
        """
        
        #top loop elements
        l1=L(value=31.6/314@u_H)
        r1=R(value=0.8@u_Ohm)
        l2=L(value=4.34/314@u_H)
        r2=R(value=2.21/.005@u_Ohm)

        #bottom loop elements
        l3=L(value=31.6/314@u_H)
        r3=R(value=0.8@u_Ohm)
        l4=L(value=4.34/314@u_H)
        r4=R(value=2.21/.005@u_Ohm)

        #conect terminal loop
        P & l1['p', 'n']  & l3['p', 'n'] & N

        #connect top loop
        l1['p'] & r1 & l2 & r2 & l1['n']

        #connect bottom loop
        l3['p'] & r3 & l4 & r4 & l3['n']

    
    def draw_me():
        """
        Lcapy schematic of the SKiDl package in this class
        
        Returns:
            returns a model based on Handbook of Electric Power Calculations, Third Edition
            By H. Wayne Beaty Section 6: SINGLE-PHASE MOTORS with values sourced
            from a stackexchange quarry (see docstring right under class)
        """
        
        schematic=kiwi.Circuit()
        
        #termanl loop
        schematic.add('W P 1;right')
        schematic.add('L1 1 3; down')
        schematic.add('W 3 3_2; down')
        schematic.add('L3 3_2 2; down=0.5')
        schematic.add('W N 2;right')
        
        #top loop
        schematic.add('R1 1 4; right')
        schematic.add('L2 4 5; right')
        schematic.add('R2 5 3_1; down')
        schematic.add('W 3 3_1; right')
        
        #bottom loop
        schematic.add('R3 3_2 6; right')
        schematic.add('L4 6 7; right')
        schematic.add('R4 7 2_1; down')
        schematic.add('W 2 2_1; right=2')
        
        
        schematic.draw()
#chapteer 2 section1 angle_phase_unwrap lambda function
#lambda function to find the unwrapped phase from complex values 

angle_phase_unwrap= lambda x: np.rad2deg(np.unwrap(np.angle(x)))
#chapteer 2 section1 power_calcs_ac class
#class to store power calculations and aid functions as staticmethods

class power_calcs_ac():
    """
    Because I am too lazy to create one that works with time
    based data as well this only works with the .ac Fourier like data
    """
    
    @staticmethod
    def instantaneous_power_calc(voltage, current):
        """
        Staticmethod to compute the instantaneous power: 
        $$p[W]=v \cdot i [V \cdot A]$$
        
        Args:
            voltage (np.float or np.complex array): the voltage in terms of 
                Fourier terms whose length must match the `current` argument
            
            current (np.float or np.complex array): the current in terms of 
                Fourier terms whose length must match the `voltage` argument
        Returns:
            returns an array with the so-called instantaneous power in Watts
        
        """
        return voltage*current
                               
    
    @staticmethod
    def averge_power_calc(voltage, current):
        """
        Staticmethod to compute the so-called average power: 
        $$P[W]=\dfrac{\Re( v \cdot i^*)}{2} [V \cdot A]$$
        
        Args:
            voltage (np.float or np.complex array): the voltage in terms of 
                Fourier terms whose length must match the `current` argument
            
            current (np.float or np.complex array): the current in terms of 
                Fourier terms whose length must match the `voltage` argument
        Returns:
            returns an array with the so-called average power in Watts
        
        """
        Power=voltage*np.conj(current)
        Power=np.real(Power)/2
        return Power
    
    @staticmethod
    def reactive_power_calc(voltage, current):
        """
        Staticmethod to compute the so-called reactive power: 
        $$Q[VAR]=\dfrac{\Im( v \cdot i^*)}{2} [V \cdot A]$$
        
        Args:
            voltage (np.float or np.complex array): the voltage in terms of 
                Fourier terms whose length must match the `current` argument
            
            current (np.float or np.complex array): the current in terms of 
                Fourier terms whose length must match the `voltage` argument
        Returns:
            returns an array with the so-called reactive power in Volt-Ampere Reactive [VAR]
        
        """
        Qpower=voltage*np.conj(current)
        Qpower=np.imag(Qpower)/2
        return Qpower
    
    @staticmethod
    def complex_power_calc(voltage, current):
        """
        Staticmethod to compute the so-called complex power: 
        $$S[VA]=\dfrac{ v \cdot i^*}{2}  [V \cdot A]$$
        
        Args:
            voltage (np.float or np.complex array): the voltage in terms of 
                Fourier terms whos length must match the `current` argument
            
            current (np.float or np.complex array): the current in terms of 
                Fourier terms whose length must match the `voltage` argument
        Returns:
            returns an array with the so-called apparent power in Volt-Ampere [VA]
        
        """
        Spower=voltage*np.conj(current)
        Spower=Spower/2
        return Spower
    
    @staticmethod
    def apparent_power_calc(voltage, current):
        """
        Staticmethod to compute the so-called apparent power: 
        $$S[VA]=\left| \dfrac{ v \cdot i^*}{2} \right| [V \cdot A]$$
        
        Args:
            voltage (np.float or np.complex array): the voltage in terms of 
                Fourier terms whose length must match the `current` argument
            
            current (np.float or np.complex array): the current in terms of 
                Fourier terms whose length must match the `voltage` argument
        Returns:
            returns an array with the so-called complex power in Volt-Ampere [VA]
        
        """
        Spower=power_calcs_ac.complex_power_calc(voltage, current)
        spower=np.abs(Spower)
        return spower
    
    @staticmethod
    def power_factor_calc(voltage, current):
        """
        Staticmethod to compute the so-called power factor via: 
        $$pf[]=\dfrac{P}{s} [W/VA]$$
        
        Args:
            voltage (np.float or np.complex array): the voltage in terms of 
                Fourier terms whos length must match the `current` argument
            
            current (np.float or np.complex array): the current in terms of 
                Fourier terms whose length must match the `voltage` argument
        Returns:
            returns an array with the so-called power factor that is unitless 
        
        """
        Ppower=power_calcs_ac.averge_power_calc(voltage, current)
        spower=power_calcs_ac.apparent_power_calc(voltage, current)
        pf=Ppower/spower
        return pf
    
    @staticmethod 
    def power_factor_angle_calc(voltage, current, deg=True):
        """
        Staticmethod to compute the so-called power factor angle via: 
        $$\phi_{pf}[deg/rad]=\arccos\left(\dfrac{P}{s}\right) [W/VA]$$
        
        Args:
            voltage (np.float or np.complex array): the voltage in terms of 
                Fourier terms whose length must match the `current` argument
            
            current (np.float or np.complex array): the current in terms of 
                Fourier terms whose length must match the `voltage` argument
            
            deg (bool; True): When True returns the power factor angle in degrees, else in radians
                
        Returns:
            returns an array with the so-called power factor angle in degrees when `deg` is true, else in radians
        
        """
        pf=power_calcs_ac.power_factor_calc(voltage, current)
        pf_ang=np.arccos(pf)
        
        if deg:
            pf_ang=np.rad2deg(pf_ang)
            
        
        return pf_ang
    
    @staticmethod
    def Cparallel_comp(new_pf_ang, voltage, current, frequancy):
        """
        Statcmethod to calculate the parallel compensating capacitor for inductive loads
        based on the following equation
        
        $$ C[F]=\frac{2 P \cdot \left( \tan{\left(\theta_{i} \right)} - \tan{\left(\theta_{f} \right)}\right)}{\omega |V|^{2}} $$

        
        Args:
            new_pf_ang (float; deg):new power factor angle goal
            
            voltage (np.float or np.complex array): the voltage in terms of 
                Fourier terms whose length must match the `current` argument
            
            current (np.float or np.complex array): the current in terms of 
                Fourier terms whose length must match the `voltage` argument
            
            frequency (np.float): the frequencies matching the current and voltage arrays
                
        Returns:
            returns an array with the needed compensation to reach the target power factor in Farads 
        
        """
        Power=power_calcs_ac.averge_power_calc(voltage, current)
        int_pf_ang=power_calcs_ac.power_factor_angle_calc(voltage, current)
        deltaQ=Power*(np.tan(np.deg2rad(int_pf_ang))-np.tan(np.deg2rad(new_pf_ang)))
        Vmag=np.abs(voltage)
        angfreq=2*np.pi*frequancy
        
        C=(2*deltaQ)/(angfreq*Vmag**2)
        
        return C
    
    @staticmethod
    def Lparallel_comp(new_pf_ang, voltage, current, frequancy):
        
        """
        Statcmethod to calculate the parallel compensating inductor for capacitive loads
        based on the following equation
        
        $$L[H]= \frac{|V|^{2}}{2 P \cdot \omega \left( \tan{\left(\theta_{i} \right)} - \tan{\left(\theta_{f} \right)}\right)} $$
        
        Args:
            new_pf_ang (float; deg):new power factor angle goal
            
            voltage (np.float or np.complex array): the voltage in terms of 
                Fourier terms whos length must match the `current` argument
            
            current (np.float or np.complex array): the current in terms of 
                Fourier terms whose length must match the `voltage` argument
            
            frequency (np.float): the frequencies matching the current and voltage arrays
                
        Returns:
            returns an array with the needed compensation to reach the target power factor in Henrys 
        
        """
        
        Power=power_calcs_ac.averge_power_calc(voltage, current)
        int_pf_ang=power_calcs_ac.power_factor_angle_calc(voltage, current)
        deltaQ=Power*(np.tan(np.deg2rad(int_pf_ang))-np.tan(np.deg2rad(new_pf_ang)))
        Vmag=np.abs(voltage)
        angfreq=2*np.pi*frequancy
        
        L=(Vmag**2)/(2*angfreq*deltaQ)
        
        return L
    
    
    @staticmethod
    def autogen_powercalcs(source_dataframe, base_term):
        """
        Static method to automatically calculate and add the dataframe the 
        instantaneous_power_calc, averge_power_calc, reactive_power_calc, complex_power_calc, apparent_power_calc, 
        power_factor_calc power_factor_angle_calc methods in this class
        
        Args:
            source_dataframe (pandas.dataframe): pandas dataframe containing the AC frequencies 
                in the index, columns in the forms of  `<base_term>_v_[V]` and `<base_term>_i_[A]`
            
            base_term (string): the based string as part of the column names of 
                `<base_term>_v_[V]` and `<base_term>_i_[A]`
                that are the source columns that will be sourced to calculate the other measurements
                
        Returns:
            adds the following:
            source_dataframe[f'{base_term}_p_[W]'], source_dataframe[f'{base_term}_P_[W]'], source_dataframe[f'{base_term}_Q_[VAR]']
            source_dataframe[f'{base_term}_S_[VA]'] ,source_dataframe[f'{base_term}_s_[VA]'], source_dataframe[f'{base_term}_pf_[]'], 
            source_dataframe[f'{base_term}_pfang_[deg]']
            
        
        """
        voltage_col=source_dataframe[f'{base_term}_v_[V]']; current_col=source_dataframe[f'{base_term}_i_[A]']
        
        source_dataframe[f'{base_term}_p_[W]']=power_calcs_ac.instantaneous_power_calc(voltage_col, current_col)
        source_dataframe[f'{base_term}_P_[W]']=power_calcs_ac.averge_power_calc(voltage_col, current_col)
        source_dataframe[f'{base_term}_Q_[VAR]']=power_calcs_ac.reactive_power_calc(voltage_col, current_col)
        source_dataframe[f'{base_term}_S_[VA]']=power_calcs_ac.complex_power_calc(voltage_col, current_col)
        source_dataframe[f'{base_term}_s_[VA]']=power_calcs_ac.apparent_power_calc(voltage_col, current_col)
        source_dataframe[f'{base_term}_pf_[]']=power_calcs_ac.power_factor_calc(voltage_col, current_col)
        source_dataframe[f'{base_term}_pfang_[deg]']=power_calcs_ac.power_factor_angle_calc(voltage_col, current_col)
    
    @staticmethod
    def make_power_plot(source_dataframe, base_term, frequncy_scale='linear'):
        """
        Staticmethod to assist in the rapid plot creation of the three common
        AC powers and the power factor angle in a single 2x2 subplot
        
        Args:
            source_dataframe (pandas.dataframe): pandas dataframe containing the AC frequencies 
                in the index, columns in the forms of 
                `<base_term>_P_[W]`, `<base_term>_Q_[VAR]`, `<base_term>_s_[VA]`, `<base_term>_pfang_[deg]`
            
            base_term (string): the based string as part of the column names of 
                `<base_term>_P_[W]`, `<base_term>_Q_[VAR]`, `<base_term>_s_[VA]`, `<base_term>_pfang_[deg]`
                that are the source columns that create the plot
            
            frequncy_scale (string; 'linear'): the frequency scale that that data was collected against, must
                be `linear` for a linear ac sweep or `decade` for a logarithm decade collection that will make the
                plot be semilogx
        
        Returns:
            creates a 2x2 plot fo the power terms 
            
        TODO:
            -add octave scale
            -make plot generation cleaner between the different frequency scales
            -make the pf_deg work for rad as well
            -make it work with dataframes missing power terms, have it make
                them on the fly
            
        """
        assert frequncy_scale in ['linear', 'decade'], 'frequncy_scale must be linear or decade'
        
        fig, [[Pplot, Qplot], [splot, pfangplot]]=plt.subplots(nrows=2, ncols=2, sharex=True)
        plots=[Pplot, Qplot, splot, pfangplot]
        x=source_dataframe.index
        x_label=source_dataframe.index.name
        
        if frequncy_scale=='linear':
            
            if f'{base_term}_P_[W]' in source_dataframe.columns:
                Pplot.plot(x, source_dataframe[f'{base_term}_P_[W]'])

            if f'{base_term}_Q_[VAR]' in source_dataframe.columns:
                Qplot.plot(x, source_dataframe[f'{base_term}_Q_[VAR]'])

            if f'{base_term}_s_[VA]' in source_dataframe.columns:
                splot.plot(x, source_dataframe[f'{base_term}_s_[VA]'])

            if f'{base_term}_pfang_[deg]' in source_dataframe.columns:
                pfangplot.plot(x, source_dataframe[f'{base_term}_pfang_[deg]'])
        
        elif frequncy_scale in ['decade']:
            if f'{base_term}_P_[W]' in source_dataframe.columns:
                Pplot.semilogx(x, source_dataframe[f'{base_term}_P_[W]'])

            if f'{base_term}_Q_[VAR]' in source_dataframe.columns:
                Qplot.semilogx(x, source_dataframe[f'{base_term}_Q_[VAR]'])

            if f'{base_term}_s_[VA]' in source_dataframe.columns:
                splot.semilogx(x, source_dataframe[f'{base_term}_s_[VA]'])

            if f'{base_term}_pfang_[deg]' in source_dataframe.columns:
                pfangplot.semilogx(x, source_dataframe[f'{base_term}_pfang_[deg]'])
        
        #set the ylabels
        Pplot.set_ylabel('P_[W]')
        Qplot.set_ylabel('Q_[VAR]')
        splot.set_ylabel('s_[VA]')
        pfangplot.set_ylabel('pf_ang_[deg]')
        
        #set the xlabel, ticklabel controls, grid for all the subplots
        for plot in plots:
            plot.set_xlabel(x_label)
            plot.ticklabel_format(useOffset=False, axis='y')
            plot.grid()
        
        
        fig.suptitle(f'"{base_term}" power plots')
        plt.tight_layout()
        
#chapteer 2 section1 two_phase_panel class
#class to create a SPICE model of a North American two-phase power panel
#with ground branch and gfi branch in addition to neutral branch

class two_phase_panel:
    """
    Class to contain the needed Skidl package (method) to invoke the circuit for simulations
    and drawing just the packaged circuit via lcapy.
    
    The circuit in hand is a model of a typical North American 2Phase home panel with resistor separated
    neutral from the ground along with open separation to a parallel gfi path for chassis ground to the earth ground
    """
    
    @subcircuit
    def SKiDl_circ(self, Load_A_posterm, Load_A_neuterm,
                   Load_B_posterm, Load_B_neuterm,
                  RaA=0@u_Ohm, RbB=0@u_Ohm, RnN=0@u_Ohm, 
                  Rearth=0@u_Ohm, Rgfi=1e16@u_Ohm):
        """
        SKiDl subcircuit of a primitive representation of a US 120 2 phase circuit
        panel with resistor separated ground from neutral and open from ground return
        
        Terminals:
            Load_A_posterm: Positive terminal for the A phase load 
            Load_A_neuterm: Neutral terminal for the A phase load 
            Load_B_posterm: Positive terminal for the B phase load 
            Load_B_neuterm: Neutral terminal for the B phase load 
        
        Args:
            RaA (float; 0@u_Ohm; ohms): Resistance between the A source and A terminal
            RbB (float; 0@u_Ohm; ohms): Resistance between the B source and B terminal
            RnN (float; 0@u_Ohm; ohms): Resistance between the N source and N terminal
            Rearth (float; 0@u_Ohm; ohms): Resistance between the earth and neutral
            Rgfi (float; 1e16@u_Ohm; ohms): Resistance between neutral and chassis ground gfi branch
        
        Returns:
            returns a SKiDl subcircuit. Some of the internal elements SkiDl objects
            are also stored in `self.panel_internals` which is a dictionary with 
            names for the objects as keys and the SkiDl object's elements as corresponding values
        """
        
        #ngspice treats net a == net A; so add underbar to save from singular matrix
        net_a=Net('a'); net_A=Net('A_')
        net_b=Net('b'); net_B=Net('B_')
        net_n=Net('n'); net_N=Net('N_')

        #create wire models to deal with voltage source parallel 
        #to inductor singularities 
        phase_a_wire=R(ref='aA', value=RaA); phase_a_wire[1, 2]+=net_a, net_A
        phase_b_wire=R(ref='bB', value=RbB); phase_b_wire[1, 2]+=net_b, net_B
        phase_n_wire=R(ref='nN', value=RnN); phase_n_wire[1, 2]+=net_n, net_N


        #bus neutral to earth
        #create wire (resistor) to isolate neutral and earth ground 
        #can be used to model floating neutral
        earth_n_sep=R(ref='gndbus', value=Rearth)
        #floating neutral monitor
        earth_ammeter=V(ref='amm_earth', dc_value=0@u_V)
        earth_n_sep[1, 2]+=net_n, earth_ammeter['p']
        earth_ammeter['n']+=gnd
        
        #gfi branch
        #isolation from netureal and chaiesse groud at load
        gfi_monter=V(ref='amm_gfi', dc_value=0@u_V)
        #gfi branch currernt monter
        chiasse_n_sep=R(ref='neu_chassie', value=Rgfi);
        net_N & gfi_monter & chiasse_n_sep & gnd
        


        #power behind the pannel
        vs_Aphase=SINEV(ref='A', ac_phase=0,  ac_magnitude=120@u_V)
        vs_Aphase['P', 'N']+=net_a, net_n
        vs_Bphase=SINEV(ref='B', ac_phase=np.deg2rad(-90),  ac_magnitude=120@u_V)
        vs_Bphase['N', 'P']+=net_b, net_n
        
        #supciruct connections
        net_A+=Load_A_posterm
        net_B+=Load_B_posterm
        net_N+=Load_A_neuterm, Load_B_neuterm
        
        #stowe away the elements to read from
        self.panel_internals={'AM_gfi':gfi_monter, 'AM_earth':earth_ammeter, 
                'APhase_Source':vs_Aphase,'BPhase_Source': vs_Bphase}
        
        
        
    
    def draw_me(self):
        """
        method to draw a representation of this subcircuit with 
        load representation with lcapy
        """
        schematic=kiwi.Circuit()
        
        #supplies
        schematic.add('VA a n; down, l=$0^{\circ}$')
        schematic.add('VB n b; down, l=$-90^{\circ}$')
        
        #earth ground leg
        schematic.add('W n n_1; left')
        schematic.add('Rgn n_1, n_2; down, f_>=Idump')
        schematic.add('AMearth n_2 0_1; down')
        schematic.add('W 0_1 0; down=0.2, ground')
        
        #wire legs
        schematic.add('RaA a A_; right=2.5, f_>=Isource')
        schematic.add('RbB b B_; right=2.5, f_>=Isource')
        schematic.add('RnN n N_; right=2.5, f<_=Ireturn')
        
        #loads
        schematic.add('RloadA A_ N_; down')
        schematic.add('RloadB N_ B_; down')
        
        #chaisse gfi branch
        schematic.add('W N_ N1; right')
        schematic.add('Ropen N1 N2; down')
        schematic.add('AMgfi N2 0_2; down, f_>=Igfi')
        schematic.add('W 0_2 0_3; down=0.2, cground')
        
        

        
        schematic.draw()

#chapteer 2 section1 three_phaseY_panel class
#class to create a SPICE model of a three-phase Y source
#with ground branch and gfi branch in addition to neutral branch

class three_phaseY_panel:
    """
    Class to contain the needed Skidl package (method) to invoke the circuit for simulations
    and drawing just the packaged circuit via lcapy.
    
    The circuit in hand is a model of a 3 phase Y (star) source circuit with resistor separated
    neutral from the ground along with an open separation to a parallel gfi path for chassis ground to the earth ground
    """
    def __init__(self):
        pass
    
    @subcircuit
    def SKiDl_circ(self, Load_A_posterm, Load_A_neuterm,
                   Load_B_posterm, Load_B_neuterm,
                   Load_C_posterm, Load_C_neuterm,
                  RaA=0@u_Ohm, RbB=0@u_Ohm, RcC=0@u_Ohm, RnN=0@u_Ohm, 
                  Rearth=0@u_Ohm, Rgfi=1e16@u_Ohm):
        """
        SKiDl subcircuit of a primitive representation of a 3 phase Y (star)  circuit
        panel with resistor separated ground from neutral and open from ground return
        
        Terminals:
            Load_A_posterm: Positive terminal for the A phase load 
            Load_A_neuterm: Neutral terminal for the A phase load 
            Load_B_posterm: Positive terminal for the B phase load 
            Load_B_neuterm: Neutral termanl for the B phase load
            Load_C_posterm: Positive terminal for the C phase load 
            Load_C_neuterm: Neutral terminal for the C phase load 
        
        Args:
            RaA (float; 0@u_Ohm; ohms): Resistance between the A source and A terminal
            RbB (float; 0@u_Ohm; ohms): Resistance between the B source and B terminal
            RcC (float; 0@u_Ohm; ohms): Resistance between the C source and C terminal
            RnN (float; 0@u_Ohm; ohms): Resistance between the N source and N terminal
            Rearth (float; 0@u_Ohm; ohms): Resistance between the earth and neutral
            Rgfi (float; 1e16@u_Ohm; ohms): Resistance between neutral and chassis ground gfi branch
        
        Returns:
            returns a SKiDl subcircuit. Some of the internal elements SkiDl objects
            are also stored in `self.panel_internals` which is a dictionary with 
            names for the objects as keys and the SkiDl object's elements as corresponding values
        """
        
        #ngspice treats net a == net A; so add underbar to save from singular matrix
        net_a=Net('a'); net_A=Net('A_')
        net_b=Net('b'); net_B=Net('B_')
        net_c=Net('c'); net_C=Net('C_')
        net_n=Net('n'); net_N=Net('N_')

        #create wire models to deal with voltage source parallel 
        #to inductor singularities 
        phase_a_wire=R(ref='aA', value=RaA); phase_a_wire[1, 2]+=net_a, net_A
        phase_b_wire=R(ref='bB', value=RbB); phase_b_wire[1, 2]+=net_b, net_B
        phase_c_wire=R(ref='cC', value=RcC); phase_c_wire[1, 2]+=net_c, net_C
        phase_n_wire=R(ref='nN', value=RnN); phase_n_wire[1, 2]+=net_n, net_N


        #bus neutral to earth
        #create wire (resistor) to isolate neutral and earth ground 
        #can be used to model floating neutral
        earth_n_sep=R(ref='gndbus', value=Rearth)
        #floating neutral monitor
        earth_ammeter=V(ref='amm_earth', dc_value=0@u_V)
        earth_n_sep[1, 2]+=net_n, earth_ammeter['p']
        earth_ammeter['n']+=gnd
        
        #gfi branch
        #isolation from netureal and chaiesse groud at load
        gfi_monter=V(ref='amm_gfi', dc_value=0@u_V)
        #gfi branch currernt monter
        chiasse_n_sep=R(ref='neu_chassie', value=Rgfi);
        net_N & gfi_monter & chiasse_n_sep & gnd
        


        #power behind the pannel
        vs_Aphase=SINEV(ref='A', ac_phase=0,  ac_magnitude=120@u_V)
        vs_Aphase['P', 'N']+=net_a, net_n
        vs_Bphase=SINEV(ref='B', ac_phase=np.deg2rad(-120),  ac_magnitude=120@u_V)
        vs_Bphase['P', 'N']+=net_b, net_n
        vs_Cphase=SINEV(ref='C', ac_phase=np.deg2rad(120),  ac_magnitude=120@u_V)
        vs_Cphase['P', 'N']+=net_c, net_n
        
        #supciruct connections
        net_A+=Load_A_posterm
        net_B+=Load_B_posterm
        net_C+=Load_C_posterm
        net_N+=Load_A_neuterm, Load_B_neuterm, Load_C_neuterm
        
        #stowe away the elements to read from
        self.panel_internals={'AM_gfi':gfi_monter, 'AM_earth':earth_ammeter, 
                'APhase_Source':vs_Aphase,'BPhase_Source': vs_Bphase, 'CPhase_Source': vs_Cphase}
        
        
        
    
    def draw_me(self):
        """
        method to draw a representation of this subcircuit with 
        load representation with lcapy
        """
        schematic=kiwi.Circuit()
        
        #supplies
        schematic.add('VA a na; down=3, l=$0^{\circ}$')
        schematic.add('W na nb; right')
        schematic.add('VB b nb; down=2, l=$-120^{\circ}$')
        schematic.add('W nb nc; right')
        schematic.add('VC c nc; down, l=$120^{\circ}$')

        
        #earth ground leg
        schematic.add('W na n_1; left')
        schematic.add('Rgn n_1, n_2; down, f_>=Idump')
        schematic.add('AMearth n_2 0_1; down')
        schematic.add('W 0_1 0; down=0.2, ground')
        
        #wire legs
        schematic.add('RaA a A_; right=2.5, f_>=Isource')
        schematic.add('RbB b B_; right=2.5, f_>=Isource')
        schematic.add('RcC c C_; right=2.5, f_>=Isource')
        schematic.add('RnN nc N_C; right=2.5, f<_=Ireturn')
        
        #loads
        schematic.add('RloadC C_ N_C; down')
        schematic.add('W N_C N_B; right')
        schematic.add('RloadB B_ N_B; down')
        schematic.add('W N_B N_A; right')
        schematic.add('RloadA A_ N_A; down')

        
        #chaisse gfi branch
        schematic.add('W N_A N1; right')
        schematic.add('Ropen N1 N2; down')
        schematic.add('AMgfi N2 0_2; down, f_>=Igfi')
        schematic.add('W 0_2 0_3; down=0.2, cground')
        
        

        
        schematic.draw()
#chapteer 2 section1 three_phaseDelta_panel class
#class to create a SPICE model of a three-phase Delta source


class three_phaseDelta_panel:
    def __init__(self):
        pass
    
    @subcircuit
    def SKiDl_circ(self, Aterm_power, Bterm_power, Cterm_power):
        """
        SKiDl subcircuit of a primitive representation of a US 120 2 phase circuit
        panel with resistor separated ground from neutral and open from ground return
        
        Terminals:
            Aterm_power: terminal for the A phase to load 
            Bterm_power: terminal for the B phase to load 
            Cterm_power: terminal for the C phase to load 
            
        Returns:
            returns a SKiDl subcircuit. Some of the internal elements SkiDl objects
            are also stored in `self.panel_internals` which is a dictionary with 
            names for the objects as keys and the SkiDl object's elements as corresponding values
        """

        net_a=Net('a'); net_b=Net('b'); net_c=Net('c')
        
        #things connected to terminal a
        vs_ABphase=SINEV(ref='AB', ac_phase=0,  ac_magnitude=120@u_V)
        rabwire=R(ref='abwire', value=0@u_Ohm)
        ragnd_open=R(ref='aopen', value=1e16@u_Ohm)
        net_a & vs_ABphase['p', 'n'] & rabwire[1,2] & net_b
        net_a & ragnd_open[1, 2] & gnd
        
        #things connected to terminal b
        vs_BCphase=SINEV(ref='BC', ac_phase=np.deg2rad(-120),  ac_magnitude=120@u_V)
        rbcwire=R(ref='bcwire', value=0@u_Ohm)
        rbgnd_open=R(ref='bopen', value=1e16@u_Ohm)
        net_b & vs_BCphase['p', 'n'] & rbcwire[1,2] & net_c
        net_b & rbgnd_open[1, 2] & gnd
        
        #things connected to terminal c
        vs_CAphase=SINEV(ref='CA', ac_phase=np.deg2rad(120),  ac_magnitude=120@u_V)
        rcawire=R(ref='cawire', value=0@u_Ohm)
        rcgnd_open=R(ref='copen', value=1e16@u_Ohm)
        net_c & vs_CAphase['p', 'n'] & rcawire[1,2] & net_a
        net_c & rcgnd_open[1, 2] & gnd


        
        # make the connections to the rest of the circuit
        net_a+=Aterm_power
        net_b+=Bterm_power
        net_c+=Cterm_power
        
        
        
        
        
        #stowe away the elements to read from
        self.panel_internals={
                'ABPhase_Source':vs_ABphase,'BCPhase_Source': vs_BCphase, 'CAPhase_Source': vs_CAphase}
        
        
        
    
    def draw_me(self):
        """
        method to draw a representation of this subcircuit with 
        load representation with lcapy
        """
        schematic=kiwi.Circuit()
        
        schematic.add('W c1 c; right')
        schematic.add('W c2 c1; right')
        schematic.add('W c3 c2; right')
        schematic.add('W vc c3; down')
        schematic.add('Vca vc va; up, l=$120^{\circ}$')
        schematic.add('Rcawire va a3; up')
        schematic.add('W a3 a2; right')
        schematic.add('W a2 a1; right')
        schematic.add('W a1 a; right')
        
        schematic.add('Ragnd a2 0; down')
        schematic.add('Rcgnd c2 0; up')
        schematic.add('Rbgnd 0 b1; right')
        schematic.add('W 0 0_1; left=0.2, ground')

        schematic.add('Vab a1 vb; down, l=$0^{\circ}$')
        schematic.add('Rabwire vb b1; down')
        
        schematic.add('Vbc b1 vc1; down, l=$-120^{\circ}$')
        schematic.add('Rbcwire vc1 c1; down')
        
        schematic.add('W b1 b; right')

        
        schematic.draw()
#chapteer 2 section 2 rc_lowpass filter class
#class with lcapy and skidl subcircuit to create an RC lowpass filter

class rc_lowpass():
    """
    holding class for SkiDl subcircuit and lcapy schematic of a
    lowpass RC filter primitive
    """
    
    def __init__(self, subcirc_ref=None, C_value=1@u_F, R_value=1@u_Ohm):
        """
        Args:
            subcirc_ref (str): reference to use for the base of the internal elements
            C_value (float; 1@u_F; Farads): the capacitance in farads for the RC capacitive element
            R_value (float; 1@u_Ohm; Ohms): the resistance in ohms for the RC resistive element

        Returns:
            None
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
        Terminals:
        term_0, term_1, term_2, term_3
        
        Terminals are defined via:
        ```
        Left_Termanals - RC_Lowpass - Right_Termanals
                             +----+  
        Postive V_i   term_0-|0  2|-term_2     Postive V_o
        Negtive V_i   term_1-|1  3|-term_3     Negtive V_o
                             +----+
        ```
        
        Args:
            return_internls (bool; False): If True return out the internal Voltage Source,
                and Resistance objects in this package
        Returns:
            Returns elements to circuit RC lowpass filter part element object and if `return_internls`
            is True will return the internal voltage and resistance objects in that order 
        """
        if self.subcirc_ref!=None:
            Cref={'ref':f'C_{self.subcirc_ref}'}
            Rref={'ref':f'R_{self.subcirc_ref}'}
        else:
            Cref={}
            Rref={}
            
        self.c=C(value=self.C_value, **Cref)
        self.r=R(value=self.R_value, **Rref)
        
        self.r[1, 2]+=term_0, self.c['p']
        self.c['p', 'n']+=term_2, term_3
        
        if return_elements:
            return self.c, self.r
    
    def lcapy_self(self, draw_me=True, with_values=True):
        """
        Creates a lcapy schematic of this classes filter that
        can be used for amongst other things: draw a basic schematic
        of this class filter, extract the transfer function, 
        exstract the 2Port Repersntation
        
        Args:
            draw_me (bool): will draw a schematic of this classes filter schematic
            
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract simply variables
            
        Return:
            the lcapys circuit object is stored in `self.schematic` and will
            have abstract sympy variable for the elements if `with_values` is False
            and will draw the schematic of just this classes filter if `draw_me` is True
        
        TODO:
            - get the Vin statement into the schematic
        
        """
        self.with_values=with_values
        
        
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 1 1_1; right=2')
        self.schematic.add('W 0 0_1; right')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add('R 0_1 2_1; right')
        self.schematic.add('C 2_1 1_1; down')
        
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 1_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')
        
        if with_values:
            self.schematic=self.schematic.subs({'R':self.R_value, 'C':self.C_value})
        
        if draw_me:
            self.schematic.draw()
    
    def get_tf(self, with_values=True, ZPK=True):
        """
        will extract the symbolic transfer function for this filter
        
        Args:
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
            ZPK (bool): if True will try to return the TF in zero-pole-gain form
                else will try to return the TF in canonical form 
                where the unity coefficient is the highest power of the denominator
        
        """
        self.lcapy_self(draw_me=False, with_values=with_values)
        self.tf=self.schematic.transfer(0, 1, 2, 3)
        
        if ZPK:
            return self.tf.ZPK()
        else:
            return self.tf.canonical()
        
    def get_twoPort(self, network_rep='Y', with_values=True):
        """
        Gets the 2Port network representation of this filter.
        The 2Port representation can be controlled
        
        Args:
            network_rep (str; 'Y'): control string for what
                representation is used to choose from:
                
                'Z': two-port Z-parameters matrix
                'Y': two-port Y-parameters matrix
                'H': two-port H-parameters matrix
                'G': two-port G-parameters matrix
                'ABCD': two-port A-parameters matrix
                'invABCD': two-port B-parameters matrix
                'S': two-port S-parameters matrix
                'T': two-port T-parameters matrix
                
                
        """
        
        self.lcapy_self(draw_me=False, with_values=with_values)
        
        #create an action dict to get the 2P rep
        rep_actions={
            'Z':lambda x: x.Zparams(0, 1, 2, 1),
            'Y':lambda x: x.Yparams(0, 1, 2, 1), 
            'H':lambda x: x.Hparams(0, 1, 2, 1),
            'G':lambda x: x.Gparams(0, 1, 2, 1),
            'ABCD':lambda x: x.Aparams(0, 1, 2, 1),
            'invABCD':lambda x: x.Bparams(0, 1, 2, 1),
            'S':lambda x: x.Sparams(0, 1, 2, 1),
            'T':lambda x: x.Tparams(0, 1, 2, 1),
        }
        
        assert network_rep in rep_actions.keys(), f'`{network_rep}` is not 2Port rep'
        
        return rep_actions[network_rep](self.schematic)
        
            
        
        
#chapteer 2 section 2 ac_ease class
#class to perform .ac simulations with a bit more grace

class ac_ease():
    """
    Class to perform AC (.ac) SPICE simulation with some grace; 
    currently limited to what pyspice and ngspice support
    
    TODO:
        - independent current sources can have their AC current measured via 
        `@I<name>[acreal]` & `@I<name>[acrimag]` not shure if this is usefull
        also trying via the sensitivity 
        - do some serious testing with ngspice directly to verify that internal
        parameters are as limited as they appear to be with .ac
        
    """
    def __init__(self, circ_netlist_obj):
        """
        Class to perform AC (.ac) SPICE simulation with some grace
        
        Args:
            circ_netlist_obj (pyspice.Spice.Netlist.Circuit): the Netlist circuit produced 
                from SKiDl's `generate_netlist()`
        
        Returns: 
            creates a table to control the ac sweep `self.fsweep_DF`
            this table will still need to be filled out before a simulation can be run with `self.do_ac_sim`
            can be filled out manually or with the helper method `self.ac_sweep_setup`
        """
        self.circ_netlist_obj=circ_netlist_obj
        self._build_table()
        
        #dic of allowed AC sweep types
        self.allowed_steptypes_map={'linear':'lin', 'decade':'dec', 'octave': 'oct'}

    
    def _build_table(self):
        """
        protected method to create `self.fsweep_DF` dataframe that stores the controls for the ac simulation
        
        TODO:
            -when pyspice accepts more things to sweep add them below
        """
        self.fsweep_DF=pd.DataFrame(columns=['Start_freq', 'Stop_Freq', 'SamplingInc', 'StepType'])
        self.fsweep_DF.at[len(self.fsweep_DF)]=[.1@u_Hz, 120@u_GHz, 10, 'decade']

    def ac_sweep_setup(self, Start_freq, Stop_Freq, SamplingInc, StepType, display_table=False):
        """
        Helper method to create the `self.fsweep_DF` to control the ac simulation
        
        Args:
            Start_freq (Hertz): starting frequency in Hertz of the ac simulation, must be less than `Stop_Freq` and can
                only be zero if `StepType='linear'`
            
            Stop_Freq (Hertz): stoping  frequency in Hertz of the ac simulation, must be greater than `Start_freq` 
            
            SamplingInc (int): number of samples per StepType interval
            
            StepType (string): string control for the ac simulation Step type.
                must be 'linear' (self-explanatory), 'decade' (log base 10 sampling interval), 
                or 'octave' (starting frequency times 2**n to create sample space a double of the starting frequency
                so that samples are pulled from 2**(n-1) and 2**(n) times the starting frequency)
            
            display_table (bool; False): when true will display the generated `self.fsweep_DF` below
             this method call in a jupyter notebook like environment
            
        TODO:
            -add display action
        """
        #check for allowed step types
        assert StepType in self.allowed_steptypes_map.keys(),  f"{StepType} is not allowed"
        #force start to non zero if sweep not linear
        if StepType != 'linear':
            if float(Start_freq)==0:
                warnings.warn('"linear" is only sweep type that can start at 0Hz,\n setting starting frequancy to 1e-1Hz')
                Start_freq=1e-1@u_Hz
                
                
        #check that stop frequency is greater than start
        assert Stop_Freq>Start_freq, 'Stop frequency must be greater then starting frequency'
        
        self.fsweep_DF.at[0]=[Start_freq, Stop_Freq, SamplingInc, StepType]
        
        if display_table:
            display(self.fsweep_DF)
    
    def _make_sim_control(self):
        """
        Internal method to extract the row information to the .ac pyspice call arguments
        Will raise a warning if the simulation start frequency is 0Hz for non-linear frequency sampling and 
        then set the start frequency to .1Hz
        
        Args:
            NONE
            
        Returns:
            `self.ac_control` what is feed into the .ac to do the simulation over a frequency
            
        """
        
        #check the control table struct
        assert (self.fsweep_DF.columns==['Start_freq', 'Stop_Freq', 'SamplingInc', 'StepType']).all(), 'Contorl Table Column structer has been altered'
        
        #will probably change this down the road
        assert len(self.fsweep_DF)==1, 'there should only be one entry in the control table'
        
        #check the sweep type
        self.fsweep_DF['StepType'][0] in self.allowed_steptypes_map.keys(), f"{self.fsweep_DF['StepType'][0]} is not allowed"
        
        #check that stop frequency is greater than start
        assert self.fsweep_DF.at[0, 'Stop_Freq']>self.fsweep_DF.at[0, 'Start_freq'], 'Stop freqauncy must be grater then starting freauncy'
        
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
        """
        Does a standard Branch and Node .ac simulation for the single filled out row in `self.fsweep_DF`
        
        Args:
            None
        
        Returns: 
            raw results are stored in `self.ac_vals`, processed results are automatically stored in 
            `self.ac_resultsNB_DF` via `self.record_ac_nodebranch`
        """
        self._make_sim_control()
        self.sim=self.circ_netlist_obj.simulator()
        self.ac_vals=self.sim.ac(**self.ac_control)
        
        self.record_ac_nodebranch()

    
    def record_ac_nodebranch(self):
        """ 
        Helper method to put .ac node branch results into a dataframe where the index is the 
        sweep frequency used in the simulation
        
        Args:
            None
        
        Returns:
            `self.ac_resultsNB_DF` which is a pandas dataframe with the index being the sweep frequency
            and the columns being the node voltages and branch currents from any available voltage sources  
            
        TODO:
            look into getting the current in any current sources
            
        """
        self.ac_resultsNB_DF=pd.DataFrame(index=self.ac_vals.frequency.as_ndarray())
        self.ac_resultsNB_DF.index.name='freq[Hz]'
        
        #get the node voltages
        for n in self.circ_netlist_obj.node_names:
            if n=='0':
                continue
            self.ac_resultsNB_DF[n+'_[V]']=self.ac_vals[n].as_ndarray()
        
        #get the current from any voltage source
        for cm in self.circ_netlist_obj.element_names:
            if 'V'==cm[0]:
                self.ac_resultsNB_DF[cm+'_[A]']=-self.ac_vals[cm].as_ndarray()
                
#chapteer 2 section 2 ac_representation_tool class
#class that converts dataframe of raw ac complex data to veries complex
#repsentations 

class ac_representation_tool:
    """
    Class to take a dataframe with AC simulation complex value data and
    represent it in various ways. raw data should come from `ac_ease.ac_resultsNB_DF`
    
    TODO:
        -get the unit renaming in the columns working with regex
    """
    
    def __init__(self, ac_sim_raw_DF):
        """
        pull in the data
        Args:
            ac_sim_raw_DF (pandas dataframe): pandas dataframe of raw data from AC simulation 
                preferbyly from `ac_ease.ac_resultsNB_DF`, index must be the simulation
                frequency and columns must be the complex data
        
        Returns: 
            None
        
        TODO:
            broaden complex assertin to include np.complex128
        """
        #write asserts for ac_sim_DF
        assert repr(type(ac_sim_raw_DF))=="<class 'pandas.core.frame.DataFrame'>", '`ac_sim_raw_DF` must be a dataframe'
        #check that all columns from raw data are complex
        assert (ac_sim_raw_DF.dtypes==np.complex64).all() or (ac_sim_raw_DF.dtypes==np.complex128).all(), 'Raw data must be complex from AC sim'
        self.ac_sim_raw_DF=ac_sim_raw_DF
    
    def make_real_imag(self):
        """
        Method to create a real and image version of the raw data

        Args: None
        
        Returns:
            real values are stored in `self.ac_sim_real_DF`; and 
            imaginary values are stored in `self.ac_sim_imag_DF`
            
        """
        
        self.ac_sim_real_DF=self.ac_sim_raw_DF.apply(np.real, axis=0)
        
        self.ac_sim_imag_DF=self.ac_sim_raw_DF.apply(np.imag, axis=0)
        
    
    def make_mag_phase(self, mag='dB', char_res=50, deg=True, phase_unwrap=True):
        """
        Method to make generate the various magnitude and phase representation of the complex data
        
        Args:
            mag (string, "dB"): control statement to specify the representation of the generated 
                magnitude data; right now only 'dB' and 'abs' are supported
            
            char_res (float; 50; ohms): the characteristic impedance for magnitude representation calculations; not 
                implemented at the moment
            
            deg (bool; True): bool control statement to represent the phase data in degrees if True; else in radians
            
            phase_unwrap (bool; True): when True and `deg` is True will represent the degrees in phased unwrapped
        
        Returns:
            magnitude data is stored in `self.ac_sim_mag_DF` and phase data is stored in `self.ac_sim_phase_DF`
        
        TODO: 
            -complete all of the magnitude conversions
        """
        #deal with the cacophony of magnitudes
        mag_conversions={
            'dB': lambda x: 10*np.log10(np.abs(x)) if x.name in ['[W]', '[VAR]', '[VA]'] else 20*np.log10(np.abs(x)), 
            'abs': lambda x: np.abs(x)
        }
        
        #check the input
        assert mag in mag_conversions.keys(), f'{mag} is not a known magnitude repsentation'
         
        self.ac_sim_mag_DF=self.ac_sim_raw_DF.apply(mag_conversions[mag], axis=0)
        
        #redo the column name units
        #get down with the regex to really do this
        if mag in ['dB']:
            self.ac_sim_mag_DF.rename(columns={i:i+'[dB]' for i in self.ac_sim_mag_DF.columns}, inplace=True)
        
        #deal with the phase
        
        if (deg==True) and (phase_unwrap==True):
            #phase unwrapped lambda function
            angle_phase_unwrap= lambda x: np.rad2deg(np.unwrap(np.angle(x)))

            self.ac_sim_phase_DF=self.ac_sim_raw_DF.apply(angle_phase_unwrap, axis=0)
        
        else:
            self.ac_sim_phase_DF=self.ac_sim_raw_DF.apply(np.angle, axis=0, deg=deg)

        #realy need that stupid regex working
        
        if deg:
            self.ac_sim_phase_DF.rename(columns={i:i+'[deg]' for i in self.ac_sim_phase_DF.columns}, inplace=True)
        else:
            self.ac_sim_phase_DF.rename(columns={i:i+'[rads]' for i in self.ac_sim_phase_DF.columns}, inplace=True)
#chapteer 2 section 2 eecomplex_plot_templets class
#class that stores templets plots for most common complex rep plots

class eecomplex_plot_templets():
    """
    Class that stores basic/common Electrical Engineering Complex value
    representation plots that may be used stand-alone or as templets in other plots 
    with refinements
    """
    def __init__(self):
        pass
    
    def bode_plot_one_templet(self, freq_data, mag_data, phase_data, ax=None, title=''):
        """
        Templet plot to make a Bode plot with Magnitude and Phase all in one
        graph using a twinx. 
        
        Args:
            freq_data (numpy array or pandas series; Hz): the sampling frequency
            
            mag_data (numpy array or pandas seres; dB): the magnitude data in decibels
            
            phase_data (numpy array or pandas series; deg unwrapped): the phase data in degrees unwrapped
            
            ax (matplotlib axis; None): If left None will create a new plot, else must
                be a matplotlib subplot axis to be added to
            
            title (str; ''): Subplot title string
            
        Returns:
            Returns a bode plot, and if an axis was passed to `ax` will be modified
            with how to plot the magnitude
            
        
        TODO:
            - figure out how to return the `ax_phase` generated internally
            - add x,y scale control
        """
        assert len(freq_data)==len(mag_data)==len(phase_data), 'freq_data, mag_data, phase_data, must all be the same length'
        
        if ax!=None:
            assert repr(type(ax))=="<class 'matplotlib.axes._subplots.AxesSubplot'>", 'ax must be a matplotlib axis'

        ax_mag=ax or plt.gca()


        #fig, ax_mag=plt.subplots()

        ax_phase=ax_mag.twinx()

        ax_mag.semilogx(freq_data, mag_data, label='mag')
        ax_phase.semilogx(freq_data, phase_data, color='green', linestyle='--', label='phase')
        

        ax_mag.set_xlabel('frequancy [Hz]')
        ax_mag.set_ylabel('[dB]')
        ax_phase.set_ylabel('[deg]')
        ax_mag.grid()
        
        #make a single legend 
        handles, labels = [(a + b) for a, b in zip(ax_mag.get_legend_handles_labels(), ax_phase.get_legend_handles_labels())]
        ax_phase.legend(handles, labels)

        if title!='':
            title=' of '+title
        ax_mag.set_title(f'Bode Plot{title}')
        
        
        
    def bode_plot_two_templet(self, freq_data, mag_data, phase_data, 
                              axs=None, title=''):
        """
        Templet plot to make a Bode plot with Magnitude and Phase in two separate subplots
        with shared x-axis
        
        Args:
            freq_data (numpy array or pandas series; Hz): the sampling frequency
            
            mag_data (numpy array or pandas seres; dB): the magnitude data in decibels
            
            phase_data (numpy array or pandas series; deg unwrapped): the phase data in degrees unwrapped
            
            axs (list of matplotlib axis; None): If left None will create a new plot, else must 
                be a list of matplotlib subplots axis to be added to where the first entry
                will be the magnitude axis, and the second will be the phase axis
            
            title (str; ''): Subplot title string
            
        Returns:
            Returns a bode plot, and if an axis was passed to `ax` will be modified
            with how to plot the magnitude
            
        
        TODO:
            - add x,y scale control
        """
        
        assert len(freq_data)==len(mag_data)==len(phase_data), 'freq_data, mag_data, phase_data, must all be the same length'
        
       
        if axs==None:
            fig, [ax_mag, ax_phase]=plt.subplots(nrows=2, sharex=True)
        else:
            assert len(axs)==2, 'there should only be two elements in axs'
            
            for i, ax in enumerate(axs):
                assert repr(type(ax))=="<class 'matplotlib.axes._subplots.AxesSubplot'>", f"element {i} in axs was not a matplotlib axis"
            ax_mag=axs[0]; ax_phase=axs[1]
            ax_mag.get_shared_x_axes().join(ax_mag, ax_phase)
        
        #fore the two axes to share x
        ax_mag.xaxis.set_tick_params(which='both', labelbottom=True)




        ax_mag.semilogx(freq_data, mag_data, label='mag')
        ax_phase.semilogx(freq_data, phase_data, color='green', linestyle='--', label='phase')



        ax_mag.set_ylabel('[dB]')
        ax_phase.set_ylabel('[deg]')
        ax_mag.grid()
        ax_phase.grid()
        
        #style the x-axis for both subplots so it's between the two
        ax_phase.set_xlabel('frequancy [Hz]')
        ax_phase.xaxis.set_label_position('top') 
        ax_phase.xaxis.set_ticks_position('top') 
        ax_phase.tick_params(labelbottom=False,labeltop=True)

        if title!='':
            title=' of '+title
        ax_mag.set_title(f'Bode Plot{title}');
        plt.tight_layout()
    
    
    def nichols_plot_templet(self, mag_data, phase_data, ax=None, title=''):
        
        """
        Templet plot to make a Nichols plot with magnitude in the y-axis and
        phase in the x-axis, with a counter arrow showing the parametric direction
        
        Args:
            
            mag_data (numpy array or pandas seres; dB): the magnitude data in decibels
            
            phase_data (numpy array or pandas series; deg unwrapped): the phase data in degrees unwrapped
            
            ax (matplotlib axis; None): If left None will create a new plot, else must
                be a matplotlib subplot axis to be added to
            
            title (str; ''): Subplot title string
            
        Returns:
            Returns a Nichols plot, and if an axis was passed to `ax` will be modified
            with the Nichols plot
            
        
        TODO:
            - add x,y scale control
        """
        
        assert len(mag_data)==len(phase_data), 'mag_data and phase_data, must all be the same length'

        if ax!=None:
            assert repr(type(ax))=="<class 'matplotlib.axes._subplots.AxesSubplot'>", 'ax must be a matplotlib axis'

        ax=ax or plt.gca()

        ax.plot(phase_data, mag_data)
        line=ax.get_lines()[0]
        eecomplex_plot_templets.add_arrow(line)

        
        #xlim
        xmin=phase_data.min()*1.1; xmax=phase_data.max()*1.1
        
        if -1*xmax<xmin:
            xmin=-1*xmax
        
        if -1*xmin>xmax:
            xmax=-1*xmin
            
        ax.set_xlim(xmin, xmax)
        

        ax.set_xlabel('[deg]'); ax.set_ylabel('[dB]')

        ax.grid()
        #ax.axhline(0, linestyle='--', linewidth=2.0, color='black')
        ax.axvline(0, linestyle='--', linewidth=2.0, color='black')
        if title!='':
            title=' of '+title
        ax.set_title(f'Nichols Plot{title}');
        

    def nyquist_plot_templet(self, real_data, imag_data, ax=None, title=''):
        """
        Templet plot to make a Nyquist plot with imaginary in the y-axis and
        real in the x-axis, with a counter arrow showing the parametric direction
        
        Args:
            
            real_data (numpy array or pandas series): the real data 
            
            imag_data (numpy array or pandas series): the imaginary data
            
            ax (matplotlib axis; None): If left None will create a new plot, else must
                be a matplotlib subplot axis to be added to
            
            title (str; ''): Subplot title string
            
        Returns:
            Returns a Nyquist plot, and if an axis was passed to `ax` will be modified
            with the Nyquist plot
            
        
        TODO:
            - add x,y scale control
        """
        assert len(real_data)==len(imag_data), 'real_data and imag_data, must all be the same length'

        if ax!=None:
            assert repr(type(ax))=="<class 'matplotlib.axes._subplots.AxesSubplot'>", 'ax must be a matplotlib axis'

        ax=ax or plt.gca()

        ax.plot(real_data, imag_data)
        line=ax.get_lines()[0]
        eecomplex_plot_templets.add_arrow(line)

        #xlim
        xmin=real_data.min()*1.1; xmax=real_data.max()*1.1

        if -1*xmax<xmin:
            xmin=-1*xmax
        
        if -1*xmin>xmax:
            xmax=-1*xmin
        ax.set_xlim(xmin, xmax)

        #ylim
        ymin=imag_data.min()*1.1; ymax=imag_data.max()*1.1

        if -1*ymax<ymin:
            ymin=-1*ymax
        
        if -1*ymin>ymax:
            ymax=-1*ymin
        ax.set_ylim(ymin, ymax)

        ax.set_xlabel('Real'); ax.set_ylabel('Imag')

        ax.grid()
        ax.axhline(0, linestyle='--', linewidth=2.0, color='black')
        ax.axvline(0, linestyle='--', linewidth=2.0, color='black')
        if title!='':
            title=' of '+title
        ax.set_title(f'Nyquist Plot{title}');
    
    @staticmethod
    def add_arrow(line, positions=None, num_positions=4, direction='right', size=15, color=None):
        """
        add an arrow to a line axis in the direction of the parametric data.

        line:       Line2D object
        positions:   list or array of index positions to draw an arrow(s) at; if None will draw at least one arrow
        num_positions: int; then number arrows to draw along the length of the line; if 1 will draw at the mean
        direction:  'left' or 'right'
        size: the size of the arrow in font-size points
        color:      if None, line color is taken.

        from: https://stackoverflow.com/questions/34017866/arrow-on-a-line-plot-with-matplotlib
        and also use:https://stackoverflow.com/questions/52042183/matplotlib-get-color-for-subplot

        
        """
        #if color is None:
        #    color = line.get_color()

        xdata = line.get_xdata()
        ydata = line.get_ydata()
        
        
        if (positions is None) and (num_positions==1):
            positions=[]
            positions[0] = xdata.mean()
        elif (positions is None) and (num_positions!=1):
            line_len=len(xdata)
            if num_positions>=line_len:
                num_positions==line_len
            positions=[xdata[int(np.ceil(i*line_len/num_positions))] for i in range(num_positions)]
        else:
            assert all(isinstance(i, int) for i in positions), 'positions must be int index positions'
        
        for pos in positions:
            # find the closest index
            start_ind = np.argmin(np.absolute(xdata - pos))
            if direction == 'right':
                end_ind = start_ind + 1
            else:
                end_ind = start_ind - 1

            line.axes.annotate('',
                xytext=(xdata[start_ind], ydata[start_ind]),
                xy=(xdata[end_ind], ydata[end_ind]),
                arrowprops=dict(arrowstyle="->", color=color),
                size=size
            )
#chapteer 2 section 2 qfilter_explorer class
#class to perform the anylsis of the filtes in Ch2 sec 2

class qfilter_explorer(ac_ease, ac_representation_tool, eecomplex_plot_templets):
    
    def __init__(self, circ, title, start_freq=.1@u_Hz, stop_freq=1@u_GHz):
        
        #do what ac_ease is supposed to do at startup
        #instainte the simulation from the circuit
        ac_ease.__init__(self, circ)
        #setup the simulation parameter with the helper method
        self.ac_sweep_setup(start_freq, stop_freq, 20, 'decade', True)
        #do the simulation
        self.do_ac_sim()
        
        
        #do what ac_representation_tool is supposed to do at startup
        #and pass in the selfs from `ac_ease`'s ac_resultsNB_DF
        ac_representation_tool.__init__(self, self.ac_resultsNB_DF)
        #generate the representations
        self.make_real_imag()
        self.make_mag_phase()
        
        
        eecomplex_plot_templets.__init__(self)


        fig=plt.figure(constrained_layout=True, figsize=(8,8))
        spec=fig.add_gridspec(2,2)

        ax_bode=fig.add_subplot(spec[0, :])
        self.bode_plot_one_templet(self.ac_sim_mag_DF.index, self.ac_sim_mag_DF['Out_[V][dB]'], self.ac_sim_phase_DF['Out_[V][deg]'], 
                          title='Out_[V]', ax=ax_bode)

        ax_nichols=fig.add_subplot(spec[1, 0])
        self.nichols_plot_templet(self.ac_sim_mag_DF['Out_[V][dB]'], self.ac_sim_phase_DF['Out_[V][deg]'], 
                          title='Out_[V]', ax=ax_nichols)

        ax_nyquist=fig.add_subplot(spec[1, 1])
        self.nyquist_plot_templet(self.ac_sim_real_DF['Out_[V]'], self.ac_sim_imag_DF['Out_[V]'], 
                          title='Out_[V]', ax=ax_nyquist)
        
        fig.suptitle(title)
    
    def symbolic_tf(self, filter_obj):
        self.symbolic_data=pd.DataFrame(index=self.ac_resultsNB_DF.index)
        
        f=self.symbolic_data.index.values
        
        #get the tf and get the data
        tf=filter_obj.get_tf()
        symbolic_data=tf.frequency_response(self.symbolic_data.index.values).astype('complex')
        
        self.symbolic_data['Out_sym_[V][dB]']=20*np.log10(np.abs(symbolic_data))
        self.symbolic_data['Out_sym_[V][deg]']=np.angle(symbolic_data, deg=True)

        fig, [ax_mag, ax_ph]=plt.subplots(nrows=2, ncols=1)
        
        self.bode_plot_two_templet(self.ac_sim_mag_DF.index, self.ac_sim_mag_DF['Out_[V][dB]'], self.ac_sim_phase_DF['Out_[V][deg]'], 
                          title='Simulated vs Symbolic Out_[V]', axs=[ax_mag, ax_ph])
        
        #add the symbolic data
        ax_mag.semilogx(f, self.symbolic_data['Out_sym_[V][dB]'], linestyle='-.' , alpha=0.5, label='symbolic')
        ax_mag.legend()
        
        ax_ph.semilogx(f, self.symbolic_data['Out_sym_[V][deg]'], color='orange' , alpha=0.75, label='symbolic')
        ax_ph.legend()    
#chapteer 2 section 2 rl_lowpass filter class
#class with lcapy and skidl subcircuit to create an RL lowpass filter

class rl_lowpass():
    """
    holding class for SkiDl subcircuit and lcapy schematic of a
    lowpass RL filter primitive
    """
    def __init__(self, subcirc_ref=None, L_value=1@u_H, R_value=1@u_Ohm):
        """
        Args:
            subcirc_ref (str): reference to use for the base of the internal elements
            L_value (float; 1@u_H; Henery): the inductance in henrys for the RL inductive element
            R_value (float; 1@u_Ohm; Ohms): the resistance in ohms for the RL resistive element

        Returns:
            None
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
        Terminals:
        term_0, term_1, term_2, term_3
        
        Terminals are defined via:
        ```
        Left_Termanals - RL_Lowpass - Right_Termanals
                             +----+  
        Postive V_i   term_0-|0  2|-term_2     Postive V_o
        Negtive V_i   term_1-|1  3|-term_3     Negtive V_o
                             +----+
        ```
        
        Args:
            return_internls (bool; False): If True return out the internal Voltage Source,
                and Resistance objects in this package
        Returns:
            Returns elements to circuit RL lowpass filter part element object and if `return_internls`
            is True will return the internal voltage and resistance objects in that order 
        """
        
        if self.subcirc_ref!=None:
            Lref={'ref':f'L_{self.subcirc_ref}'}
            Rref={'ref':f'R_{self.subcirc_ref}'}
        else:
            Lref={}
            Rref={}
        
        self.l=L(value=self.L_value, **Lref)
        self.r=R(value=self.R_value, **Rref)
        
        self.l['p', 'n']+=term_0, term_2
        self.r[1, 2]+=self.l['n'], term_1
        
        if return_elements:
            return self.l, self.r
    
    def lcapy_self(self, draw_me=True, with_values=True):
        """
        Creates a lcapy schematic of this classes filter that
        can be used for amongst other things: draw a basic schematic
        of this class filter, extract the transfer function, 
        exstract the 2Port Repersntation
        
        Args:
            draw_me (bool): will draw a schematic of this classes filter schematic
            
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
        Return:
            the lcapys circuit object is stored in `self.schematic` and will
            have abstract sympy variable for the elements if `with_values` is False
            and will draw the schematic of just this classes filter if `draw_me` is True
        
        TODO:
            - get the Vin statement into the schematic
        
        """
        
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 1 1_1; right=2')
        self.schematic.add('W 0 0_1; right')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add(f'L 0_1 2_1; right, l=L{str(self.L_value)}')
        self.schematic.add(f'R 2_1 1_1; down, l=R{str(self.R_value)}')
        
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 1_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')
        
        if with_values:
            self.schematic=self.schematic.subs({'R':self.R_value, 'L':self.L_value})
        
        if draw_me:
            self.schematic.draw()

    
    def get_tf(self, with_values=True, ZPK=True):
        """
        will extract the symbolic transfer function for this filter
        
        Args:
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
            ZPK (bool): if True will try to return the TF in zero-pole-gain form
                else will try to return the TF in canonical form 
                where the unity coefficient is the highest power of the denominator
        
        """
        self.lcapy_self(draw_me=False, with_values=with_values)
        self.tf=self.schematic.transfer(0, 1, 2, 3)
        
        if ZPK:
            return self.tf.ZPK()
        else:
            return self.tf.canonical()
        
    def get_twoPort(self, network_rep='Y', with_values=True):
        """
        Gets the 2Port network representation of this filter.
        The 2Port representation can be controlled
        
        Args:
            network_rep (str; 'Y'): control string for what
                representation is used to choose from:
                
                'Z': two-port Z-parameters matrix
                'Y': two-port Y-parameters matrix
                'H': two-port H-parameters matrix
                'G': two-port G-parameters matrix
                'ABCD': two-port A-parameters matrix
                'invABCD': two-port B-parameters matrix
                'S': two-port S-parameters matrix
                'T': two-port T-parameters matrix
                
                
        """
        
        self.lcapy_self(draw_me=False, with_values=with_values)
        
        #create an action dict to get the 2P rep
        rep_actions={
            'Z':lambda x: x.Zparams(0, 1, 2, 1),
            'Y':lambda x: x.Yparams(0, 1, 2, 1), 
            'H':lambda x: x.Hparams(0, 1, 2, 1),
            'G':lambda x: x.Gparams(0, 1, 2, 1),
            'ABCD':lambda x: x.Aparams(0, 1, 2, 1),
            'invABCD':lambda x: x.Bparams(0, 1, 2, 1),
            'S':lambda x: x.Sparams(0, 1, 2, 1),
            'T':lambda x: x.Tparams(0, 1, 2, 1),
        }
        
        assert network_rep in rep_actions.keys(), f'`{network_rep}` is not 2Port rep'
        
        return rep_actions[network_rep](self.schematic)
        
#chapteer 2 section 2 rc_highpass filter class
#class with lcapy and skidl subcircuit to create an RC highpass filter

class rc_highpass():
    """
    holding class for SkiDl subcircuit and lcapy schematic of a
    highpass RC filter primitive
    """
    
    def __init__(self, subcirc_ref=None, C_value=1@u_F, R_value=1@u_Ohm):
        """
        Args:
            subcirc_ref (str): reference to use for the base of the internal elements
            C_value (float; 1@u_F; Farads): the capacitance in farads for the RC capacitive element
            R_value (float; 1@u_Ohm; Ohms): the resistance in ohms for the RC resistive element

        Returns:
            None
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
        Terminals:
        term_0, term_1, term_2, term_3
        
        Terminals are defined via:
        ```
        Left_Termanals - RC_highpass - Right_Termanals
                             +----+  
        Postive V_i   term_0-|0  2|-term_2     Postive V_o
        Negtive V_i   term_1-|1  3|-term_3     Negtive V_o
                             +----+
        ```
        
        Args:
            return_internls (bool; False): If True return out the internal Voltage Source,
                and Resistance objects in this package
        Returns:
            Returns elements to circuit RC highpass filter part element object and if `return_internls`
            is True will return the internal voltage and resistance objects in that order 
        """
        if self.subcirc_ref!=None:
            Cref={'ref':f'C_{self.subcirc_ref}'}
            Rref={'ref':f'R_{self.subcirc_ref}'}
        else:
            Cref={}
            Rref={}
            
        self.c=C(value=self.C_value, **Cref)
        self.r=R(value=self.R_value, **Rref)
        
        self.c['p', 'n']+=term_0, term_2
        self.r[1, 2]+=self.c['n'], term_1
        
        if return_elements:
            return self.c, self.r
    
    def lcapy_self(self, draw_me=True, with_values=True):
        """
        Creates a lcapy schematic of this classes filter that
        can be used for amongst other things: draw a basic schematic
        of this class filter, extract the transfer function, 
        exstract the 2Port Repersntation
        
        Args:
            draw_me (bool): will draw a schematic of this classes filter schematic
            
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
        Return:
            the lcapys circuit object is stored in `self.schematic` and will
            have abstract sympy variable for the elements if `with_values` is False
            and will draw the schematic of just this classes filter if `draw_me` is True
        
        TODO:
            - get the Vin statement into the schematic
        
        """
        
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 1 1_1; right=2')
        self.schematic.add('W 0 0_1; right')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add('C 0_1 2_1; right')
        self.schematic.add('R 2_1 1_1; down')
        
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 1_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')
        
        if with_values:
            self.schematic=self.schematic.subs({'R':self.R_value, 'C':self.C_value})
        
        if draw_me:
            self.schematic.draw()

    def get_tf(self, with_values=True, ZPK=True):
        """
        will extract the symbolic transfer function for this filter
        
        Args:
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
            ZPK (bool): if True will try to return the TF in zero-pole-gain form
                else will try to return the TF in canonical form 
                where the unity coefficient is the highest power of the denominator
        
        """
        self.lcapy_self(draw_me=False, with_values=with_values)
        self.tf=self.schematic.transfer(0, 1, 2, 3)
        
        if ZPK:
            return self.tf.ZPK()
        else:
            return self.tf.canonical()
        
    def get_twoPort(self, network_rep='Y', with_values=True):
        """
        Gets the 2Port network representation of this filter.
        The 2Port representation can be controlled
        
        Args:
            network_rep (str; 'Y'): control string for what
                representation is used to choose from:
                
                'Z': two-port Z-parameters matrix
                'Y': two-port Y-parameters matrix
                'H': two-port H-parameters matrix
                'G': two-port G-parameters matrix
                'ABCD': two-port A-parameters matrix
                'invABCD': two-port B-parameters matrix
                'S': two-port S-parameters matrix
                'T': two-port T-parameters matrix
                
                
        """
        
        self.lcapy_self(draw_me=False, with_values=with_values)
        
        #create an action dict to get the 2P rep
        rep_actions={
            'Z':lambda x: x.Zparams(0, 1, 2, 1),
            'Y':lambda x: x.Yparams(0, 1, 2, 1), 
            'H':lambda x: x.Hparams(0, 1, 2, 1),
            'G':lambda x: x.Gparams(0, 1, 2, 1),
            'ABCD':lambda x: x.Aparams(0, 1, 2, 1),
            'invABCD':lambda x: x.Bparams(0, 1, 2, 1),
            'S':lambda x: x.Sparams(0, 1, 2, 1),
            'T':lambda x: x.Tparams(0, 1, 2, 1),
        }
        
        assert network_rep in rep_actions.keys(), f'`{network_rep}` is not 2Port rep'
        
        return rep_actions[network_rep](self.schematic)
            
        
#chapteer 2 section 2 rl_highpass filter class
#class with lcapy and skidl subcircuit to create an RL highpass filter

class rl_highpass():
    """
    holding class for SkiDl subcircuit and lcapy schematic of a
    highpass RL filter primitive
    """
    
    def __init__(self, subcirc_ref=None, L_value=1@u_H, R_value=1@u_Ohm):
        """
        Args:
            subcirc_ref (str): reference to use for the base of the internal elements
            L_value (float; 1@u_F; Henerys): the inductance in farads for the RL inductive element
            R_value (float; 1@u_Ohm; Ohms): the resistance in ohms for the RL resistive element

        Returns:
            None
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
        Terminals:
        term_0, term_1, term_2, term_3
        
        Terminals are defined via:
        ```
        Left_Termanals - RL_highpass - Right_Termanals
                             +----+  
        Postive V_i   term_0-|0  2|-term_2     Postive V_o
        Negtive V_i   term_1-|1  3|-term_3     Negtive V_o
                             +----+
        ```
        
        Args:
            return_internls (bool; False): If True return out the internal Voltage Source,
                and Resistance objects in this package
        Returns:
            Returns elements to circuit RL highpass filter part element object and if `return_internls`
            is True will return the internal voltage and resistance objects in that order 
        """
        if self.subcirc_ref!=None:
            Lref={'ref':f'L_{self.subcirc_ref}'}
            Rref={'ref':f'R_{self.subcirc_ref}'}
        else:
            Lref={}
            Rref={}
        
        self.l=L(value=self.L_value, **Lref)
        self.r=R(value=self.R_value, **Rref)
        
        self.r[1, 2]+=term_0, self.l['p']
        self.l['p', 'n']+=term_2, term_3
        
        
        
        
        if return_elements:
            return self.l, self.r
    
    def lcapy_self(self, draw_me=True, with_values=True):
        """
        Creates a lcapy schematic of this classes filter that
        can be used for amongst other things: draw a basic schematic
        of this class filter, extract the transfer function, 
        exstract the 2Port Repersntation
        
        Args:
            draw_me (bool): will draw a schematic of this classes filter schematic
            
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
        Return:
            the lcapys circuit object is stored in `self.schematic` and will
            have abstract sympy variable for the elements if `with_values` is False
            and will draw the schematic of just this classes filter if `draw_me` is True
        
        TODO:
            - get the Vin statement into the schematic
        
        """
        
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 1 1_1; right=2')
        self.schematic.add('W 0 0_1; right')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add(f'R 0_1 2_1; right')
        self.schematic.add(f'L 2_1 1_1; down')
        
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 1_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')

        
        if with_values:
            self.schematic=self.schematic.subs({'R':self.R_value, 'L':self.L_value})
        
        if draw_me:
            self.schematic.draw()
    
    def get_tf(self, with_values=True, ZPK=True):
        """
        will extract the symbolic transfer function for this filter
        
        Args:
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
            ZPK (bool): if True will try to return the TF in zero-pole-gain form
                else will try to return the TF in canonical form 
                where the unity coefficient is the highest power of the denominator
        
        """
        self.lcapy_self(draw_me=False, with_values=with_values)
        self.tf=self.schematic.transfer(0, 1, 2, 3)
        
        if ZPK:
            return self.tf.ZPK()
        else:
            return self.tf.canonical()
        
    def get_twoPort(self, network_rep='Y', with_values=True):
        """
        Gets the 2Port network representation of this filter.
        The 2Port representation can be controlled
        
        Args:
            network_rep (str; 'Y'): control string for what
                representation is used to choose from:
                
                'Z': two-port Z-parameters matrix
                'Y': two-port Y-parameters matrix
                'H': two-port H-parameters matrix
                'G': two-port G-parameters matrix
                'ABCD': two-port A-parameters matrix
                'invABCD': two-port B-parameters matrix
                'S': two-port S-parameters matrix
                'T': two-port T-parameters matrix
                
                
        """
        
        self.lcapy_self(draw_me=False, with_values=with_values)
        
        #create an action dict to get the 2P rep
        rep_actions={
            'Z':lambda x: x.Zparams(0, 1, 2, 1),
            'Y':lambda x: x.Yparams(0, 1, 2, 1), 
            'H':lambda x: x.Hparams(0, 1, 2, 1),
            'G':lambda x: x.Gparams(0, 1, 2, 1),
            'ABCD':lambda x: x.Aparams(0, 1, 2, 1),
            'invABCD':lambda x: x.Bparams(0, 1, 2, 1),
            'S':lambda x: x.Sparams(0, 1, 2, 1),
            'T':lambda x: x.Tparams(0, 1, 2, 1),
        }
        
        assert network_rep in rep_actions.keys(), f'`{network_rep}` is not 2Port rep'
        
        return rep_actions[network_rep](self.schematic)
        
            
#chapteer 2 section 2 rlc_series_lowpass filter class
#class with lcapy and skidl subcircuit to create an RLC series lowpass filter

class rlc_series_lowpass():
    """
    holding class for SkiDl subcircuit and lcapy schematic of a
    lowpass series RLC filter primitive
    """
    
    def __init__(self, subcirc_ref=None, L_value=1@u_H, C_value=1@u_F, R_value=1@u_Ohm):
        """
        Args:
            subcirc_ref (str): reference to use for the base of the internal elements
            L_value (float; 1@u_F; Henerys): the inductance in farads for the RLC inductive element
            C_value (float; 1@u_F; Farads): the capacitance in farads for the RLC capacitive element
            R_value (float; 1@u_Ohm; Ohms): the resistance in ohms for the RLC resistive element

        Returns:
            None
        TODO:
            -add assertions
        """
        #add assertions
        self.subcirc_ref=subcirc_ref
        self.L_value=L_value
        self.C_value=C_value
        self.R_value=R_value
        

    @subcircuit
    def SKiDl(self, term_0, term_1, term_2, term_3, return_elements=False):
        """
        Terminals:
        term_0, term_1, term_2, term_3
        
        Terminals are defined via:
        ```
        Left_Termanals - RLC_s_Lowpass - Right_Termanals
                             +----+  
        Postive V_i   term_0-|0  2|-term_2     Postive V_o
        Negtive V_i   term_1-|1  3|-term_3     Negtive V_o
                             +----+
        ```
        
        Args:
            return_internls (bool; False): If True return out the internal Voltage Source,
                and Resistance objects in this package
        Returns:
            Returns elements to circuit RLC series lowpass filter part element object and if `return_internls`
            is True will return the internal voltage and resistance objects in that order 
        """
        if self.subcirc_ref!=None:
            Lref={'ref':f'L_{self.subcirc_ref}'}
            Cref={'ref':f'C_{self.subcirc_ref}'}
            Rref={'ref':f'R_{self.subcirc_ref}'}
        else:
            Lref={}
            Cref={}
            Rref={}
        
        self.l=L(value=self.L_value, **Lref)
        self.c=C(value=self.C_value, **Cref)
        self.r=R(value=self.R_value, **Rref)
        
        self.r[1, 2]+=term_0, self.l[1]
        self.l[2]+=term_2, self.c['p']
        self.c['n']+=term_1, term_3
        
        
        if return_elements:
            return self.c, self.r, self.l
    
    def lcapy_self(self, draw_me=True, with_values=True):
        """
        Creates a lcapy schematic of this classes filter that
        can be used for amongst other things: draw a basic schematic
        of this class filter, extract the transfer function, 
        exstract the 2Port Repersntation
        
        Args:
            draw_me (bool): will draw a schematic of this classes filter schematic
            
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
        Return:
            the lcapys circuit object is stored in `self.schematic` and will
            have abstract sympy variable for the elements if `with_values` is False
            and will draw the schematic of just this classes filter if `draw_me` is True
        
        TODO:
            - get the Vin statement into the schematic
        
        """
        
        
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 0 0_1; right')
        self.schematic.add('W 1 1_1; right=3')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add('R 0_1 N1; right')
        self.schematic.add('L N1 2_1; right')
        self.schematic.add('C 2_1 1_1; down')
        
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 1_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')

        
        if with_values:
            self.schematic=self.schematic.subs({'R':self.R_value, 'L':self.L_value, 'C':self.C_value})
        
        if draw_me:
            self.schematic.draw()
    
    def get_tf(self, with_values=True, ZPK=True):
        """
        will extract the symbolic transfer function for this filter
        
        Args:
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
            ZPK (bool): if True will try to return the TF in zero-pole-gain form
                else will try to return the TF in canonical form 
                where the unity coefficient is the highest power of the denominator
        
        """
        self.lcapy_self(draw_me=False, with_values=with_values)
        self.tf=self.schematic.transfer(0, 1, 2, 3)
        
        if ZPK:
            return self.tf.ZPK()
        else:
            return self.tf.canonical()
        
    def get_twoPort(self, network_rep='Y', with_values=True):
        """
        Gets the 2Port network representation of this filter.
        The 2Port representation can be controlled
        
        Args:
            network_rep (str; 'Y'): control string for what
                representation is used to choose from:
                
                'Z': two-port Z-parameters matrix
                'Y': two-port Y-parameters matrix
                'H': two-port H-parameters matrix
                'G': two-port G-parameters matrix
                'ABCD': two-port A-parameters matrix
                'invABCD': two-port B-parameters matrix
                'S': two-port S-parameters matrix
                'T': two-port T-parameters matrix
                
                
        """
        
        self.lcapy_self(draw_me=False, with_values=with_values)
        
        #create an action dict to get the 2P rep
        rep_actions={
            'Z':lambda x: x.Zparams(0, 1, 2, 1),
            'Y':lambda x: x.Yparams(0, 1, 2, 1), 
            'H':lambda x: x.Hparams(0, 1, 2, 1),
            'G':lambda x: x.Gparams(0, 1, 2, 1),
            'ABCD':lambda x: x.Aparams(0, 1, 2, 1),
            'invABCD':lambda x: x.Bparams(0, 1, 2, 1),
            'S':lambda x: x.Sparams(0, 1, 2, 1),
            'T':lambda x: x.Tparams(0, 1, 2, 1),
        }
        
        assert network_rep in rep_actions.keys(), f'`{network_rep}` is not 2Port rep'
        
        return rep_actions[network_rep](self.schematic)
#chapteer 2 section 2 rlc_series_highpass filter class
#class with lcapy and skidl subcircuit to create an RLC series highpass filter

class rlc_series_highpass():
    """
    holding class for SkiDl subcircuit and lcapy schematic of a
    highpass series RLC filter primitive
    """
    
    def __init__(self, subcirc_ref=None, L_value=1@u_H, C_value=1@u_F, R_value=1@u_Ohm):
        """
        Args:
            subcirc_ref (str): reference to use for the base of the internal elements
            L_value (float; 1@u_F; Henerys): the inductance in farads for the RLC inductive element
            C_value (float; 1@u_F; Farads): the capacitance in farads for the RLC capacitive element
            R_value (float; 1@u_Ohm; Ohms): the resistance in ohms for the RLC resistive element

        Returns:
            None
        TODO:
            -add assertions
        """
        #add assertions
        self.subcirc_ref=subcirc_ref
        self.L_value=L_value
        self.C_value=C_value
        self.R_value=R_value
        

    @subcircuit
    def SKiDl(self, term_0, term_1, term_2, term_3, return_elements=False):
        """
        Terminals:
        term_0, term_1, term_2, term_3
        
        Terminals are defined via:
        ```
        Left_Termanals - RLC_s_highpass - Right_Termanals
                             +----+  
        Postive V_i   term_0-|0  2|-term_2     Postive V_o
        Negtive V_i   term_1-|1  3|-term_3     Negtive V_o
                             +----+
        ```
        
        Args:
            return_internls (bool; False): If True return out the internal Voltage Source,
                and Resistance objects in this package
        Returns:
            Returns elements to circuit RLC series highpass filter part element object and if `return_internls`
            is True will return the internal voltage and resistance objects in that order 
        """
        if self.subcirc_ref!=None:
            Lref={'ref':f'L_{self.subcirc_ref}'}
            Cref={'ref':f'C_{self.subcirc_ref}'}
            Rref={'ref':f'R_{self.subcirc_ref}'}
        else:
            Lref={}
            Cref={}
            Rref={}
        
        self.l=L(value=self.L_value, **Lref)
        self.c=C(value=self.C_value, **Cref)
        self.r=R(value=self.R_value, **Rref)
        
        self.r[1, 2]+=term_0, self.c['p']
        self.c['n']+=term_2, self.l[1]
        self.l[2]+=term_1, term_3
        
        
        if return_elements:
            return self.c, self.r, self.l
    
    def lcapy_self(self, draw_me=True, with_values=True):
        """
        Creates a lcapy schematic of this classes filter that
        can be used for amongst other things: draw a basic schematic
        of this class filter, extract the transfer function, 
        exstract the 2Port Repersntation
        
        Args:
            draw_me (bool): will draw a schematic of this classes filter schematic
            
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
        Return:
            the lcapys circuit object is stored in `self.schematic` and will
            have abstract sympy variable for the elements if `with_values` is False
            and will draw the schematic of just this classes filter if `draw_me` is True
        
        TODO:
            - get the Vin statement into the schematic
        
        """
        
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 0 0_1; right')
        self.schematic.add('W 1 1_1; right=3')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add(f'R 0_1 N1; right')
        self.schematic.add(f'C N1 2_1; right')
        self.schematic.add(f'L 2_1 1_1; down')
        
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 1_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')
        
        if with_values:
            self.schematic=self.schematic.subs({'R':self.R_value, 'L':self.L_value, 'C':self.C_value})
        
        if draw_me:
            self.schematic.draw()
    
    def get_tf(self, with_values=True, ZPK=True):
        """
        will extract the symbolic transfer function for this filter
        
        Args:
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
            ZPK (bool): if True will try to return the TF in zero-pole-gain form
                else will try to return the TF in canonical form 
                where the unity coefficient is the highest power of the denominator
        
        """
        self.lcapy_self(draw_me=False, with_values=with_values)
        self.tf=self.schematic.transfer(0, 1, 2, 3)
        
        if ZPK:
            return self.tf.ZPK()
        else:
            return self.tf.canonical()
        
    def get_twoPort(self, network_rep='Y', with_values=True):
        """
        Gets the 2Port network representation of this filter.
        The 2Port representation can be controlled
        
        Args:
            network_rep (str; 'Y'): control string for what
                representation is used to choose from:
                
                'Z': two-port Z-parameters matrix
                'Y': two-port Y-parameters matrix
                'H': two-port H-parameters matrix
                'G': two-port G-parameters matrix
                'ABCD': two-port A-parameters matrix
                'invABCD': two-port B-parameters matrix
                'S': two-port S-parameters matrix
                'T': two-port T-parameters matrix
                
                
        """
        
        self.lcapy_self(draw_me=False, with_values=with_values)
        
        #create an action dict to get the 2P rep
        rep_actions={
            'Z':lambda x: x.Zparams(0, 1, 2, 1),
            'Y':lambda x: x.Yparams(0, 1, 2, 1), 
            'H':lambda x: x.Hparams(0, 1, 2, 1),
            'G':lambda x: x.Gparams(0, 1, 2, 1),
            'ABCD':lambda x: x.Aparams(0, 1, 2, 1),
            'invABCD':lambda x: x.Bparams(0, 1, 2, 1),
            'S':lambda x: x.Sparams(0, 1, 2, 1),
            'T':lambda x: x.Tparams(0, 1, 2, 1),
        }
        
        assert network_rep in rep_actions.keys(), f'`{network_rep}` is not 2Port rep'
        
        return rep_actions[network_rep](self.schematic)
#chapteer 2 section 2 rlc_series_bandpass filter class
#class with lcapy and skidl subcircuit to create an RLC series bandpass filter

class rlc_series_bandpass():
    """
    holding class for SkiDl subcircuit and lcapy schematic of a
    bandpass series RLC filter primitive
    """
    
    def __init__(self, subcirc_ref=None, L_value=1@u_H, C_value=1@u_F, R_value=1@u_Ohm):
        """
        Args:
            subcirc_ref (str): reference to use for the base of the internal elements
            L_value (float; 1@u_F; Henerys): the inductance in farads for the RLC inductive element
            C_value (float; 1@u_F; Farads): the capacitance in farads for the RLC capacitive element
            R_value (float; 1@u_Ohm; Ohms): the resistance in ohms for the RLC resistive element

        Returns:
            None
        TODO:
            -add assertions
        """
        #add assertions
        self.subcirc_ref=subcirc_ref
        self.L_value=L_value
        self.C_value=C_value
        self.R_value=R_value
        

    @subcircuit
    def SKiDl(self, term_0, term_1, term_2, term_3, return_elements=False):
        """
        Terminals:
        term_0, term_1, term_2, term_3
        
        Terminals are defined via:
        ```
        Left_Termanals - RLC_s_bandpass - Right_Termanals
                             +----+  
        Postive V_i   term_0-|0  2|-term_2     Postive V_o
        Negtive V_i   term_1-|1  3|-term_3     Negtive V_o
                             +----+
        ```
        
        Args:
            return_internls (bool; False): If True return out the internal Voltage Source,
                and Resistance objects in this package
        Returns:
            Returns elements to circuit RLC series bandpass filter part element object and if `return_internls`
            is True will return the internal voltage and resistance objects in that order 
        """
        if self.subcirc_ref!=None:
            Lref={'ref':f'L_{self.subcirc_ref}'}
            Cref={'ref':f'C_{self.subcirc_ref}'}
            Rref={'ref':f'R_{self.subcirc_ref}'}
        else:
            Lref={}
            Cref={}
            Rref={}
        
        self.l=L(value=self.L_value, **Lref)
        self.c=C(value=self.C_value, **Cref)
        self.r=R(value=self.R_value, **Rref)
        
        self.l[1, 2]+=term_0, self.c['p']
        self.c['n']+=term_2, self.r[1]
        self.r[2]+=term_1, term_3
        
        
        if return_elements:
            return self.c, self.r, self.l
    
    def lcapy_self(self, draw_me=True, with_values=True):
        """
        Creates a lcapy schematic of this classes filter that
        can be used for amongst other things: draw a basic schematic
        of this class filter, extract the transfer function, 
        exstract the 2Port Repersntation
        
        Args:
            draw_me (bool): will draw a schematic of this classes filter schematic
            
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
        Return:
            the lcapys circuit object is stored in `self.schematic` and will
            have abstract sympy variable for the elements if `with_values` is False
            and will draw the schematic of just this classes filter if `draw_me` is True
        
        TODO:
            - get the Vin statement into the schematic
        
        """
        
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 0 0_1; right')
        self.schematic.add('W 1 1_1; right=3')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add(f'L 0_1 N1; right')
        self.schematic.add(f'C N1 2_1; right')
        self.schematic.add(f'R 2_1 1_1; down')
        
       
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 1_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')
        
        if with_values:
            self.schematic=self.schematic.subs({'R':self.R_value, 'L':self.L_value, 'C':self.C_value})
        
        if draw_me:
            self.schematic.draw()
    
    def get_tf(self, with_values=True, ZPK=True):
        """
        will extract the symbolic transfer function for this filter
        
        Args:
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
            ZPK (bool): if True will try to return the TF in zero-pole-gain form
                else will try to return the TF in canonical form 
                where the unity coefficient is the highest power of the denominator
        
        """
        self.lcapy_self(draw_me=False, with_values=with_values)
        self.tf=self.schematic.transfer(0, 1, 2, 3)
        
        if ZPK:
            return self.tf.ZPK()
        else:
            return self.tf.canonical()
        
    def get_twoPort(self, network_rep='Y', with_values=True):
        """
        Gets the 2Port network representation of this filter.
        The 2Port representation can be controlled
        
        Args:
            network_rep (str; 'Y'): control string for what
                representation is used to choose from:
                
                'Z': two-port Z-parameters matrix
                'Y': two-port Y-parameters matrix
                'H': two-port H-parameters matrix
                'G': two-port G-parameters matrix
                'ABCD': two-port A-parameters matrix
                'invABCD': two-port B-parameters matrix
                'S': two-port S-parameters matrix
                'T': two-port T-parameters matrix
                
                
        """
        
        self.lcapy_self(draw_me=False, with_values=with_values)
        
        #create an action dict to get the 2P rep
        rep_actions={
            'Z':lambda x: x.Zparams(0, 1, 2, 1),
            'Y':lambda x: x.Yparams(0, 1, 2, 1), 
            'H':lambda x: x.Hparams(0, 1, 2, 1),
            'G':lambda x: x.Gparams(0, 1, 2, 1),
            'ABCD':lambda x: x.Aparams(0, 1, 2, 1),
            'invABCD':lambda x: x.Bparams(0, 1, 2, 1),
            'S':lambda x: x.Sparams(0, 1, 2, 1),
            'T':lambda x: x.Tparams(0, 1, 2, 1),
        }
        
        assert network_rep in rep_actions.keys(), f'`{network_rep}` is not 2Port rep'
        
        return rep_actions[network_rep](self.schematic)
#chapteer 2 section 2 rlc_series_bandstop filter class
#class with lcapy and skidl subcircuit to create an RLC series bandstop filter

class rlc_series_bandstop():
    """
    holding class for SkiDl subcircuit and lcapy schematic of a
    bandstop series RLC filter primitive
    """
    
    def __init__(self, subcirc_ref=None, L_value=1@u_H, C_value=1@u_F, R_value=1@u_Ohm):
        """
        Args:
            subcirc_ref (str): reference to use for the base of the internal elements
            L_value (float; 1@u_F; Henerys): the inductance in farads for the RLC inductive element
            C_value (float; 1@u_F; Farads): the capacitance in farads for the RLC capacitive element
            R_value (float; 1@u_Ohm; Ohms): the resistance in ohms for the RLC resistive element

        Returns:
            None
        TODO:
            -add assertions
        """
        #add assertions
        self.subcirc_ref=subcirc_ref
        self.L_value=L_value
        self.C_value=C_value
        self.R_value=R_value
        

    @subcircuit
    def SKiDl(self, term_0, term_1, term_2, term_3, return_elements=False):
        """
        Terminals:
        term_0, term_1, term_2, term_3
        
        Terminals are defined via:
        ```
        Left_Termanals - RLC_s_bandstop - Right_Termanals
                             +----+  
        Postive V_i   term_0-|0  2|-term_2     Postive V_o
        Negtive V_i   term_1-|1  3|-term_3     Negtive V_o
                             +----+
        ```
        
        Args:
            return_internls (bool; False): If True return out the internal Voltage Source,
                and Resistance objects in this package
        Returns:
            Returns elements to circuit RLC series bandstop filter part element object and if `return_internls`
            is True will return the internal voltage and resistance objects in that order 
        """
        if self.subcirc_ref!=None:
            Lref={'ref':f'L_{self.subcirc_ref}'}
            Cref={'ref':f'C_{self.subcirc_ref}'}
            Rref={'ref':f'R_{self.subcirc_ref}'}
        else:
            Lref={}
            Cref={}
            Rref={}
        
        self.l=L(value=self.L_value, **Lref)
        self.c=C(value=self.C_value, **Cref)
        self.r=R(value=self.R_value, **Rref)
        
        self.r[1, 2]+=term_0, term_2
        self.l[1, 2]+=self.r[2], self.c['p']
        self.c['n']+=term_1, term_3
        
        
        if return_elements:
            return self.c, self.r, self.l
    
    def lcapy_self(self, draw_me=True, with_values=True):
        """
        Creates a lcapy schematic of this classes filter that
        can be used for amongst other things: draw a basic schematic
        of this class filter, extract the transfer function, 
        exstract the 2Port Repersntation
        
        Args:
            draw_me (bool): will draw a schematic of this classes filter schematic
            
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
        Return:
            the lcapys circuit object is stored in `self.schematic` and will
            have abstract sympy variable for the elements if `with_values` is False
            and will draw the schematic of just this classes filter if `draw_me` is True
        
        TODO:
            - get the Vin statement into the schematic
        
        """
        
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 0 0_1; right')
        self.schematic.add('W 1 1_1; right=2')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add(f'R 0_1 2_1; right')
        self.schematic.add(f'L 2_1 N1; down')
        self.schematic.add(f'C N1 1_1; down')
        
       
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 1_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')
        
        if with_values:
            self.schematic=self.schematic.subs({'R':self.R_value, 'L':self.L_value, 'C':self.C_value})
        
        if draw_me:
            self.schematic.draw()
    
    def get_tf(self, with_values=True, ZPK=True):
        """
        will extract the symbolic transfer function for this filter
        
        Args:
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
            ZPK (bool): if True will try to return the TF in zero-pole-gain form
                else will try to return the TF in canonical form 
                where the unity coefficient is the highest power of the denominator
        
        """
        self.lcapy_self(draw_me=False, with_values=with_values)
        self.tf=self.schematic.transfer(0, 1, 2, 3)
        
        if ZPK:
            return self.tf.ZPK()
        else:
            return self.tf.canonical()
        
    def get_twoPort(self, network_rep='Y', with_values=True):
        """
        Gets the 2Port network representation of this filter.
        The 2Port representation can be controlled
        
        Args:
            network_rep (str; 'Y'): control string for what
                representation is used to choose from:
                
                'Z': two-port Z-parameters matrix
                'Y': two-port Y-parameters matrix
                'H': two-port H-parameters matrix
                'G': two-port G-parameters matrix
                'ABCD': two-port A-parameters matrix
                'invABCD': two-port B-parameters matrix
                'S': two-port S-parameters matrix
                'T': two-port T-parameters matrix
                
                
        """
        
        self.lcapy_self(draw_me=False, with_values=with_values)
        
        #create an action dict to get the 2P rep
        rep_actions={
            'Z':lambda x: x.Zparams(0, 1, 2, 1),
            'Y':lambda x: x.Yparams(0, 1, 2, 1), 
            'H':lambda x: x.Hparams(0, 1, 2, 1),
            'G':lambda x: x.Gparams(0, 1, 2, 1),
            'ABCD':lambda x: x.Aparams(0, 1, 2, 1),
            'invABCD':lambda x: x.Bparams(0, 1, 2, 1),
            'S':lambda x: x.Sparams(0, 1, 2, 1),
            'T':lambda x: x.Tparams(0, 1, 2, 1),
        }
        
        assert network_rep in rep_actions.keys(), f'`{network_rep}` is not 2Port rep'
        
        return rep_actions[network_rep](self.schematic)
#chapteer 2 section 2 rlc_parallel_lowpass filter class
#class with lcapy and skidl subcircuit to create an RLC parallel lowpass filter

class rlc_parallel_lowpass():
    """
    holding class for SkiDl subcircuit and lcapy schematic of a
    lowpass parallel RLC filter primitive
    """
    
    def __init__(self, subcirc_ref=None, L_value=1@u_H, C_value=1@u_F, R_value=1@u_Ohm):
        """
        Args:
            subcirc_ref (str): reference to use for the base of the internal elements
            L_value (float; 1@u_F; Henerys): the inductance in farads for the RLC inductive element
            C_value (float; 1@u_F; Farads): the capacitance in farads for the RLC capacitive element
            R_value (float; 1@u_Ohm; Ohms): the resistance in ohms for the RLC resistive element

        Returns:
            None
        TODO:
            -add assertions
        """
        #add assertions
        self.subcirc_ref=subcirc_ref
        self.L_value=L_value
        self.C_value=C_value
        self.R_value=R_value
        

    @subcircuit
    def SKiDl(self, term_0, term_1, term_2, term_3, return_elements=False):
        """
        Terminals:
        term_0, term_1, term_2, term_3
        
        Terminals are defined via:
        ```
        Left_Termanals - RLC_s_Lowpass - Right_Termanals
                             +----+  
        Postive V_i   term_0-|0  2|-term_2     Postive V_o
        Negtive V_i   term_1-|1  3|-term_3     Negtive V_o
                             +----+
        ```
        
        Args:
            return_internls (bool; False): If True return out the internal Voltage Source,
                and Resistance objects in this package
        Returns:
            Returns elements to circuit RLC parallel lowpass filter part element object and if `return_internls`
            is True will return the internal voltage and resistance objects in that order 
        """
        if self.subcirc_ref!=None:
            Lref={'ref':f'L_{self.subcirc_ref}'}
            Cref={'ref':f'C_{self.subcirc_ref}'}
            Rref={'ref':f'R_{self.subcirc_ref}'}
        else:
            Lref={}
            Cref={}
            Rref={}
        
        self.l=L(value=self.L_value, **Lref)
        self.c=C(value=self.C_value, **Cref)
        self.r=R(value=self.R_value, **Rref)
        
        self.l[1, 2]+=term_0, term_2
        self.c['p', 'n']+=term_2, term_1
        self.r[1, 2]+=term_2, term_3
        self.c['n']+=self.r[2]
        
        
        
        if return_elements:
            return self.c, self.r, self.l
    
    def lcapy_self(self, draw_me=True, with_values=True):
        """
        Creates a lcapy schematic of this classes filter that
        can be used for amongst other things: draw a basic schematic
        of this class filter, extract the transfer function, 
        exstract the 2Port Repersntation
        
        Args:
            draw_me (bool): will draw a schematic of this classes filter schematic
            
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
        Return:
            the lcapys circuit object is stored in `self.schematic` and will
            have abstract sympy variable for the elements if `with_values` is False
            and will draw the schematic of just this classes filter if `draw_me` is True
        
        TODO:
            - get the Vin statement into the schematic
        
        """
        
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 0 0_1; right')
        self.schematic.add('W 1 1_1; right=2')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add(f'L 0_1 2_2; right')
        self.schematic.add(f'C 2_2 1_1; down')
        self.schematic.add('W 2_2 2_1; right')
        self.schematic.add('W 1_1 3_1; right')
        self.schematic.add(f'R 2_1 3_1; down')

       
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 3_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')
        
        if with_values:
            self.schematic=self.schematic.subs({'R':self.R_value, 'L':self.L_value, 'C':self.C_value})
        
        if draw_me:
            self.schematic.draw()
    
    def get_tf(self, with_values=True, ZPK=True):
        """
        will extract the symbolic transfer function for this filter
        
        Args:
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
            ZPK (bool): if True will try to return the TF in zero-pole-gain form
                else will try to return the TF in canonical form 
                where the unity coefficient is the highest power of the denominator
        
        """
        self.lcapy_self(draw_me=False, with_values=with_values)
        self.tf=self.schematic.transfer(0, 1, 2, 3)
        
        if ZPK:
            return self.tf.ZPK()
        else:
            return self.tf.canonical()
        
    def get_twoPort(self, network_rep='Y', with_values=True):
        """
        Gets the 2Port network representation of this filter.
        The 2Port representation can be controlled
        
        Args:
            network_rep (str; 'Y'): control string for what
                representation is used to choose from:
                
                'Z': two-port Z-parameters matrix
                'Y': two-port Y-parameters matrix
                'H': two-port H-parameters matrix
                'G': two-port G-parameters matrix
                'ABCD': two-port A-parameters matrix
                'invABCD': two-port B-parameters matrix
                'S': two-port S-parameters matrix
                'T': two-port T-parameters matrix
                
                
        """
        
        self.lcapy_self(draw_me=False, with_values=with_values)
        
        #create an action dict to get the 2P rep
        rep_actions={
            'Z':lambda x: x.Zparams(0, 1, 2, 1),
            'Y':lambda x: x.Yparams(0, 1, 2, 1), 
            'H':lambda x: x.Hparams(0, 1, 2, 1),
            'G':lambda x: x.Gparams(0, 1, 2, 1),
            'ABCD':lambda x: x.Aparams(0, 1, 2, 1),
            'invABCD':lambda x: x.Bparams(0, 1, 2, 1),
            'S':lambda x: x.Sparams(0, 1, 2, 1),
            'T':lambda x: x.Tparams(0, 1, 2, 1),
        }
        
        assert network_rep in rep_actions.keys(), f'`{network_rep}` is not 2Port rep'
        
        return rep_actions[network_rep](self.schematic)
#chapteer 2 section 2 rlc_parallel_highpass filter class
#class with lcapy and skidl subcircuit to create an RLC parallel highpass filter

class rlc_parallel_highpass():
    """
    holding class for SkiDl subcircuit and lcapy schematic of a
    highpass parallel RLC filter primitive
    """
    
    def __init__(self, subcirc_ref=None, L_value=1@u_H, C_value=1@u_F, R_value=1@u_Ohm):
        """
        Args:
            subcirc_ref (str): reference to use for the base of the internal elements
            L_value (float; 1@u_F; Henerys): the inductance in farads for the RLC inductive element
            C_value (float; 1@u_F; Farads): the capacitance in farads for the RLC capacitive element
            R_value (float; 1@u_Ohm; Ohms): the resistance in ohms for the RLC resistive element

        Returns:
            None
        TODO:
            -add assertions
        """
        #add assertions
        self.subcirc_ref=subcirc_ref
        self.L_value=L_value
        self.C_value=C_value
        self.R_value=R_value
        

    @subcircuit
    def SKiDl(self, term_0, term_1, term_2, term_3, return_elements=False):
        """
        Terminals:
        term_0, term_1, term_2, term_3
        
        Terminals are defined via:
        ```
        Left_Termanals - RLC_s_highpass - Right_Termanals
                             +----+  
        Postive V_i   term_0-|0  2|-term_2     Postive V_o
        Negtive V_i   term_1-|1  3|-term_3     Negtive V_o
                             +----+
        ```
        
        Args:
            return_internls (bool; False): If True return out the internal Voltage Source,
                and Resistance objects in this package
        Returns:
            Returns elements to circuit RLC parallel highpass filter part element object and if `return_internls`
            is True will return the internal voltage and resistance objects in that order 
        """
        if self.subcirc_ref!=None:
            Lref={'ref':f'L_{self.subcirc_ref}'}
            Cref={'ref':f'C_{self.subcirc_ref}'}
            Rref={'ref':f'R_{self.subcirc_ref}'}
        else:
            Lref={}
            Cref={}
            Rref={}
        
        self.l=L(value=self.L_value, **Lref)
        self.c=C(value=self.C_value, **Cref)
        self.r=R(value=self.R_value, **Rref)
        
        self.c['p', 'n']+=term_0, term_2
        self.r[1, 2]+=term_2, term_1
        self.l[1, 2]+=term_2, term_3
        self.l[2]+=self.r[2]
        
        
        
        if return_elements:
            return self.c, self.r, self.l
    
    def lcapy_self(self, draw_me=True, with_values=True):
        """
        Creates a lcapy schematic of this classes filter that
        can be used for amongst other things: draw a basic schematic
        of this class filter, extract the transfer function, 
        exstract the 2Port Repersntation
        
        Args:
            draw_me (bool): will draw a schematic of this classes filter schematic
            
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
        Return:
            the lcapys circuit object is stored in `self.schematic` and will
            have abstract sympy variable for the elements if `with_values` is False
            and will draw the schematic of just this classes filter if `draw_me` is True
        
        TODO:
            - get the Vin statement into the schematic
        
        """
        
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 0 0_1; right')
        self.schematic.add('W 1 1_1; right=2')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add(f'C 0_1 2_2; right')
        self.schematic.add(f'R 2_2 1_1; down')
        self.schematic.add('W 2_2 2_1; right')
        self.schematic.add('W 1_1 3_1; right')
        self.schematic.add(f'L 2_1 3_1; down')
        
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 3_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')
        
        if with_values:
            self.schematic=self.schematic.subs({'R':self.R_value, 'L':self.L_value, 'C':self.C_value})
        
        if draw_me:
            self.schematic.draw()
    
    def get_tf(self, with_values=True, ZPK=True):
        """
        will extract the symbolic transfer function for this filter
        
        Args:
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
            ZPK (bool): if True will try to return the TF in zero-pole-gain form
                else will try to return the TF in canonical form 
                where the unity coefficient is the highest power of the denominator
        
        """
        self.lcapy_self(draw_me=False, with_values=with_values)
        self.tf=self.schematic.transfer(0, 1, 2, 3)
        
        if ZPK:
            return self.tf.ZPK()
        else:
            return self.tf.canonical()
        
    def get_twoPort(self, network_rep='Y', with_values=True):
        """
        Gets the 2Port network representation of this filter.
        The 2Port representation can be controlled
        
        Args:
            network_rep (str; 'Y'): control string for what
                representation is used to choose from:
                
                'Z': two-port Z-parameters matrix
                'Y': two-port Y-parameters matrix
                'H': two-port H-parameters matrix
                'G': two-port G-parameters matrix
                'ABCD': two-port A-parameters matrix
                'invABCD': two-port B-parameters matrix
                'S': two-port S-parameters matrix
                'T': two-port T-parameters matrix
                
                
        """
        
        self.lcapy_self(draw_me=False, with_values=with_values)
        
        #create an action dict to get the 2P rep
        rep_actions={
            'Z':lambda x: x.Zparams(0, 1, 2, 1),
            'Y':lambda x: x.Yparams(0, 1, 2, 1), 
            'H':lambda x: x.Hparams(0, 1, 2, 1),
            'G':lambda x: x.Gparams(0, 1, 2, 1),
            'ABCD':lambda x: x.Aparams(0, 1, 2, 1),
            'invABCD':lambda x: x.Bparams(0, 1, 2, 1),
            'S':lambda x: x.Sparams(0, 1, 2, 1),
            'T':lambda x: x.Tparams(0, 1, 2, 1),
        }
        
        assert network_rep in rep_actions.keys(), f'`{network_rep}` is not 2Port rep'
        
        return rep_actions[network_rep](self.schematic)
#chapteer 2 section 2 rlc_parallel_bandpass filter class
#class with lcapy and skidl subcircuit to create an RLC parallel bandpass filter

class rlc_parallel_bandpass():
    """
    holding class for SkiDl subcircuit and lcapy schematic of a
    bandpass parallel RLC filter primitive
    """
    
    def __init__(self, subcirc_ref=None, L_value=1@u_H, C_value=1@u_F, R_value=1@u_Ohm):
        """
        Args:
            subcirc_ref (str): reference to use for the base of the internal elements
            L_value (float; 1@u_F; Henerys): the inductance in farads for the RLC inductive element
            C_value (float; 1@u_F; Farads): the capacitance in farads for the RLC capacitive element
            R_value (float; 1@u_Ohm; Ohms): the resistance in ohms for the RLC resistive element

        Returns:
            None
        TODO:
            -add assertions
        """
        #add assertions
        self.subcirc_ref=subcirc_ref
        self.L_value=L_value
        self.C_value=C_value
        self.R_value=R_value
        

    @subcircuit
    def SKiDl(self, term_0, term_1, term_2, term_3, return_elements=False):
        """
        Terminals:
        term_0, term_1, term_2, term_3
        
        Terminals are defined via:
        ```
        Left_Termanals - RLC_s_bandpass - Right_Termanals
                             +----+  
        Postive V_i   term_0-|0  2|-term_2     Postive V_o
        Negtive V_i   term_1-|1  3|-term_3     Negtive V_o
                             +----+
        ```
        
        Args:
            return_internls (bool; False): If True return out the internal Voltage Source,
                and Resistance objects in this package
        Returns:
            Returns elements to circuit RLC parallel bandpass filter part element object and if `return_internls`
            is True will return the internal voltage and resistance objects in that order 
        """
        if self.subcirc_ref!=None:
            Lref={'ref':f'L_{self.subcirc_ref}'}
            Cref={'ref':f'C_{self.subcirc_ref}'}
            Rref={'ref':f'R_{self.subcirc_ref}'}
        else:
            Lref={}
            Cref={}
            Rref={}
        
        self.l=L(value=self.L_value, **Lref)
        self.c=C(value=self.C_value, **Cref)
        self.r=R(value=self.R_value, **Rref)
        
        self.r[1, 2]+=term_0, term_2
        self.c['p', 'n']+=term_2, term_1
        self.l[1, 2]+=term_2, term_3
        self.l[2]+=self.c['n']
        
        
        
        if return_elements:
            return self.c, self.r, self.l
    
    def lcapy_self(self, draw_me=True, with_values=True):
        """
        Creates a lcapy schematic of this classes filter that
        can be used for amongst other things: draw a basic schematic
        of this class filter, extract the transfer function, 
        exstract the 2Port Repersntation
        
        Args:
            draw_me (bool): will draw a schematic of this classes filter schematic
            
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
        Return:
            the lcapys circuit object is stored in `self.schematic` and will
            have abstract sympy variable for the elements if `with_values` is False
            and will draw the schematic of just this classes filter if `draw_me` is True
        
        TODO:
            - get the Vin statement into the schematic
        
        """
        
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 0 0_1; right')
        self.schematic.add('W 1 1_1; right=2')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add(f'R 0_1 2_2; right')
        self.schematic.add(f'C 2_2 1_1; down')
        self.schematic.add('W 2_2 2_1; right')
        self.schematic.add('W 1_1 3_1; right')
        self.schematic.add(f'L 2_1 3_1; down')
        
        self.schematic.add('W 2_1 2; right=1.5')
        self.schematic.add('W 3_1 3; right=1.5')
        self.schematic.add('P2 2 3; down, v=V_o')
        
        if with_values:
            self.schematic=self.schematic.subs({'R':self.R_value, 'L':self.L_value, 'C':self.C_value})
        
        if draw_me:
            self.schematic.draw()
    
    def get_tf(self, with_values=True, ZPK=True):
        """
        will extract the symbolic transfer function for this filter
        
        Args:
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
            ZPK (bool): if True will try to return the TF in zero-pole-gain form
                else will try to return the TF in canonical form 
                where the unity coefficient is the highest power of the denominator
        
        """
        self.lcapy_self(draw_me=False, with_values=with_values)
        self.tf=self.schematic.transfer(0, 1, 2, 3)
        
        if ZPK:
            return self.tf.ZPK()
        else:
            return self.tf.canonical()
        
    def get_twoPort(self, network_rep='Y', with_values=True):
        """
        Gets the 2Port network representation of this filter.
        The 2Port representation can be controlled
        
        Args:
            network_rep (str; 'Y'): control string for what
                representation is used to choose from:
                
                'Z': two-port Z-parameters matrix
                'Y': two-port Y-parameters matrix
                'H': two-port H-parameters matrix
                'G': two-port G-parameters matrix
                'ABCD': two-port A-parameters matrix
                'invABCD': two-port B-parameters matrix
                'S': two-port S-parameters matrix
                'T': two-port T-parameters matrix
                
                
        """
        
        self.lcapy_self(draw_me=False, with_values=with_values)
        
        #create an action dict to get the 2P rep
        rep_actions={
            'Z':lambda x: x.Zparams(0, 1, 2, 1),
            'Y':lambda x: x.Yparams(0, 1, 2, 1), 
            'H':lambda x: x.Hparams(0, 1, 2, 1),
            'G':lambda x: x.Gparams(0, 1, 2, 1),
            'ABCD':lambda x: x.Aparams(0, 1, 2, 1),
            'invABCD':lambda x: x.Bparams(0, 1, 2, 1),
            'S':lambda x: x.Sparams(0, 1, 2, 1),
            'T':lambda x: x.Tparams(0, 1, 2, 1),
        }
        
        assert network_rep in rep_actions.keys(), f'`{network_rep}` is not 2Port rep'
        
        return rep_actions[network_rep](self.schematic)
#chapteer 2 section 2 rlc_parallel_bandstop filter class
#class with lcapy and skidl subcircuit to create an RLC parallel bandstop filter

class rlc_parallel_bandstop():
    """
    holding class for SkiDl subcircuit and lcapy schematic of a
    bandstop parallel RLC filter primitive
    """
    
    def __init__(self, subcirc_ref=None, L_value=1@u_H, C_value=1@u_F, R_value=1@u_Ohm):
        """
        Args:
            subcirc_ref (str): reference to use for the base of the internal elements
            L_value (float; 1@u_F; Henerys): the inductance in farads for the RLC inductive element
            C_value (float; 1@u_F; Farads): the capacitance in farads for the RLC capacitive element
            R_value (float; 1@u_Ohm; Ohms): the resistance in ohms for the RLC resistive element

        Returns:
            None
        TODO:
            -add assertions
        """
        #add assertions
        self.subcirc_ref=subcirc_ref
        self.L_value=L_value
        self.C_value=C_value
        self.R_value=R_value
        

    @subcircuit
    def SKiDl(self, term_0, term_1, term_2, term_3, return_elements=False):
        """
        Terminals:
        term_0, term_1, term_2, term_3
        
        Terminals are defined via:
        ```
        Left_Termanals - RLC_s_bandstop - Right_Termanals
                             +----+  
        Postive V_i   term_0-|0  2|-term_2     Postive V_o
        Negtive V_i   term_1-|1  3|-term_3     Negtive V_o
                             +----+
        ```
        
        Args:
            return_internls (bool; False): If True return out the internal Voltage Source,
                and Resistance objects in this package
        Returns:
            Returns elements to circuit RLC parallel bandstop filter part element object and if `return_internls`
            is True will return the internal voltage and resistance objects in that order 
        """
        if self.subcirc_ref!=None:
            Lref={'ref':f'L_{self.subcirc_ref}'}
            Cref={'ref':f'C_{self.subcirc_ref}'}
            Rref={'ref':f'R_{self.subcirc_ref}'}
        else:
            Lref={}
            Cref={}
            Rref={}
        
        self.l=L(value=self.L_value, **Lref)
        self.c=C(value=self.C_value, **Cref)
        self.r=R(value=self.R_value, **Rref)
        
        self.c['p', 'n']+=term_0, term_2
        self.l[1, 2]+=term_0, term_2
        self.r[1, 2]+=term_2, term_3
        self.r[2]+=term_1
        
        
        
        if return_elements:
            return self.c, self.r, self.l
    
    def lcapy_self(self, draw_me=True, with_values=True):
        """
        Creates a lcapy schematic of this classes filter that
        can be used for amongst other things: draw a basic schematic
        of this class filter, extract the transfer function, 
        exstract the 2Port Repersntation
        
        Args:
            draw_me (bool): will draw a schematic of this classes filter schematic
            
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
        Return:
            the lcapys circuit object is stored in `self.schematic` and will
            have abstract sympy variable for the elements if `with_values` is False
            and will draw the schematic of just this classes filter if `draw_me` is True
        
        TODO:
            - get the Vin statement into the schematic
        
        """
        
        self.schematic=kiwi.Circuit()
        self.schematic.add('W 0 0_1; right')
        self.schematic.add('W 1 3_1; right=3.5')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        self.schematic.add('W 0_1, 0_2; up=.3')
        self.schematic.add(f'L 0_2 2_3; right')
        self.schematic.add('W 2_3 2_2; down=.3')
        
        self.schematic.add('W 0_1, 0_3; down=.3')
        self.schematic.add(f'C 0_3 2_4; right')
        self.schematic.add('W 2_4 2_2; up=.3')

        self.schematic.add(f'R 2_1 3_1; down')

        
        self.schematic.add('W 2_2 2_1; right')
        self.schematic.add('W 2_1 2; right')
        self.schematic.add('W 3_1 3; right')
        self.schematic.add('P2 2 3; down, v=V_o')
        
        if with_values:
            self.schematic=self.schematic.subs({'R':self.R_value, 'L':self.L_value, 'C':self.C_value})
        
        if draw_me:
            self.schematic.draw()
    
    def get_tf(self, with_values=True, ZPK=True):
        """
        will extract the symbolic transfer function for this filter
        
        Args:
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
            ZPK (bool): if True will try to return the TF in zero-pole-gain form
                else will try to return the TF in canonical form 
                where the unity coefficient is the highest power of the denominator
        
        """
        self.lcapy_self(draw_me=False, with_values=with_values)
        self.tf=self.schematic.transfer(0, 1, 2, 3)
        
        if ZPK:
            return self.tf.ZPK()
        else:
            return self.tf.canonical()
        
    def get_twoPort(self, network_rep='Y', with_values=True):
        """
        Gets the 2Port network representation of this filter.
        The 2Port representation can be controlled
        
        Args:
            network_rep (str; 'Y'): control string for what
                representation is used to choose from:
                
                'Z': two-port Z-parameters matrix
                'Y': two-port Y-parameters matrix
                'H': two-port H-parameters matrix
                'G': two-port G-parameters matrix
                'ABCD': two-port A-parameters matrix
                'invABCD': two-port B-parameters matrix
                'S': two-port S-parameters matrix
                'T': two-port T-parameters matrix
                
                
        """
        
        self.lcapy_self(draw_me=False, with_values=with_values)
        
        #create an action dict to get the 2P rep
        rep_actions={
            'Z':lambda x: x.Zparams(0, 1, 2, 1),
            'Y':lambda x: x.Yparams(0, 1, 2, 1), 
            'H':lambda x: x.Hparams(0, 1, 2, 1),
            'G':lambda x: x.Gparams(0, 1, 2, 1),
            'ABCD':lambda x: x.Aparams(0, 1, 2, 1),
            'invABCD':lambda x: x.Bparams(0, 1, 2, 1),
            'S':lambda x: x.Sparams(0, 1, 2, 1),
            'T':lambda x: x.Tparams(0, 1, 2, 1),
        }
        
        assert network_rep in rep_actions.keys(), f'`{network_rep}` is not 2Port rep'
        
        return rep_actions[network_rep](self.schematic)
#chapteer 2 section 2 lc_balanced_allpass_lowfreq_lattice_filt filter class
#class with lcapy and skidl subcircuit to create an lc all-pass filter

class lc_balanced_allpass_lowfreq_lattice_filt():
    """
    holding class for SkiDl subcircuit and lcapy schematic of a
    bandstop parallel RLC filter primitive
    """
    
    def __init__(self, subcirc_ref=None, L1_value=1@u_H, L2_value=1@u_H,  C1_value=1@u_F, C2_value=1@u_F):
        """
        Args:
            subcirc_ref (str): reference to use for the base of the internal elements
            L1_value (float; 1@u_F; Henerys): the inductance in farads for the top inductive element
            L2_value (float; 1@u_F; Henerys): the inductance in farads for the bottom inductive element
            C1_value (float; 1@u_F; Farads): the capacitance in farads for the term 2 to term 1 capacitive element
            C2_value (float; 1@u_F; Farads): the capacitance in farads for the term 0 to term 2 capacitive element

        Returns:
            None
        TODO:
            -add assertions
        """
        #add assertions
        self.subcirc_ref=subcirc_ref
        self.L1_value=L1_value
        self.L2_value=L2_value
        self.C1_value=C1_value
        self.C2_value=C2_value
        

    @subcircuit
    def SKiDl(self, term_0, term_1, term_2, term_3, return_elements=False):
        """
        Terminals:
        term_0, term_1, term_2, term_3
        
        Terminals are defined via:
        ```
        Left_Termanals - RLC_s_bandstop - Right_Termanals
                             +----+  
        Postive V_i   term_0-|0  2|-term_2     Postive V_o
        Negtive V_i   term_1-|1  3|-term_3     Negtive V_o
                             +----+
        ```
        
        Args:
            return_internls (bool; False): If True return out the internal Voltage Source,
                and Resistance objects in this package
        Returns:
            Returns elements to circuit RLC parallel bandstop filter part element object and if `return_internls`
            is True will return the internal voltage and resistance objects in that order 
        """
        if self.subcirc_ref!=None:
            L1ref={'ref':f'L_{self.subcirc_ref}1'}
            L2ref={'ref':f'L_{self.subcirc_ref}2'}

            C1ref={'ref':f'C_{self.subcirc_ref}1'}
            C2ref={'ref':f'C_{self.subcirc_ref}2'}
        else:
            L1ref={}
            L2ref={}
            
            C1ref={}
            C2ref={}
        
        self.l1=L(value=self.L1_value, **L1ref)
        self.l2=L(value=self.L2_value, **L2ref)

        self.c1=C(value=self.C1_value, **C1ref)
        self.c2=C(value=self.C2_value, **C2ref)
        
        
        term_0+=self.l1[1], self.c2['p']
        term_1+=self.l2[1], self.c1['n']
        term_2+=self.l1[2], self.c1['p']
        term_3+=self.l2[2], self.c2['n']   
        
        
        
        
        if return_elements:
            return self.c, self.r, self.l
    
    def lcapy_self(self, draw_me=True, with_values=True):
        """
        Creates a lcapy schematic of this classes filter that
        can be used for amongst other things: draw a basic schematic
        of this class filter, extract the transfer function, 
        exstract the 2Port Repersntation
        
        Args:
            draw_me (bool): will draw a schematic of this classes filter schematic
            
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
        Return:
            the lcapys circuit object is stored in `self.schematic` and will
            have abstract sympy variable for the elements if `with_values` is False
            and will draw the schematic of just this classes filter if `draw_me` is True
        
        TODO:
            - must draw this better
            - get the Vin statement into the schematic
            
        
        """
        
        self.schematic=kiwi.Circuit()
        
        self.schematic=kiwi.Circuit()
        #self.schematic.add('W 0 0; right')
        #self.schematic.add('W 1 3; right=3.5')
        #self.schematic.add('P1 0 1; down, v=V_i')
        
        #It ant pretty but it works
        self.schematic.add('W 0 0_3; right')
        self.schematic.add(f'L1 0_3 2_2; right')
        self.schematic.add('W 2_2 2; right')

        self.schematic.add('W 0 0_2; rotate=-45')
        self.schematic.add(f'C2 0_2 3; rotate=-45')
        
        self.schematic.add(f'C1 2 1_2; rotate=225')
        self.schematic.add(f'W 1_2 1; rotate=225')
        
        self.schematic.add('W 1 1_3; right')
        self.schematic.add(f'L2 1_3 3_2; right=2')
        self.schematic.add('W 3_2 3; right')

        
        if with_values:
            self.schematic=self.schematic.subs({'L1':self.L1_value, 'L2':self.L2_value, 
                                                'C1':self.C1_value, 'C2':self.C2_value})
        
        if draw_me:
            self.schematic.draw()
    
    def get_tf(self, with_values=True, ZPK=True):
        """
        will extract the symbolic transfer function for this filter
        
        Args:
            with_values (bool): will push the filters element values into
                the resulting lcapy circuit object in place of abstract sympy variables
            
            ZPK (bool): if True will try to return the TF in zero-pole-gain form
                else will try to return the TF in canonical form 
                where the unity coefficient is the highest power of the denominator
        
        """
        self.lcapy_self(draw_me=False, with_values=with_values)
        self.tf=self.schematic.transfer(0, 1, 2, 3)
        
        if ZPK:
            return self.tf.ZPK()
        else:
            return self.tf.canonical()
        
    def get_twoPort(self, network_rep='Y', with_values=True):
        """
        Gets the 2Port network representation of this filter.
        The 2Port representation can be controlled
        
        Args:
            network_rep (str; 'Y'): control string for what
                representation is used to choose from:
                
                'Z': two-port Z-parameters matrix
                'Y': two-port Y-parameters matrix
                'H': two-port H-parameters matrix
                'G': two-port G-parameters matrix
                'ABCD': two-port A-parameters matrix
                'invABCD': two-port B-parameters matrix
                'S': two-port S-parameters matrix
                'T': two-port T-parameters matrix
                
                
        """
        
        self.lcapy_self(draw_me=False, with_values=with_values)
        
        #create an action dict to get the 2P rep
        rep_actions={
            'Z':lambda x: x.Zparams(0, 1, 2, 1),
            'Y':lambda x: x.Yparams(0, 1, 2, 1), 
            'H':lambda x: x.Hparams(0, 1, 2, 1),
            'G':lambda x: x.Gparams(0, 1, 2, 1),
            'ABCD':lambda x: x.Aparams(0, 1, 2, 1),
            'invABCD':lambda x: x.Bparams(0, 1, 2, 1),
            'S':lambda x: x.Sparams(0, 1, 2, 1),
            'T':lambda x: x.Tparams(0, 1, 2, 1),
        }
        
        assert network_rep in rep_actions.keys(), f'`{network_rep}` is not 2Port rep'
        
        return rep_actions[network_rep](self.schematic)

#chapteer 2 section 4 pz_ease class
#class to perform .pz simulations with a bit more grace 
#with some additional built-in analysis tools

class pz_ease(ac_representation_tool, eecomplex_plot_templets):
    def __init__(self, circ_netlist_obj, K=1.0):
        """
        Class to perform Pole Zero (.pz) SPICE simulation with grace
        
        Args:
            circ_netlist_obj (pyspice.Spice.Netlist.Circuit): the Netlist circuit produced 
                from SKiDl's `generate_netlist()`
            
            K (float/int; 1): the gain; must be manually put in or found from .tf analysis
        
        Returns: 
            
        """
        self.circ_netlist_obj=circ_netlist_obj
        
        assert (type(K)==float) or (type(K)==int), 'K must be a float or int'
        self.K=K
        
        
        #dic of allowed pz control statements
        self.allowed_control_statments={'voltage':'vol', 'current':'cur',
                                       'pole-zero':'pz', 'zeros':'zer', 'poles':'pol'}

    
    def pz_def_ports(self, port_0_pos_term, port_0_neg_term, port_1_pos_term, port_1_neg_term, display_table=False):
        """
        Method to set the Port terminals for the two-port section of the circuit under test
        where all inputs must be nodes in the circuit under test
        
        Terminals:
        port_0_pos_term, port_0_neg_term, port_1_pos_term, port_1_neg_term
        
        Port & Terminals are defined via:
        ```
                        Left_Port - Two-Port Section under Test - Right_Port
                                        +-------------+  
        Postive Port0   port_0_pos_term-| DUT Section |-port_1_pos_term     Postive Port1
        Negtive Port0   port_0_neg_term-|             |-port_1_neg_term     Negtive Port1
                                        +-------------+
        ```
        Args:
            display_table (bool; False): when true will display the generated `self.control_df` below
                this method call in a jupyter notebook like environment
            
        
        Returns:
            Settings are recoded in `self.control_df` rows: `'port_0_terms+-'` & `'port_1_terms+-'`
        """
        
        assert port_0_pos_term in self.circ_netlist_obj.node_names, f'`{port_0_pos_term}` is not a node in the circuit under test'
        self.port_0_pos_term=port_0_pos_term
        
        assert port_0_neg_term in self.circ_netlist_obj.node_names, f'`{port_0_neg_term}` is not a node in the circuit under test'
        self.port_0_neg_term=port_0_neg_term
        
        assert port_1_pos_term in self.circ_netlist_obj.node_names, f'`{port_1_pos_term}` is not a node in the circuit under test'
        self.port_1_pos_term=port_1_pos_term
        
        assert port_1_neg_term in self.circ_netlist_obj.node_names, f'`{port_1_neg_term}` is not a node in the circuit under test'
        self.port_1_neg_term=port_1_neg_term
        
        #record the results in table
        self._build_control_table(display_table)
        
    
    def pz_mode_set(self, tf_type='voltage', pz_acu='pole-zero', display_table=False):
        """
        Method to set the pole-zero analysis controls
        
        Args:
            tf_type (str; 'voltage'): the tf for wich the poles and zeros fit to
                if `voltage` the tf is of the form V_o/V_i else if `current` in the form of
                V_o/I_i
            
            pz_acu (str; 'pole-zero'): if `pole-zero` will attempt to get all the poles and zeros for the
                specfied transfer function; else if `zeros` or `poles` will get just the respective zeros
                or poles 
            
            display_table (bool; False): when true will display the generated `self.control_df` below
                this method call in a jupyter notebook like environment
        
        Returns:
            Settings are recoded in `self.control_df` rows: `'tf_type'` & `'acqui_mode'`
        
        
        """
        assert tf_type in self.allowed_control_statments.keys(), f'`{tf_type}` is not `voltage` or `current`'
        self.tf_type=tf_type
        
        assert pz_acu in self.allowed_control_statments.keys(), f'`{pz_acu}` is not `pole-zero` or `poles` or `zeros`'
        self.pz_acu=pz_acu
        
        #record the results in table
        self._build_control_table(display_table)
        
    def _build_control_table(self, display_table=True):
        """
        Internal method to build a pz control table to display pz simulation settings
        
        Args:
            display_table (bool; True): when true will display the generated `self.control_df` below
                this method call in a jupyter notebook like environment
        
        Returns:
            creates dataframe table `self.control_df` that records pz simulation controls
            if `display_table` is true will force showing under jupyter notebook cell
            
        """
        
        self.control_df=pd.DataFrame(columns=['value'], 
                                     index=['tf_type', 
                                            'acqui_mode',
                                            'port_0_terms+-',
                                            'port_1_terms+-'
                                           ])
        if hasattr(self, 'tf_type'):
            self.control_df.at['tf_type']=self.tf_type
        
        if hasattr(self, 'pz_acu'):
            self.control_df.at['acqui_mode']=self.pz_acu
        
        if hasattr(self, 'port_0_pos_term') and hasattr(self, 'port_0_neg_term') :
            self.control_df.at['port_0_terms+-', 'value']=[self.port_0_pos_term,  self.port_0_neg_term]
        
        if hasattr(self, 'port_1_pos_term') and hasattr(self, 'port_1_neg_term') :
            self.control_df.at['port_1_terms+-', 'value']=[self.port_1_pos_term,  self.port_1_neg_term]
        
        self.control_df.index.name='pz_sim_control'
        
        if display_table:
            display(self.control_df)
        
    
    def do_pz_sim(self, display_table=False):
        """
        Method to perform the pole-zero simulation based on values stored in self.control_df
        If the simulation does not converge will give a warning with a basic debug action
        but will set `self.pz_values` to empty dict.
        
        TODO:
            - add simulation kwargs
            - flush out exception handling
        """
        
        attriputs_to_check=['port_0_pos_term', 'port_0_neg_term', 'port_1_pos_term', 'port_1_neg_term', 
                           'tf_type', 'pz_acu']
        
        for i in attriputs_to_check:
            if hasattr(self, i):
                pz_is_go=True
            else:
                pz_is_go=False
                warnings.warn(f'{i} has not been set; pole-zero simulation will not procdede till set')
                
        if pz_is_go:
            self.sim=self.circ_netlist_obj.simulator()
            #I cant catch the warning when it hangs so going to have to do this
            self.pz_values={}

            try:
                self.pz_values=self.sim.polezero(
                    node1=self.port_0_pos_term, 
                    node2=self.port_0_neg_term, 
                    node3=self.port_1_pos_term, 
                    node4=self.port_1_neg_term, 
                    tf_type=self.allowed_control_statments[self.tf_type], 
                    pz_type=self.allowed_control_statments[self.pz_acu]
                )
                
                self._record_pz_results(display_table)
            
            except pspice.Spice.NgSpice.Shared.NgSpiceCommandError:
                self.pz_values={}
                warnings.warn("""PZ analysis did not converge with the current setting:
                start by changing the tf type (self.tf_type) and pz acusisiton type (self.pz_acu) """)
            
    
    def _record_pz_results(self, display_table=True):
        """
        Internal method to record the PZ results to a dataframe
        
        Args:
            display_table (bool; True): when true will display the generated `self.control_df` below
                this method call in a jupyter notebook like environment
        
        Returns:
            creates dataframe table `self.pz_results_DF` that records pz simulation results
            if `display_table` is true will force showing under jupyter notebook cell
            
        """
        self.pz_results_DF=pd.DataFrame(columns=['Type', 'Values'])
        
        if hasattr(self.pz_values, 'nodes'):
            for k, v in self.pz_values.nodes.items():
                self.pz_results_DF.at[len(self.pz_results_DF)]=k, v.as_ndarray()[0]
        
        if display_table:
            display(self.pz_results_DF)
            
        
    def get_pz_sym_tf(self, dec_round=None, overload_K=None):
        """
        Method to get the symbolic transfer function via lacpy
        
        Args:
            dec_round (int; None): contorl to `np.around`'s `decimals` argument
                if left `None` np.around will not be used
            
            overload_K (float/int; None): if not `None` will overload the DC
                gain constant stored in `self.K`
        
        Returns:
            if `self.pz_results_DF` exists return the symbolic transfer function in the s
            dominan in `self.sym_tf`
        """
        if overload_K!=None:
            assert (type(overload_K)==float) or (type(overload_K)==int), 'K must be a float or int'
            self.K=overload_K
        
        if hasattr(self, 'pz_results_DF')!=True:
            warnings.warn('no poles/zero recorded run `self.do_pz_sim`')
        else:
            zeros_B=np.empty(0)
            poles_A=np.empty(0)
            for index, row in self.pz_results_DF.iterrows():
                if 'zero' in row['Type']:
                    zeros_B=np.hstack((zeros_B, row['Values']))
                elif 'pole' in row['Type']:
                    poles_A=np.hstack((poles_A, row['Values']))
            
            if dec_round!=None:
                zeros_B=np.around(zeros_B, dec_round)
                poles_A=np.around(poles_A, dec_round)
            
            self.zeros_B=zeros_B; self.poles_A=poles_A
            #wish I didn't have to do this
            zeros_B=zeros_B.tolist(); poles_A=poles_A.tolist()  
            
            #use lcapy to get the symbolic tf 
            self.sym_tf=kiwi.zp2tf(zeros_B, poles_A, K=self.K)
            #use simplify because if in pzk it does weird things with j that
            #lambdfy has issues with
            self.sym_tf=self.sym_tf.simplify()
            
    def plot_pz_loc(self, ax=None, title='', unitcircle=False):
        """
        uses lcapy's `plot_pole_zero` in https://github.com/mph-/lcapy/blob/6e42983d6b77954e694057d61045bd73d17b4616/lcapy/plot.py#L12
        to plot the poles and zero locations on a Nyquist chart
        
        Args:
            axs (list of matplotlib axis; None): If left None will create a new plot, else must 
                be a list of matplotlib subplots axis to be added to where the first entry
                will be the magnitude axis, and the second will be the phase axis
            
            title (str; ''): Subplot title string
            
            unitcircle (bool; False): when True will plot the unit circle on the resulting plot
        
        Returns:
            Returns a real-imag with Poles and Zero map, and if an axis was passed to `ax` will be modified
            with the pole-zero map
        
        
        """
        axs=ax or plt.gca()
        if hasattr(self, 'sym_tf')==False:
            warnings.warn("""Trying to get symbolic transfer function from `self.get_pz_sym_tf`
            thus you will get what you get""")
            self.get_pz_sym_tf()
        
        self.sym_tf.plot(axes=axs, unitcircle=unitcircle, 
                         #wish there be a better way to do this
                        label="'X'=pole; 'O'=zero")
        
        axs.axhline(0, linestyle='--', linewidth=2.0, color='black')
        axs.axvline(0, linestyle='--', linewidth=2.0, color='black')
        axs.set_xlabel('Real'); axs.set_ylabel('Imag')
        
        axs.legend()
        
        if title!='':
            title=' of '+title
        axs.set_title(f'Pole-Zero locations plot{title}');

    
    
    def plot_3d_laplce(self, title=''):
        """
        Creates 3d plots of the Laplace space of the transfer function, one for the mag
        the other for the phase in degrees unwrapped
        
        Args:
                title (str; ''): Subplot title string
        
        Returns:
            returns a 3d plot with the left subplot being the  mag and the right being the phase
        
        TODO:
            - get the freaking color bar into a clean location when working with 3d plots
            - merge phase as color into mag see the physics video by eugene on Laplace 
        """
        if hasattr(self, 'sym_tf')==False:
            warnings.warn("""Trying to get symbolic transfer function from `self.get_pz_sym_tf`
            thus you will get what you get""")
            self.get_pz_sym_tf()
        
        #import the additnal matplotlib featuers for 3d
        from mpl_toolkits.mplot3d import Axes3D
        from matplotlib import cm
        
        #stole this off lcapy's plot_pole_zero
        #https://github.com/mph-/lcapy/blob/7c4225f2159aa33398dac481041ed538169b7058/lcapy/plot.py
        
        
        #check self.sys_tf is good to be used
        sys_tf_syms=self.sym_tf.symbols
        assert len(sys_tf_syms)==1 and ('s' in sys_tf_syms.keys()), 'trasfer function must be laplce form and only have `s` as a free symbol'

        #lambdfy the tf
        sys_tf_lam=sym.lambdify(kiwi.s, self.sym_tf.canonical(), 'numpy', dummify=False)

        #get the plot bounds
        #stole this off lcapy's plot_pole_zero
        #https://github.com/mph-/lcapy/blob/7c4225f2159aa33398dac481041ed538169b7058/lcapy/plot.py

        poles = self.sym_tf.poles()
        zeros = self.sym_tf.zeros()
        try:
            p = np.array([p.cval for p in poles.keys()])
            z = np.array([z.cval for z in zeros.keys()])
        except ValueError:
            raise TypeError('Cannot get poles and zeros of `self.sym_tf')
        a = np.hstack((p, z))
        x_min = a.real.min()
        x_max = a.real.max()
        y_min = a.imag.min()
        y_max = a.imag.max()

        x_extra, y_extra = 3.0, 3.0

        # This needs tweaking for better bounds.
        if len(a) >= 2:
            x_extra, y_extra = 0.1 * (x_max - x_min), 0.1 * (y_max - y_min)
        if x_extra == 0:
            x_extra += 1.0
        if y_extra == 0:
            y_extra += 1.0

        x_min -= 0.5 * x_extra
        x_max += 0.5 * x_extra
        y_min -= 0.5 * y_extra
        y_max += 0.5 * y_extra

        #the input domain
        RealRange=np.linspace(x_min, x_max, 100); ImagRange=np.linspace(y_min, y_max, 100)
        sr, si=np.meshgrid(RealRange, ImagRange)
        s_num=sr+1j*si

        #plot this
        fig = plt.figure()
        
        #mag 3d plot
        ax3d_mag = fig.add_subplot(121, projection='3d')

        XmagPlot=ax3d_mag.plot_surface(sr, si, np.abs(sys_tf_lam(s_num)), alpha=0.5,
                                 cmap=cm.coolwarm, antialiased=False)
        ax3d_mag.set_xlabel(r'$\sigma$'); ax3d_mag.set_ylabel(r'$j\omega$'), ax3d_mag.set_zlabel(r'$|X|$')
        fig.colorbar(XmagPlot, shrink=0.5, aspect=5)

        #phase 3d plot
        ax3d_phase = fig.add_subplot(122, projection='3d')

        XphasePlot=ax3d_phase.plot_surface(sr, si, angle_phase_unwrap(sys_tf_lam(s_num)), alpha=0.5,
                                 cmap=cm.coolwarm, antialiased=False)
        ax3d_phase.set_xlabel(r'$\sigma$'); ax3d_phase.set_ylabel(r'$j\omega$'), ax3d_phase.set_zlabel(r'$ang(X)$')
        fig.colorbar(XphasePlot, shrink=0.5, aspect=5)

        plt.tight_layout()
        
        if title!='':
            title=' of '+title
        ax3d_mag.set_title(f'3D Mag Laplace plot{title}');
        ax3d_phase.set_title(f'3D Phase_deg Laplace plot{title}');

    
    
    def scipy_pzk(self, dec_round=None, overload_K=None):
        """
        Method to create to generate the Numerator (a) and Denominator (b)
        arrays to use within the scipy signal framework see 
        https://docs.scipy.org/doc/scipy/reference/signal.html
        &
        https://docs.scipy.org/doc/scipy/reference/generated/scipy.signal.zpk2tf.html#scipy.signal.zpk2tf
        
        Args:
            dec_round (int; None): contorl to `np.around`'s `decimals` argument
                if left `None` np.around will not be used
            
            overload_K (float/int; None): if not `None` will overload the DC
                gain constant stored in `self.K`
        
        Returns:
            if `self.pz_results_DF` exists tries to return the coefficients for an lti
            filter for use in scipy's signal library in a dictionary in `self.scipy_tfcoef`
        """
        if overload_K!=None:
            assert (type(overload_K)==float) or (type(overload_K)==int), 'K must be a float or int'
            self.K=overload_K
        
        if hasattr(self, 'pz_results_DF')!=True:
            warnings.warn('no poles/zero recorded run `self.do_pz_sim`')
        else:
            zeros_B=np.empty(0)
            poles_A=np.empty(0)
            for index, row in self.pz_results_DF.iterrows():
                if 'zero' in row['Type']:
                    zeros_B=np.hstack((zeros_B, row['Values']))
                elif 'pole' in row['Type']:
                    poles_A=np.hstack((poles_A, row['Values']))
            
            if dec_round!=None:
                zeros_B=np.around(zeros_B, dec_round)
                poles_A=np.around(poles_A, dec_round)
                
            b, a=scipy_zpk2tf(zeros_B, poles_A, self.K)
            self.scipy_tfcoef={'b':b, 'a':a}
    
    def get_sym_freq_resp(self, freq_vec=None):
        """
        method to get the symbolic transfer function  response 
        
        Args:
            freq_vec (numpy array or pandas series; None): frequencies
            to generate results from the generated symbolic transfer function
            if `None` will come up with a freancy vector inside via `np.logspace(-1, 12, 12*10*10)`
        
        Returns:
            will store frequency vector as the index in `self.pz_symresults_DF`
            where the response of the symbolic transfer function will be stored
            in `self.pz_symresults_DF['pzk_symres_[V]']` from there will
            use the inheritance from `ac_representation_tool` to then 
            generate `self.ac_sim_real_DF['pzk_symres_[V]'], self.ac_sim_imag_DF['pzk_symres_[V]']`, ect
            
        """
        if hasattr(self, 'sym_tf')==False:
            warnings.warn("""Trying to get symbolic transfer function from `self.get_pz_sym_tf`
            thus you will get what you get""")
            self.get_pz_sym_tf()
            
        if freq_vec==None:
            self.freq_vec=np.logspace(-1, 12, 12*10*10)
        
        self.pz_symresults_DF=pd.DataFrame(index=self.freq_vec)
        self.pz_symresults_DF.index.name='freq[Hz]'
        
        self.pz_symresults_DF['pzk_symres_[V]']=self.sym_tf.frequency_response(self.pz_symresults_DF.index.values).astype('complex')
        
        ac_representation_tool.__init__(self, self.pz_symresults_DF)
        
        #just make the reps
        self.make_real_imag()
        self.make_mag_phase()
    
    def plot_nyquist_with_pz(self, ax=None, title='', unitcircle=False):
        """
        plotting utility to plot the pole-zero location on top of a Nyquist
        a plot of the response of the symbolic transfer function
        
        Args:
            ax (matplotlib axis; None): If left None will create a new plot, else must
                be a matplotlib subplot axis to be added to
            
            title (str; ''): Subplot title string
            
        Returns:
            Returns a Nyquist plot with the pole-zero locations superimposed on top of, 
            and if an axis was passed to `ax` will be modified
            
            
        
        """
        axs=ax or plt.gca()
        if hasattr(self, 'ac_sim_real_DF')==False:
            warnings.warn("""Trying to get the real imag data from `self.get_sym_freq_resp`, 
            you get what you get""")
            self.get_sym_freq_resp()
        
        
        self.nyquist_plot_templet(self.ac_sim_real_DF['pzk_symres_[V]'], self.ac_sim_imag_DF['pzk_symres_[V]'], 
                           ax=axs)
        self.plot_pz_loc(ax=axs, unitcircle=unitcircle)

        
        if title!='':
            title=' of '+title
        axs.set_title(f'Nyquist with Pole-Zero locations plot{title}');
    
    def plot_bode_from_pz(self, ax=None, title=''):
        """
        plotting utility to plot single mag and phase on a single bode plot 
        for the response of the symbolic transfer function
        
        Args:
            ax (matplotlib axis; None): If left None will create a new plot, else must
                be a matplotlib subplot axis to be added to
            
            title (str; ''): Subplot title string
            
        Returns:
            Returns a bode plot from the symbolic transfer function, 
            and if an axis was passed to `ax` will be modified
        """
        axs=ax or plt.gca()
        if hasattr(self, 'ac_sim_real_DF')==False:
            warnings.warn("""Trying to get the mag phase data from `self.get_sym_freq_resp`, 
            you get what you get""")
            self.get_sym_freq_resp()
        
        
        
        self.bode_plot_one_templet(self.ac_sim_mag_DF.index, self.ac_sim_mag_DF['pzk_symres_[V][dB]'], self.ac_sim_phase_DF['pzk_symres_[V][deg]'], 
                           ax=axs)
        
        if title!='':
            title=' of '+title
        axs.set_title(f'Bode plot from Pole-Zero anylsis plot{title}');
    
    def plot_nichols_from_pz(self, ax=None, title=''):
        """
        plotting utility to plot the Nichols chart
        for the response of the symbolic transfer function
        
        Args:
            ax (matplotlib axis; None): If left None will create a new plot, else must
                be a matplotlib subplot axis to be added to
            
            title (str; ''): Subplot title string
            
        Returns:
            Returns a Nichols chart plot from the symbolic transfer function, 
            and if an axis was passed to `ax` will be modified
        """
        axs=ax or plt.gca()
        
        if hasattr(self, 'ac_sim_real_DF')==False:
            warnings.warn("""Trying to get the mag phase data from `self.get_sym_freq_resp`, 
            you get what you get""")
            self.get_sym_freq_resp()
        
        
        
        self.nichols_plot_templet(self.ac_sim_mag_DF['pzk_symres_[V][dB]'], self.ac_sim_phase_DF['pzk_symres_[V][deg]'], 
                           ax=axs)
        
        if title!='':
            title=' of '+title
        axs.set_title(f'Nichols plot from Pole-Zero anylsis plot{title}');
    

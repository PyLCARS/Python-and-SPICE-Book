#Library import statements

from skidl.pyspice import *
#can you say cheeky 
import PySpice as pspice
#becouse it's written by a kiwi you know
import lcapy as kiwi

import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import warnings

from IPython.display import YouTubeVideo, display

import traceback

#chapter 1 section 1 op_results_collect class
#used for basic operating point simulation analysis, does not support internal parameters

class op_results_collect():
    """
    Basic class to get pyspice operating_point (ngspice `.op`) simulation results
    for simulated circuit's global elements, no internal parameters just yet, that's coming
    """
    
    def __init__(self, op_sim_circ, display_results=False):
        """
        Basic class to get pyspice operating_point (ngspice `.op`) simulation results
        for simulated circuit's global elements
        
        Args:
            op_sim_circ (pspice.Spice.Netlist.Circuit): the Netlist circuit produced 
            from SKiDl's `generate_netlist()`
            
            display_results (bool; False): option to have the simulation results
                stored in `self.results_df` automatic displayed from a jupyter notebook
                ell
        
        Returns:
            will create a simulation in `self.op_sim`, raw results of dc operating 
            point simulation will be stored in `self.op_sim_results`, the tablized
            results will be stored in pandas dataframe `self.results_df`
        
        TODO: 
            - add kwargs to the simulator
            - add an assertion that only a pyspice netlist circuit obj can
                    be passed into op_sim_circ
            
        """
        #need to add assertions for op_sim_circ ==pspice.Spice.Netlist.Circuit
        
        #store the circuit internally
        self.op_sim_circ=op_sim_circ
        #create the sim obj
        self.op_sim=self.op_sim_circ.simulator()
        #do the sim
        self.op_sim_results=self.op_sim.operating_point()
        #get the results in a table
        self.get_values(display_results)
    
    def create_results_df(self):
        """
        Creates an empty dataframe of the results
        Args:
            none
        Returns:
            self.results_df: a pandas dataframe 
        """
        #could be done in self.__init__ but OOP
        #create pandas DF container
        self.results_df=pd.DataFrame(columns=['Circ_Item', 'Item_Type', 'Value', 'Unit'])
    
    def get_values(self, display_results=False):
        """
        Retrieves elements and nodes from the circuit and then gets the corresponding 
        value from the simulation
        
        Args:
            display_results (bool; False): if True will use ipython display to show 
            self.results_df just below the cell
        
        Returns:
            self.results_df: populated dataframe of simulation results
        
        """
        #insure the dataframe is constructed to fill
        self.create_results_df()
        
        #get node voltages
        for n in self.op_sim_circ.node_names:
            #skip ground
            if n=='0':
                continue
            
            #add row and set value
            self.results_df.at[len(self.results_df), :]=[n, 'Node_Voltage', self.op_sim_results[n].as_ndarray()[0], 'V']
        
        #get the current from  any voltage sources
        for cm in self.op_sim_circ.element_names:
            #only look at voltage sources since SPICE uses them for current reading
            if 'V'==cm[0]:
                #the negative sign on the results will be explained in alsdkjfaldsjlk
                self.results_df.at[len(self.results_df), :]=[cm, 'Branch_Curr', -self.op_sim_results[cm].as_ndarray()[0], 'A']
        
        #gets rid of the generic index and replace it with items are index
        self.results_df.set_index('Circ_Item', drop=True, append=False, inplace=True, verify_integrity=False)

        
        #prety display of results
        if display_results:
            display(self.results_df)
        
        return self.results_df
 
 
#chapter 1 section 2 get_skidl_spice_ref function
#used for getting the name of the element as it would appear in a
#generated netlist

def get_skidl_spice_ref(skidle_element):
    """
    Helper function to retrieve SKiDL element name as appears in the final netlist
    
    Args:
        skidle_element (skidl.Part.Part): SKiDl part to get the netlist name from
    
    Returns:
        returns a string with the netlist name of `skidle_element`, or throws an
        error if `skidle_element` is not a SKiDl part
    
    """
    assert repr(type(skidle_element))=="<class 'skidl.part.Part'>", '`skidle_element` must be a SKiDl part'
    
    if skidle_element.ref_prefix!=skidle_element.ref[0]:
        return skidle_element.ref_prefix+skidle_element.ref
    else:
        return skidle_element.ref

#chapter 1 section 2 dc_cs2vs function
#creates a voltage source element to the current source based on the 
#value if the input DC current element and it's parallel resistor
#via the source transformation rules

def dc_cs2vs(dc_cs, cs_par_res):
    """
    Create a new equivalent voltage source to the current source with parallel resistance
    Args:
        dc_cs (SKiDl current source): the current source to transform to a voltage source
        cs_par_res (SKiDl resistor): the parrel resistor to the current source to be transformed
    
    Returns:
        returns an equivalent DC voltage source element to the current source based on the 
        value if the input DC current element and it's parallel resistor via the source transformation rules 
    
    TODO:
        -figure out how to do assertion check that cs_par_res is in parallel to dc_cs
        
    Future:
        -make into subcircuit with net inputs to automatically add the new source and resistance to the circuit
    """
    
    #do assertion checks to insure that passed in objects are what they are supposed to be
    assert dc_cs.ref_prefix=='I', '<dc_cs> was not a current source'
    assert cs_par_res.ref_prefix=='R', '<cs_par_res> was not a resistor'
    
    
    old_maxpower=float((dc_cs.dc_value**2)*cs_par_res.value)
    new_vs_val=float(dc_cs.dc_value*cs_par_res.value)
    assert np.around(float(new_vs_val*dc_cs.dc_value), 6)==np.around(old_maxpower, 6), "Um, something wrong since before and after max power not equal"
    
    new_vs_ref=dc_cs.ref
    if new_vs_ref[0]!='I':
        new_vs_ref='I'+new_vs_ref
    
    
    new_vs_ref=f"V{new_vs_ref[1:]}_f_{new_vs_ref}"
    print(new_vs_ref)
    
    
    eq_vs=V(dc_value=new_vs_val@u_V, ref=new_vs_ref)
    warnings.warn(f"""New voltage source values: {new_vs_val} [V] with max aviabel power {old_maxpower} [W] \n transformed creation statment will be like: \n`(gnd & <eq_vs>['n', 'p'] & <cs_par_res>)`""")
    
    return eq_vs

#chapter 1 section 2 dc_vs2cs function
#creats current source element to the voltage source based on the 
#value if the input DC current element and it's series resistor
#via the source transformation rules

def dc_vs2cs(dc_vs, vs_ser_res):
    """
    Create a new equivalent current source to voltage source with series resistance
    Args:
        dc_vs (SKiDl voltage source): the voltage source to transform to a current source
        vs_ser_res (SKiDl resistor): the serries resistor to the voltage source to be transformed
    
    Returns:
    
    TODO:
        -figure out how to do assertion check that vs_ser_res is in serries to dc_vs
    
    Future:
        -make into subcircuit with net inputs to automatically add the new source and resistance to the circuit
    """
    
    #do assertion checks to insure that passed in objects are what they are supposed to be
    assert dc_vs.ref_prefix=='V', '<dc_vs> was not a voltage source'
    assert vs_ser_res.ref_prefix=='R', '<vs_ser_res> was not a resistor'
    
    old_maxpower=float((dc_vs.dc_value**2)/vs_ser_res.value)
    new_cs_val=float(dc_vs.dc_value/vs_ser_res.value)
    assert np.around(float(new_cs_val*dc_vs.dc_value), 6)==np.around(old_maxpower, 6), "Um, something wrong since before and after max power not equal"


    new_cs_ref=dc_vs.ref
    if new_cs_ref[0]!='V':
        new_cs_ref='V'+new_cs_ref
    
    
    new_cs_ref=f"I{new_cs_ref[1:]}_f_{new_cs_ref}"
    #print(new_cs_ref)
    
    
    eq_cs=I(dc_value=new_cs_val@u_A, ref=new_cs_ref)# might still need this: , circuit=vs_ser_res.circuit)

    warnings.warn(f"""New current source values: {new_cs_val} [A] with max aviabel power {old_maxpower} [W] \n transformed creation statment will be like:\n `(gnd & <eq_cs>['n', 'p'] | <vs_ser_res>)` \n""")
    
    return eq_cs

#chapter 1 section 2 op_internal_ivp class
# class to get both the branch currents and node voltages,
# along with the internal parameters values for 
# DC operating point analysis

class op_internal_ivp():
    """
    Class for creating a SPICE simulation on the whole circuits internal parameters for current, voltage, and power
    for dc operating point (.op) simulations. Will only collect internal parameters and not global voltage and currents 
    of the circuit
    TODO:
        Make this inheritable from op_results_collect 
    """
    def __init__(self, op_sim_circ, display_results=False):
        """
        Basic class to get pyspice operating_point (ngspice `.op`) simulation results
        for internal parameters for Resistors, Current Source, Voltage Source current, voltage, power
        respectively
        
        
        Args:
            op_sim_circ (pspice.Spice.Netlist.Circuit): the Netlist circuit produced 
            from SKiDl's `generate_netlist()`
            
            display_results (bool; False): option to have the simulation results
                stored in `self.results_df` automatically displayed from a jupyter notebook
                ell
        
        Returns:
            will create a simulation in `self.op_sim`, raw results of dc operating 
            point simulation will be stored in `self.op_sim_results`, the tablized
            results will be stored in pandas dataframe `self.results_df`
        
        TODO: 
            - add kwargs to the simulator
            - add an assertion that only a pyspice netlist circuit obj can
                    be passed into op_sim_circ
            
        """
        #need to add assertions for op_sim_circ ==pspice.Spice.Netlist.Circuit
        
        #store the circuit internally
        self.op_sim_circ=op_sim_circ
        #create the sim obj
        self.op_sim=self.op_sim_circ.simulator()
        
        #store bool to display results dataframe
        self.display_results=display_results
        
        #create the internal parameters to save
        self.create_savable_items()
        #run the sim for .op for internal parameters and record results
        self.sim_and_record()
        
    def create_savable_items(self):
        """
        Helper method to create a listing of internal parameters and the table of the results.
        Right now only creates savable internal parameters for:
        Linear Dependent Voltage Sources: current, power
        Linear Dependent Current Sources: current, voltage, power
        Standard Resistor: current, voltage, power
        Linear Dependent Current Sources: current, voltage, power
        VCCS: current, voltage, power
        VCVS: current, voltage, power
        CCVS: current, voltage, power
        CCCS:currrent, voltage, power
        
        See ngspice manual typically chapter 31 "Model and Device Parameters"
        for more deitals about deice intiernal parmaters

        """
        self.results_df=pd.DataFrame(columns=['Circ_Item', 'Item_Type', 'Value', 'Unit'])
        self.results_df.set_index('Circ_Item', drop=True, append=False, inplace=True, verify_integrity=False)
        
        #helper row creation statement
        def add_row(index, unit): 
            self.results_df.at[index, ['Item_Type', 'Unit']]=['internal', unit]
            
        
        for e in self.op_sim_circ.element_names:
            """
            Ref: ngspice documentation chapter 31 (typically): Model and Device Parameters
                
            """
            #resistors
            if e[0]=='R':
                add_row(f'@{e}[i]', 'A')
                add_row(f'@{e}[p]', 'W')
            
            #independnt current source
            elif e[0]=='I':
                add_row(f'@{e}[c]', 'A')
                add_row(f'@{e}[v]', 'V')
                add_row(f'@{e}[p]', 'W')
            
            #independ Voltage source
            elif e[0]=='V':
                add_row(f'@{e}[i]', 'A')
                add_row(f'@{e}[p]', 'W')
            
            #controlled sources
            elif e[0] in ['F', 'H', 'G', 'E']:
                add_row(f'@{e}[i]', 'A')
                add_row(f'@{e}[v]', 'V')
                add_row(f'@{e}[p]', 'W')
            
            else:
                warnings.warn(f"Circ Element {e} is not setup to have internals read, skiping")
                
    def sim_and_record(self):
        """
        run the .op simulation and get the internal values 
        
        Args: None
        Returns:
            `self.internal_opsim_res` store the raw results of the .op for internal pamtyers
            whereas `self.results_df` stores the pandas dataframe of internal parameters results
        
        TODO: figure out how to do this at the same time as the main node branch sim
            this doing separately is getting old
        """
        save_items=list(self.results_df.index)
        self.op_sim.save_internal_parameters(*save_items)
        self.internal_opsim_res=self.op_sim.operating_point()
        
        for save in save_items:
            self.results_df.at[save, 'Value']=self.internal_opsim_res[save].as_ndarray()[0]
        
        if self.display_results:
            print('.op sim internal parmter results')
            display(self.results_df)

#chapter 1 section 3 dc_ease class
# class to perform easily perform DC sweep simulations
#gets both branch currents/node voltages and also internal parameters

class dc_ease():
    """
    Class to perform DC (.dc) SPICE simulation with some grace; only does single variable
    and is currently limited to what +variable pyspice supports
    
    TODO:
        make so that it can iterate over all filled in entries in `self.sweep_DF` and
        exports data to xarray
        
    """
    def __init__(self, circ_netlist_obj):
        """
        Class to perform DC (.dc) SPICE simulation with some grace; only does single variable
        and is currently limited to what variable pyspice supports
        
        Args:
            circ_netlist_obj (pspice.Spice.Netlist.Circuit): the Netlist circuit produced 
                from SKiDl's `generate_netlist()`
        
        Returns: 
            creates a table of variables that could be swept in `self.sweep_DF`
            this table will still need to be filled out for a single variable to be swept
        
        TODO: add kwargs to pass to sim creation

        """
        #need to add assertions for op_sim_circ ==pspice.Spice.Netlist.Circuit
        self.circ_netlist_obj=circ_netlist_obj
        
        self._build_table()
    
    def _build_table(self):
        """
        protected method to create `self.sweep_DF` dataframe of all things in the
        the circuit that can be swept; currently that is only voltage and current sources
        
        TODO:
            -when pyspice accepts more things to sweep add them below
        """
        
        self.sweep_DF=pd.DataFrame(columns=['Element', 'Start', 'Stop', 'Step'])
        
        for e in self.circ_netlist_obj.element_names:
            if e[0] in ['V', 'I']:
                self.sweep_DF.at[len(self.sweep_DF), 'Element']=e
        
        #convert ELmenet to index
        self.sweep_DF.set_index('Element', drop=True, append=False, inplace=True, verify_integrity=False)
    

        
    def clean_table(self):
        """
        Helper method to clean `self.sweep_DF` of any rows that are not to going to be swept
        """
        self.sweep_DF.dropna(axis=0, how='any', inplace=True)
        
    def addback_elements_2table(self):
        """
        Helper method to return empty rows to `self.sweep_DF` to manually then add sweep info to
        uses `self._build_table` to rebuild the table
        """
        for e in self.circ_netlist_obj.element_names:
            if e not in self.sweep_DF.index:
                self.sweep_DF.at[e, :]=np.nan
    
    #so ngspice is only 2d sweeps and so need figurer this out to do 2d sets
    #but will, for now, do sets of 1d
    #def _make_sim_control(self):
    #    self.clean_table()
    #    
    #    #for now dc sim through pyspice only support voltage and current sources
    #    self.dc_contorl={row[0]: slice(row[1], row[2], row[3]) for row in self.sweep_DF.itertuples() if row[0][0] in ['V', 'I']}
    
    
    def _make_sim_control(self, varible):
        """
        Internal method to extract the row information to the .dc kargs and slice argument
        
        Args:
            varible (str): the index in `self.sweep_DF` to exstract the info from
            
        Returns:
            `self.dc_contorl` what is feed into the .dc to do the simulation for a single var
            
        """
        #clean the table
        self.clean_table()
        #get the location of the variable's index
        try:
            row=self.sweep_DF.index.get_loc(varible)
            #get said row
            row=self.sweep_DF.iloc[row]
            #need to fix this in pyspice to make it fully ngspice capale
            if row.name[0] not in ['I', 'V']:
                raise KeyError
            
            #build the control
            self.dc_contorl={row.name:slice(row.Start, row.Stop, row.Step)}
            
            return 0
            
        except KeyError:
            print(f' Variable: {varible} is not in or able self.sweep_DF to be swept ')
            
            return 1
    
    def do_dc_sim(self, varible):
        """
        Does a standard Branch and Node .dc simulation for a single filled out row in `self.sweep_DF`
        
        Args:
            variable (str): the index in `self.sweep_DF` that corresponds to the thing to sweep in the .dc simulation
        
        Returns: 
            raw results are stored in `self.dc_vals`, processed results are automatically stored in 
            `self.dc_resultsNB_DF` via `self.record_dc_nodebranch`
        """
        
        
        if self._make_sim_control(varible)==0:
            self.sim=self.circ_netlist_obj.simulator()
            self.dc_vals=self.sim.dc(**self.dc_contorl)
            
            self.record_dc_nodebranch(varible)
        
        else:
            pass
                
        
    
    def record_dc_nodebranch(self, varible):
        """ 
        Helper method to put .dc node branch results into a dataframe where the index is the  variable
        
        Args:
            variable (str): used to set the index name
        
        Returns:
            `self.dc_resultsNB_DF` which is a pandas dataframe with the index being the sweep and the columns being
            the node voltages and branch currents from any available voltage sources
            
        """
        self.dc_resultsNB_DF=pd.DataFrame(index=self.dc_vals.sweep.as_ndarray())
        
        self.dc_resultsNB_DF.index.name=varible
        
        #get node voltages:
        
        for n in self.circ_netlist_obj.node_names:
            if n=='0':
                continue
            self.dc_resultsNB_DF[n+'_[V]']=self.dc_vals[n].as_ndarray()
        
        #get the current from any voltage sources:
        for cm in self.circ_netlist_obj.element_names:
            if 'V'==cm[0]:
                self.dc_resultsNB_DF[cm+'_[A]']=-self.dc_vals[n].as_ndarray()
    
    def do_dc_intsim(self, varible):
        """
        Does a .dc simulation the internal variables for  a single filled out row in `self.sweep_DF`
        Right now the internal variables gathered are:
        current, voltage power for: VCCS VCVS, CCVS, CCCS, Independent current source 
        current and power for Resistors and independent voltage sources
        
        Args:
            variable (str): the index in `self.sweep_DF` that corresponds to the thing to sweep in the .dc simulation
        
        Returns: 
            raw results are stored in `self.dc_resultsINT_DF`, processed results are automatically stored in 
            `self.dc_resultsNB_DF` via `self.record_dc_internals`
        """
        
        
        self.save_internals=[]
        
        for e in self.circ_netlist_obj.element_names:
            if e[0]=='R':
                self.save_internals+=[f'@{e}[i]', f'@{e}[p]']
            
            elif e[0]=='I':
                self.save_internals+=[f'@{e}[c]', f'@{e}[v]', f'@{e}[p]']
            
            elif e[0]=='V':
                self.save_internals+=[f'@{e}[i]', f'@{e}[p]']
            
            elif e[0] in ['F', 'H', 'G', 'E']:
                self.save_internals+=[f'@{e}[i]', f'@{e}[v]', f'@{e}[p]']
        
        if self._make_sim_control(varible)==0:
            self.sim=self.circ_netlist_obj.simulator()
            self.sim.save_internal_parameters(*self.save_internals)
            self.dc_vals=self.sim.dc(**self.dc_contorl)
            
            self.record_dc_internals(varible)
    
    def record_dc_internals(self, varible):
        """ 
        Helper method to put .dc internal variable results into a dataframe where the index is the  variable
        
        Args:
            variable (str): used to set the index name
        
        Returns:
            `self.dc_resultsINT_DF` which is a pandas dataframe with the index being the sweep and the columns being
            the saved internal variables
            
        """
        self.dc_resultsINT_DF=pd.DataFrame(index=self.dc_vals.sweep.as_ndarray())
        self.dc_resultsINT_DF.index.name=varible

        
        for i in self.save_internals:
            if i[1]=='I' and i[-2]=='c':
                result_name=i[1:]+'_[A]'
            
            elif i[-2]=='i':
                result_name=i[1:]+'_[A]'
            
            elif i[-2]=='v':
                result_name=i[1:]+'_[V]'
            
            elif i[-2]=='p':
                result_name=i[1:]+'_[W]'
            
            #deal with the fliped current for voltage sources even internally
            if i[1]=='V' and i[-2]=='i':
                self.dc_resultsINT_DF[result_name]=-self.dc_vals[i].as_ndarray()
            else:
                self.dc_resultsINT_DF[result_name]=self.dc_vals[i].as_ndarray()
            
    
    def quick_plot(self, nodebranch_int_cont='nb'):
        """
        Creates a quick plot of all columns with respect to the index for the data stored in 
        `self.dc_resultsNB_DF` or `self.dc_resultsINT_DF`. Note that for larger cirucits
        these plots can have a ridicules number of subplots so yeah there is a reason for this
        methods name
        
        Args:
            nodebranch_int_cont (str; 'nb'): control word for wither to plot the node branch data ('nb')
                or the internal data ('int')
        
        Returns:
            Creates a rudimentary plot sharing a common x axis correspond the index
            for the data source and subplot row for each column in the data source
        
        """
        
        assert nodebranch_int_cont in ['nb', 'int'], f'{nodebranch_int_cont} in not a allowed control statment'
        #add a check for existence
        
        if nodebranch_int_cont=='nb':
            plotDF=self.dc_resultsNB_DF
        elif nodebranch_int_cont=='int':
            plotDF=self.dc_resultsINT_DF
        
        
        plotDF.plot(subplots=True, sharex=True, grid=True)
        # won't work for a data source with lots of columns
        plt.tight_layout()

#chapter 1 section 3 real_dcVs subcircuit function
# SKiDl subcirucit to create a dc voltage source with 
# added series resistor


@subcircuit
def real_dcVs(global_ref, pos_term, neg_term, starting_V=1@u_V, starting_R=50@u_Ohm, 
           return_internls=False):
    """
    SKiDl subcircuit to create a simple non-ideal DC voltage source
    
    Args:
        global_ref (str): reference to use for the base of the internal elements
        
        pos_term (SKiDl net or pin): positive terminal of the nonideal voltage source
            to connect to the rest of the circuit
        
        neg_term (SKiDl net or pin): negative terminal of the nonideal voltage source
            to connect to the rest of the circuit 
        
        starting_V (float; 1; Volts):the intial DC voltage to set the internal ideal
            the voltage source in this package to
        
        starting_R (float; 50; Ohm): the initial resistance to set the internal
            serial resistance to the ideal voltage source in this subcircuit to
        
        return_internls (bool; False): If True return out the internal Voltage Source,
            and Resistance objects in this package
    
    Returns:
        Returns it's self a SKiDl part element object and if `return_internls`
        is True will return the internal voltage and resistance objects in that order 
    """
    vs=V(ref=f'V_{global_ref}', dc_value=starting_V)
    rs=R(ref=f'R_{global_ref}', value=starting_R)
    
    vs['p', 'n']+=rs[1], neg_term
    rs[2]+=pos_term
    
    if return_internls:
        return vs, rs

#chapter 1 section 3 real_dcIs subcircuit function
# SKiDl subcirucit to create a dc current source with 
# added parallel resistor

@subcircuit
def real_dcIs(global_ref, pos_term, neg_term, starting_I=1@u_A, starting_R=50@u_Ohm, 
           return_internls=False):
    """
    SKiDl subcircuit to create a simple non-ideal DC current source.
    Where the positive terminal is to the top of the arrow of a schematically drawn
    current source
    
    Args:
        global_ref (str): reference to use for the base of the internal elements
        
        pos_term (SKiDl net or pin): positive terminal of the nonideal voltage source
            to connect to the rest of the circuit 
        
        neg_term (SKiDl net or pin): negative terminal of the nonideal voltage source
            to connect to the rest of the circuit 
        
        starting_V (float; 1; Amps): the initial DC current to set the internal ideal
            the current source in this package to
        
        starting_R (float; 50; Ohm): the initial resistance to set the internal
            parral resistance to the ideal current source in this subcircuit to
        
        return_internls (bool; False): If True return out the internal Voltage Source,
            and Resistance objects in this package
    
    Returns:
        Returns it's self a SKiDl part element object and if `return_internls`
        is True will return the internal voltage and resistance objects in that order 
    """
    cs=I(ref=f'I_{global_ref}', dc_value=starting_I)
    rp=R(ref=f'R_{global_ref}', value=starting_R)
    
    cs['p', 'n'] | rp[2, 1]
    rp[1, 2]+=pos_term, neg_term
    
    if return_internls:
        return cs, rp

#chapter 1 section 4 easy_tf class
# class to make using the SPICE .tf function so much less of a pain
# has methods to all four case usages

class easy_tf():
    """
    Wraper tool for ngspice's to make performing trnafer function (.tf) anyslis
    a lot less of a headeach
    """
    def __init__(self, circ_netlist_obj):
        """
        Intialztoin method
        
        Args:
            circ_netlist_obj (pspice.Spice.Netlist.Circuit): the Netlist circuit produced 
                from SKiDl's `generate_netlist()`
        
        Returns:
            creates a simation from the netlist in `self.sim`
        
        TODO: add kwargs to pass to sim creation

        """
        
        self.circ_netlist_obj=circ_netlist_obj
        #create the sim
        self.sim=self.circ_netlist_obj.simulator()
    
    
    def dc_voltage_gain(self, input_voltage_source, output_node_pos=None , output_node_neg='0', output_voltag_source=None):
        """
        method to perfom .tf to measur an input voltage souce to either a node pair
        wich is eqivlint to the SPICE statoemt
        ```
        .tf(V(<node_1>, <node_2>) , <voltage_souce_in>)
        ```
        where `<node_2>` can be ommited where then seconed node is always the ground node
        or a seconed output voltage souce can be used in wich case this meothd is then eqivlitnt to 
        ```
        .tf(<voltage_souce_out> , <voltage_souce_in>)
        ```
        it will retun the input resitince looking from the input souce to the output measure, the ouput resitince looking
        from the ouput meassue to the input voltage souce, and finaly the voltage gain of the input voltage souce to ouput meassure
        
        Args:
            input_voltage_source (SKiDl Linear Voltage Souce obj): the input voltage souce given as a SKiDl linear depened voltage souce like
                object. 
            
            output_node_pos (SKiDL net obj or SPICE net name str; None): If using a voltage souce as the ouput measurment this may be left as None
                else must be filled in with eather a SPICE net name as a string or a SKiDl net complint object where the SPICE
                net name can be exstracted from
            
            output_node_neg (SKiDL net obj or SPICE net name str; '0'): If using a voltage souce as the ouput measurment this may be left as '0' else
                if using measuing from a node pair must be eather a SPICE net name as a string or a SKiDl net complint object where the SPICE
                net name can be exstracted from. If left '0' the measured seconed node of the node pair will be the ground net '0'
            
            output_voltag_source (SKiDl Linear Voltage Souce obj; None): If meauing the ouput from a node pair (`output_node_pos` & `output_node_neg`) this may be 
                left as None. Else when measuing ouput voltge from a voltage souce provined a SKiDl linear depened voltage souce like
                object.
        
        Returns:
            results from this usge of .tf to measre the voltage to voltage gain will be stored in `self.vg_results`
            
        """
        #do assertion check that the input is a voltage source and 
        assert input_voltage_source.ref_prefix=='V', 'voltage gain must be sourced from a voltage source'
        input_src=get_skidl_spice_ref(input_voltage_source)
        
        
        
        #deal with source output
        if output_voltag_source!=None:

            if (output_node_pos!=None) & (output_node_neg!='0'):
                warnings.warn( 'If output measurment is from a `V` source, nodes node output location will be ignored')
            #creat the output statment
            output_src=f'V({get_skidl_spice_ref(output_voltag_source)})'
            #run the .tf
            self.tf_res=self.sim.tf(output_src, input_src)
            
            
            output_src=output_src[2:-1].replace(", ", "-")
            #create storage
            self.vg_results=pd.DataFrame(columns=['loc', 'value', 'units'])
            
            #store results
            for k, v in self.tf_res.nodes.items():
                if 'Transfer_function' in k:
                    self.vg_results.at['DC_Vgain']=[f'{output_src}/{input_src}', v.as_ndarray()[0], '[V/V]']
                elif 'Input_impedance' in k:
                    self.vg_results.at['Input_resistance']=[input_src, v.as_ndarray()[0], '[Ohm]']
                elif 'output_impedance' in k:
                    self.vg_results.at['Output_resistance']=[output_src, v.as_ndarray()[0], '[Ohm]']
                else:
                    warnings.warn(f'unexspected addintal outputs `{k, v}`')
            
        #deal with node ouput
        else:
            #check the pos node
            assert output_node_pos!=None, 'output node must be specfied when not using a `V` source'
            #TODO need to make this assertion stronger to the cirucit under test
            assert repr(type(output_node_pos))=="<class 'skidl.Net.Net'>" or isinstance(output_node_pos, str), '`output_node_pos` must be a SKiDl net or a node string' 
            if repr(type(output_node_pos))=="<class 'skidl.Net.Net'>":
                pos_node=node(output_node_pos)
            else:
                pos_node=output_node_pos
            
            #check the neg node
            #TODO need to make this assertion stronger to the cirucit under test
            assert repr(type(output_node_neg))=="<class 'skidl.Net.Net'>" or isinstance(output_node_neg, str), '`output_node_neg` must be a SKiDl net or a node string' 
            if repr(type(output_node_neg))=="<class 'skidl.Net.Net'>":
                neg_node=node(output_node_neg)
            else:
                neg_node=output_node_neg
                
            #formulate the output statment
            output_src=f'V({pos_node}, {neg_node})'
            
            #run the tf
            self.tf_res=self.sim.tf(output_src, input_src)
            
            self.vg_results=pd.DataFrame(columns=['loc', 'value', 'units'])
            
            output_src='('+output_src[2:-1].replace(", ", "-")+')'
            #store results
            for k, v in self.tf_res.nodes.items():
                if 'Transfer_function' in k:
                    self.vg_results.at['DC_Vgain']=[f'{output_src}/{input_src}', v.as_ndarray()[0], '[V/V]']
                elif 'Input_impedance' in k:
                    self.vg_results.at['Input_resistance']=[input_src, v.as_ndarray()[0], '[Ohm]']
                elif 'output_impedance' in k:
                    self.vg_results.at['Output_resistance']=[output_src, v.as_ndarray()[0], '[Ohm]']
                else:
                    warnings.warn(f'unexspected addintal outputs `{k, v}`')

        
            
        

        
    
    def dc_current_gain(self, input_current_source, output_ammeter):
        """
        method to perfom .tf to measur a input current souce to ouput ammeter
        
        wich is eqivlint to the SPICE statoemt
        ```
        .tf(I(<output_ammeter>) , <input_current_source>)
        ```
        it will retun the input resitince looking from the input souce to the output measure, the ouput resitince looking
        from the ouput meassue to the input current souce, and finaly the current gain of the input current souce to ouput meassure
        
        Args:
            input_current_source (SKiDl Linear Current Souce obj): the input current souce given as a SKiDl linear depened current souce like
                object. 
           
            output_ammeter (SKiDl Linear Voltage Souce obj; None): This may be a already esisitng linear voltage souce or a 0V linear voltage souce 
                SPICE ammeater add speficly to meassure the current
        
        Returns:
            results from this usge of .tf to measre the current to current gain will be stored in `self.ig_results`
            
        """
        assert input_current_source.ref_prefix=='I', 'current gain must be sourced from a current source'
        assert output_ammeter.ref_prefix=='V', 'current gain must be meassured from a voltage source'
        
        input_src=get_skidl_spice_ref(input_current_source)
        output_src=get_skidl_spice_ref(output_ammeter)
        
        self.tf_res=self.sim.tf(f'I({output_src})', input_src)
        
        self.ig_results=pd.DataFrame(columns=['loc', 'value', 'units'])
            
        for k, v in self.tf_res.nodes.items():
            if 'Transfer_function' in k:
                self.ig_results.at['DC_Igain']=[f'{output_src}/{input_src}', v.as_ndarray()[0], '[A/A]']
            elif 'Input_impedance' in k:
                self.ig_results.at['Input_resistance']=[input_src, v.as_ndarray()[0], '[Ohm]']
            elif 'Output_impedance' in k:
                self.ig_results.at['Output_resistance']=[output_src, v.as_ndarray()[0], '[Ohm]']
            else:
                warnings.warn(f'unexspected addintal outputs `{k, v}`')


    
    def dc_transresistance(self, input_current_source, output_node_pos=None , output_node_neg='0', output_voltag_source=None):
        
        """
        method to perfom .tf to fine the output voltge to input current.
        via measuring from an input current souce to either a node pair
        wich is eqivlint to the SPICE statoemt
        ```
        .tf(V(<node_1>, <node_2>) , <current_souce_in>)
        ```
        where `<node_2>` can be ommited where then seconed node is always the ground node
        or a seconed output voltage souce can be used in wich case this meothd is then eqivlitnt to 
        ```
        .tf(<voltage_souce_out> , <current_souce_in>)
        ```
        it will retun the input resitince looking from the input souce to the output measure, the ouput resitince looking
        from the ouput meassue to the input voltage souce, and finaly the transriseiince [A/V] gain of the input current souce to ouput meassure
        
        Args:
            input_voltage_source (SKiDl Linear Current Souce obj): the input current souce given as a SKiDl linear depened current souce like
                object. 
            
            output_node_pos (SKiDL net obj or SPICE net name str; None): If using a voltage souce as the ouput measurment this may be left as None
                else must be filled in with eather a SPICE net name as a string or a SKiDl net complint object where the SPICE
                net name can be exstracted from
            
            output_node_neg (SKiDL net obj or SPICE net name str; '0'): If using a voltage souce as the ouput measurment this may be left as '0' else
                if using measuing from a node pair must be eather a SPICE net name as a string or a SKiDl net complint object where the SPICE
                net name can be exstracted from. If left '0' the measured seconed node of the node pair will be the ground net '0'
            
            output_voltag_source (SKiDl Linear Voltage Souce obj; None): If meauing the ouput from a node pair (`output_node_pos` & `output_node_neg`) this may be 
                left as None. Else when measuing ouput voltge from a voltage souce provined a SKiDl linear depened voltage souce like
                object.
        
        Returns:
            results from this usge of .tf to measre the current to voltage gain will be stored in `self.tr_results`
            
        """
        assert input_current_source.ref_prefix=='I', 'transresistance gain must be sourced from a current source'
        input_src=get_skidl_spice_ref(input_current_source)
        
        #deal with node or source output
        if output_voltag_source!=None:

            if (output_node_pos!=None) & (output_node_neg!='0'):
                warnings.warn( 'If output measurment is from a `V` source, nodes node output location will be ignored')

            output_src=f'V({get_skidl_spice_ref(output_voltag_source)})'
            self.tf_res=self.sim.tf(output_src, input_src)
            
            
            output_src=output_src[2:-1].replace(", ", "-")
            
            self.tr_results=pd.DataFrame(columns=['loc', 'value', 'units'])
            
            for k, v in self.tf_res.nodes.items():
                if 'Transfer_function' in k:
                    self.tr_results.at['DC_Transresistance']=[f'{output_src}/{input_src}', v.as_ndarray()[0], '[V/A]']
                elif 'Input_impedance' in k:
                    self.tr_results.at['Input_resistance']=[input_src, v.as_ndarray()[0], '[Ohm]']
                elif 'output_impedance' in k:
                    self.tr_results.at['Output_resistance']=[output_src, v.as_ndarray()[0], '[Ohm]']
                else:
                    warnings.warn(f'unexspected addintal outputs `{k, v}`')
            
        
        else:
            assert output_node_pos!=None, 'output node must be specfied when not using a `V` source'
            #TODO need to make this assertion stronger to the cirucit under test
            assert repr(type(output_node_pos))=="<class 'skidl.Net.Net'>" or isinstance(output_node_pos, str), '`output_node_pos` must be a SKiDl net or a node string' 
            if repr(type(output_node_pos))=="<class 'skidl.Net.Net'>":
                pos_node=node(output_node_pos)
            else:
                pos_node=output_node_pos
            
            #TODO need to make this assertion stronger to the cirucit under test
            assert repr(type(output_node_neg))=="<class 'skidl.Net.Net'>" or isinstance(output_node_neg, str), '`output_node_neg` must be a SKiDl net or a node string' 
            if repr(type(output_node_neg))=="<class 'skidl.Net.Net'>":
                neg_node=node(output_node_neg)
            else:
                neg_node=output_node_neg
                
                
            output_src=f'V({pos_node}, {neg_node})'
            self.output_src=output_src
            
            self.tf_res=self.sim.tf(output_src, input_src)
            
            self.tr_results=pd.DataFrame(columns=['loc', 'value', 'units'])
            
            output_src='('+output_src[2:-1].replace(", ", "-")+')'
        
            for k, v in self.tf_res.nodes.items():
                if 'Transfer_function' in k:
                    self.tr_results.at['DC_Transresistance']=[f'{output_src}/{input_src}', v.as_ndarray()[0], '[V/A]']
                elif 'Input_impedance' in k:
                    self.tr_results.at['Input_resistance']=[input_src, v.as_ndarray()[0], '[Ohm]']
                elif 'output_impedance' in k:
                    self.tr_results.at['Output_resistance']=[output_src, v.as_ndarray()[0], '[Ohm]']
                else:
                    warnings.warn(f'unexspected addintal outputs `{k, v}`')
    
    def dc_transconductance(self, input_voltage_source, output_ammeter):
        """
        method to perfom .tf to measur a input voltage souce to an ouput ammeater
        
        wich is eqivlint to the SPICE statoemt
        ```
        .tf(I(<output_ammeter>) , <input_voltage_source>)
        ```
        
        
        it will retun the input resitince looking from the input souce to the output measure, the ouput resitince looking
        from the ouput meassue to the input current souce, and finaly the transconductance gain of the input current souce to ouput meassure
        
        Args:
            input_voltage_source (SKiDl Linear Voltage Souce obj): the input Voltage souce given as a SKiDl linear depened Voltage souce like
                object. 
           
            output_ammeter (SKiDl Linear Voltage Souce obj; None): This may be a already esisitng linear voltage souce or a 0V linear voltage souce 
                SPICE ammeater add speficly to meassure the current
        
        Returns:
            results from this usge of .tf to measre the current to current gain will be stored in `self.tc_results`
            
        """
        assert input_voltage_source.ref_prefix=='V', 'transconductance gain must be sourced from a voltage source'
        
        input_src=get_skidl_spice_ref(input_voltage_source)

        assert output_ammeter.ref_prefix=='V', 'transconductance gain must be meassured from a voltage source'

        output_src=f'I({get_skidl_spice_ref(output_ammeter)})'
            
        self.tf_res=self.sim.tf(output_src, input_src)


        output_src=output_src[2:-1].replace(", ", "-")

        self.tc_results=pd.DataFrame(columns=['loc', 'value', 'units'])

        for k, v in self.tf_res.nodes.items():
            if 'Transfer_function' in k:
                self.tc_results.at['DC_Transconductance']=[f'{output_src}/{input_src}', v.as_ndarray()[0], '[A/V]']
            elif 'Input_impedance' in k:
                self.tc_results.at['Input_resistance']=[input_src, v.as_ndarray()[0], '[Ohm]']
            elif 'Output_impedance' in k:
                self.tc_results.at['Output_resistance']=[output_src, v.as_ndarray()[0], '[Ohm]']
            else:
                warnings.warn(f'unexspected addintal outputs `{k, v}`')
    
        
#chapter 1 section 5 Thevenin class
# class that is a automated testsuite for finding the 
# DC Thevenin voltage and resistince of a port

class Thevenin(dc_ease):
    """
    Tool for finding the DC Thevenin equivalent voltage and resistance of a one port
    linear circuit with an open test port. The DUT must contain all linear elements.
    This class inheritance from the class `dc_ease` which is a tool to wrap and conduct .dc
    SPICE simulation with SkiDl and pyspice
    """
    
    def __init__(self, port_pos, port_neg):
        """
        A new instantiation unique to this class but utilizes `dc_ease` internally.
        Will add a 1A current source to the open port and will then generate the netlist
        to simulate against within this class
        
        Args:
            port_pos (SkiDl net): A SKiDl net that makes up the positive side of 
                the DUT's open port to test
            
            port_neg (SkiDl net): A SKiDl net that makes up the negative side of 
                the DUT's open port to test
        
        Returns:
            adds a 1A test source: `self.ithev` and creates that netlist to test:
            `self.circ_netlist_obj`
        
        TODO:
            -add assertions so that only a SKiDl net obj can only be passed in
            -add an assertion to see if the port is truly open
        """
        
        #add a current source to get Thevenin at port
        self.ithev=I(ref='Ithev', dc_value=1@u_A)
        self.ithev['n', 'p']+=port_pos, port_neg
        self.ithev_ref=self.ithev.ref
        
       
        #call generate netlist to create circ like we would if using `dc_ease`
        #by its self and simultaneously pass it into while invoking `dc_ease`'s
        # own instatation method
        super().__init__(generate_netlist())
        #print out the resulting circuit for debugging
        print('circuit and Thevenin finding Isource')
        print(self.circ_netlist_obj)
    
    
    def find_thev(self):
        """
        method to conduct the dc sweep of `self.ithev` automatically and find the 
        resulting Thevenin voltage and resistance
        
        Args: NONE
        
        Returns:
            see `dc_ease`'s `record_dc_nodebranch` for a full list of returns. But specifically will
            return `self.thevenin_sweep_data` which is the pandas' object with the
            returned sweep data with only the necessary data for finding the Vth and Rth.
            `self.rthev` and `self.vthev` which are Vth and Rth as floats respectively and
            `self.thev_results` which is a panda dataframe presenting Vth and Rth in a table
        
        TODO:
            
        """
        #set the sweep table
        self.sweep_DF.at[self.ithev_ref]=[0, 1, 0.1]
        self.clean_table()
        
        #do the sweep
        self.do_dc_sim(self.ithev_ref)
        #get the results 
        self.record_dc_nodebranch(self.ithev_ref)
        
        #reduce the data
        self.thevenin_sweep_data=self.dc_resultsNB_DF[node(self.ithev['n'])+'_[V]']
        
        #perform the 1d polyfit
        self.rthev, self.vthev=np.polyfit(self.thevenin_sweep_data.index, self.thevenin_sweep_data.values, 1)
        
        #make a pandas dataframe table, because it's the nice thing to do
        self.thev_results=pd.DataFrame(index=['Rthev', 'Vthev'], columns=['Values', 'Units'])
        self.thev_results['Units']=['[Ohm]', '[V]']
        self.thev_results['Values']=[self.rthev, self.vthev]
        
        
    def display_thev(self):
        """
        Auxiliary function to plot our Thevenin findings, useful for presentions
        and debugging
        
        Args: None
        
        Returns:
            Creates a plot and with a table to the side with the Thevenin equivalent finding
            with a table to the side summering our finding
        
        TODO:
            - add a check to make sure self.find_thev has been done
            
        """
        #create the plot "canvas"
        fig, axis=plt.subplots(nrows=1, ncols=1)
        
        #make a line plot on canvas of sweep data 
        self.thevenin_sweep_data.plot(grid=True, 
        xlabel=self.ithev_ref+'_[A]', ylabel=node(self.ithev['n'])+'_[V]', 
                                     ax=axis)
        
        #create "annotation" to hold self.thev_results in and to display next to the plot
        Anotation=['Theven Port Results:\n']
        for row in self.thev_results.itertuples():
            Anotation.append('   '+f'{row.Index}: {row.Values:.3f} {row.Units}' +'\n')
        
        Anotation=''.join(Anotation)

        fig.text(1.1,.4, Anotation, fontsize=14, transform=fig.transFigure)
        fig.suptitle('Thevin current sweep results')

#chapter 1 section 5 Norton class
# class that is a automated testsuite for finding the 
# DC Norton current and resistince of a port

class Norton(dc_ease):
    """
    Tool for finding the DC Norton equivalent current and resistance of a one port
    linear circuit with an open test port. The DUT must contain all linear elements.
    This class inheritance from the class `dc_ease` which is a tool to wrap and conduct .dc
    SPICE simulation with SkiDl and pyspice
    """
    
    def __init__(self, port_pos, port_neg):
        """
        A new instantiation unique to this class but utilizes `dc_ease` internally.
        Will add a 1V voltage source to the open port and will then generate the netlist
        to simulate against within this class
        
        Args:
            port_pos (SkiDl net): A SKiDl net that makes up the positive side of 
                the DUT's open port to test
            
            port_neg (SkiDl net): A SKiDl net that makes up the negative side of 
                the DUT's open port to test
        
        Returns:
            adds a 1V test source: `self.vnort` and creates that netlist to test:
            `self.circ_netlist_obj`
        
        TODO:
            -add assertions to so that only a SKiDl net obj can only be passed in
            -add an assertion to see if the port is truly open
        """
        
        #add voltage source to get Norton at port
        self.vnort=V(ref='Vnort', dc_value=1@u_V)
        self.vnort['p', 'n']+=port_pos, port_neg
        self.vnort_ref=self.vnort.ref
        
       
        #call generate netlist to create circ like we would if using `dc_ease`
        #by its self and simultaneously pass it into while invoking `dc_ease`'s
        # own instatation method
        super().__init__(generate_netlist())
        #print out the resulting circuit for debugging
        print('circuit and Thevenin finding Isource')
        print(self.circ_netlist_obj)
    
    
    def find_nort(self):
        """
        method to conduct the dc sweep of `self.ithev` automatically and find the 
        resulting Norton current and resistance
        
        Args: NONE
        
        Returns:
            see `dc_ease`'s `record_dc_internals` for a full list of returns. But specifically will
            return `self.norton_sweep_data` which is the pandas' object with the
            returned sweep data with only the necessary data for finding the Vth and Rth.
            `self.rnort` (`self.gnort`) and `self.inort` which are In and Rn as floats respectively and
            `self.nort_results` which is a pandas dataframe presenting In and Rn in a table
        
        TODO:
            
        """
        #set the sweep table
        self.sweep_DF.at[self.vnort_ref]=[0, 1, 0.1]
        self.clean_table()
        
        #do the sweep
        self.do_dc_intsim(self.vnort_ref)
        #get the results 
        self.record_dc_internals(self.vnort_ref)
        
        #reduce the data
        self.norton_sweep_data=self.dc_resultsINT_DF[self.vnort_ref+'[i]_[A]']
        
        #perform the 1d polyfit
        self.gnort, self.inort=np.polyfit(self.norton_sweep_data.index, self.norton_sweep_data.values, 1)
        #flip from conductance to resistance
        self.rnort=1/self.gnort
        
        #make a pandas dataframe table, because it's the nice thing to do
        self.nort_results=pd.DataFrame(index=['Rnort', 'Inort'], columns=['Values', 'Units'])
        self.nort_results['Units']=['[Ohm]', '[A]']
        self.nort_results['Values']=[self.rnort, self.inort]
        
        
    def display_nort(self):
        """
        Auxiliary function to plot our Norton findings, useful for presentions
        and debugging
        
        Args: None
        
        Returns:
            Creates a plot and with a table to the side with the Norton equivalent finding
            with a table to the side summering our finding
        
        TODO:
            - add a check to make sure self.find_nort  has been done
            
        """
        #create the plot "canvas"
        fig, axis=plt.subplots(nrows=1, ncols=1)
        
        #make a line plot on canvas of sweep data 
        self.norton_sweep_data.plot(grid=True, 
        xlabel=self.vnort_ref+'_[V]', ylabel=self.vnort_ref+'[i]_[A]', 
                                     ax=axis)
        
        #create "annotation" to hold self.nort_results in and to display next to the plot
        Anotation=['Norton Port Results:\n']
        for row in self.nort_results.itertuples():
            Anotation.append('   '+f'{row.Index}: {row.Values:.3f} {row.Units}' +'\n')
        
        Anotation=''.join(Anotation)

        fig.text(1.1,.4, Anotation, fontsize=14, transform=fig.transFigure)
        fig.suptitle('Norten voltage sweep results')

#chapter 1 section 6 findIntersection function
#Assist function to find the intersection of two functions

#from https://glowingpython.blogspot.com/2011/05/hot-to-find-intersection-of-two.html 
#load fslove from scipy's optimize module
from scipy.optimize import fsolve

#helper function to find the intersection of two functions with an initial guess
def findIntersection(fun1,fun2,x0):
    """
    Aid function to find the intersection point of two curves
    from: https://glowingpython.blogspot.com/2011/05/hot-to-find-intersection-of-two.html 
    
    Args:
        func1(function or class): the first function whose curve is 
            used to find the intersection of the two curves
        
        func2(function or class): the second function whose curve is 
            used to find the intersection of the two curves
        
        x0 (float); initial guess of the intersection of the two functions
    
    Returns:
        Returns array of float that are the intersections of the two functions, 
        this is not very robust and thus one should read `fsolve`'s documentation 
        for caveats  of usage
    """
    return fsolve(lambda x : fun1(x) - fun2(x),x0)

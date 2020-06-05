#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Wed Jan 31 02:00:37 2018

@author: cbkeller
"""
#def configure_isobaric(perplexdir, scratchdir, composition, elements = ['SIO2','TIO2','AL2O3','FEO','MGO','CAO','NA2O','K2O','H2O'], index = 1, P = 10000, T_range = [500+273.15, 1500+273.15], dataset = 'hp11ver.dat', solution_phases = 'O(HP)\nOpx(HP)\nOmph(GHP)\nGt(HP)\noAmph(DP)\ncAmph(DP)\nT\nB\nChl(HP)\nBio(TCC)\nMica(CF)\nCtd(HP)\nIlHm(A)\nSp(HP)\nSapp(HP)\nSt(HP)\nfeldspar_B\nDo(HP)\nF\n', excludes = 'ts\nparg\ngl\nged\nfanth\ng\n'):

# Import some useful packages
import os # os.system lets us access the command line
import numpy as np # For np.array
import pandas as pd # Pandas, for importing PerpleX text file output as data frames

def configure(meltspath, scratchdir, batchstring,T_range, P_range,
              dT=-10, dP=0, index=1,mass_retained=0.001,min_fractionation=0.005,
              version='pMELTS',mode='isobaric',fo2path='FMQ',elements=[],composition=[],
              trace_elements=[],trace_composition=[],phases_fractionated=[],phases_suppressed = [],
              fractionate_solids=True,fractionate_water=True,continuous_fractionation_melting=False,
              output_celsius=True,save_output=True,do_trace_elements=False,water_is_trace=False,
              **kwargs):
    """

    Parameters
    ----------
    meltspath: str
        Path to MELTS or pMELTS
    scratchdir: str
        A scratch directory path
    composition: np.ndarray
        An array of compositions corresponding to either wt% or mol%. Unsure.
    elements: tuple
        An array of elements to consider
    batchstring:
        not sure.
    T_range: float or tuple(float)
        A temperature or temperature range to apply. Unsure if in celcius or kelvin
    P_range: float
        A pressure or pressure range to apply in MPa
    dT: float
        temperature increments to explore (C? K?)
    dP: float
        pressure increments to explore (MPa?)
    min_fractionation: float
        threshold at which melt is removed
    index: int
        index?
    version: str
        unsure.
    mode:
        define mode as
        isothermal, isobaric, isentropic, isenthalpic, isochoric, geothermal or PTPath
    fo2path: str
        set oxygen fugacity path/buffer. can be any of:
        'FMQ','IW','COH','FMQ','NNO','HM','None'
    fractionate_solids: bool
        fractionate solids?
    fractionate_water: bool
        if should fractionate water
    continuous_fractionation_melting: bool
        if should continuously fractionate during melting
    do_trace_elements: bool
        if should fractionate trace elements
    water_is_trace: bool
        if should treat water as a trace element
    output_celsius: bool
        if perple_x should output celsius
    save_output: bool
        if perple_x should save output
    mass_retained: float
        mass retained during fractionation. default is 1e-3
    min_fractionation: float
        minimal mass needed to fractionate
    phases_fractionated: list
        a list of phases to be fractionated. ex: ['olivine','spinel']
    phases_suppressed: list
        a list of phases to suppress. ex: ['leucite']
    trace_elements: list
        a list of trace elements
    trace_composition: list
        a list of trace element starting composition
    Returns
    -------

    """

    celsius_out   ='' if output_celsius     else celsius_out='!'
    save_all      ='' if save_output        else save_all='!'
    frack_solids  ='' if fractionate_solids else frack_solids = '!'
    frack_water   ='' if fractionate_water  else frack_water='!'
    cont_melt     ='' if continuous_fractionation_melting else cont_melt ='!'
    water_trace   ='' if water_is_trace     else water_trace = '!'
    do_trace      ='' if do_trace_elements  else do_trace='!'

    Pmax=90000
    Pmin=2
    Tmax=3000
    Tmin=450
    # Simulation number (for folder, etc)
    
    ########################## end Default Settings ############################
    
    # Guess if intention is for calculation to end at Tf or Pf as a min or max
    if T_range[1]<T_range[0]:
        Tmin=T_range[1]
    if T_range[1]>T_range[0]:
        Tmax=T_range[1]
    if P_range[1]<P_range[0]:
        Pmin=P_range[1]
    if P_range[1]>P_range[0]:
        Pmax=P_range[1]
    
    # Normalize starting composition
    composition = np.array(composition)/sum(composition)*100.0
    
    # output prefixectory name
    prefix = scratchdir + os.sep + '{index}'
    os.system(f'rm -rf {prefix}; mkdir -p {prefix}') # Ensure directory is empty
    
    # Make .melts file containing the starting composition you want to run
    # simulations on
    with open(prefix+os.sep + 'sc.melts','w') as melts:
        for element,comp in zip(elements,composition):
            melts.write(f'Initial Composition: {element} {comp:4f}\n')
        for trace, comp in zip(trace_elements,trace_composition):
            melts.write(f'Initial Trace: {trace} {comp:4f}\n')
        
        melts.write(f'Initial Temperature: {T_range[0]:2f}\n'+
                 f'Initial Pressure: {P_range[0]:2f}\n'+
                 f'log fo2 Path: {fo2path}\n')
    
        for phase in phases_fractionated:
            melts.write(f'Fractionate: {phase}\n')
        for phase in phases_suppressed:
            melts.write(f'Suppress: {phase}\n')
    
    # Make melts_env file to specify type of MELTS calculation
    with open(f'{prefix}'+os.sep + 'melts_env.txt','w') as new_file:
        new_file.write('! *************************************\n!  Python-generated environment file\n! *************************************\n\n' +
        '! this variable chooses MELTS or pMELTS; for low-pressure use MELTS\n' +
        f'ALPHAMELTS_VERSION		{version}\n\n' +
        '! do not use this unless fO2 anomalies at the solidus are a problem\n' +
        '!ALPHAMELTS_ALTERNATIVE_FO2	true\n\n' +
        '! use this if you want to buffer fO2 for isentropic, isenthalpic or isochoric mode\n! e.g. if you are doing isenthalpic AFC\n' +
        '!ALPHAMELTS_IMPOSE_FO2		true\n\n' +
        '! use if you want assimilation and fractional crystallization (AFC)\n' +
        '!ALPHAMELTS_ASSIMILATE		true\n\n' +
        '! isothermal, isobaric, isentropic, isenthalpic, isochoric, geothermal or PTPath\n' + \
        f'ALPHAMELTS_MODE			{mode}\n'  + \
        '!ALPHAMELTS_PTPATH_FILE		ptpath.txt\n\n' + \
        f'! need to set DELTAP for polybaric paths; DELTAT for isobaric paths\nALPHAMELTS_DELTAP	{dP:0f}\n' + \
        f'ALPHAMELTS_DELTAT	{dT:1f}\n'  + \
        f'ALPHAMELTS_MAXP		{Pmax:0f}\n' + \
        f'ALPHAMELTS_MINP		{Pmin:0f}\n'  + \
        f'ALPHAMELTS_MAXT		{Tmax:0f}\n'  + \
        f'ALPHAMELTS_MINT		{Tmin:0f}\n\n'  + \
        '! this one turns on fractional crystallization for all solids\n! use Fractionate: in the melts file instead for selective fractionation\n' + \
        f'{frack_solids}ALPHAMELTS_FRACTIONATE_SOLIDS	true\n'  + \
        f'{frack_solids}ALPHAMELTS_MASSIN		{mass_retained:0.f}\n\n'  + \
        '! free water is unlikely but can be extracted\n' + \
        f'{frack_water}ALPHAMELTS_FRACTIONATE_WATER	true\n' + \
        f'{frack_water}ALPHAMELTS_MINW			0.005\n\n'  + \
        '! the next one gives an output file that is always updated, even for single calculations\n' + \
        f'{save_all}ALPHAMELTS_SAVE_ALL		true\n'+ \
        '!ALPHAMELTS_SKIP_FAILURE		true\n\n' + \
        '! this option converts the output temperature to celcius, like the input\n' + \
        f'{celsius_out}ALPHAMELTS_CELSIUS_OUTPUT	true\n\n'  + \
        '! the next two turn on and off fractional melting\n' + \
        f'{cont_melt}ALPHAMELTS_CONTINUOUS_MELTING	true\n'  + \
        f'{cont_melt}ALPHAMELTS_MINF			{min_fractionation:6.4f}\n'  + \
        f'{cont_melt}ALPHAMELTS_INTEGRATE_FILE	integrate.txt\n\n'  + \
        '! the next two options refer to the trace element engine\n' + \
        f'{do_trace}ALPHAMELTS_DO_TRACE		true\n'+ \
        f'{water_trace}ALPHAMELTS_DO_TRACE_H2O		true\n')

    # Make a batch file to run the above .melts file starting from the liquidus
    with open('{prefix}' + os.sep + 'batch.txt', 'w') as batch:
        batch.write(batchstring)
    
    # Run the command
    # Edit the following line(s to make sure you have a correct path to the 'run_alphamelts.command' perl script
    os.system('cd ' + prefix  + '; ' + meltspath + ' -f melts_env.txt -b batch.txt')
    return


# Get modal phase proportions, return as pandas DataFrame
def query(scratchdir, index=1):
    prefix = scratchdir + f'out{index}/'  # path to data files
    file = prefix + 'Phase_main_tbl.txt'
    return pd.read_csv(file)

# Get modal phase proportions, return as pandas DataFrame
def query_modes(scratchdir, index=1):
    prefix = scratchdir + f'out{index}/'  # path to data files
    n_header_lines = 1
    
    # Read results and return them if possible
    try:
        data = pd.read_csv(prefix + 'Phase_mass_tbl.txt', delim_whitespace=True, header=n_header_lines)
                # Ensure columns are numeric

        for c in data.columns:
            if data[c].dtype!='float64':
                data[c] = np.genfromtxt(data[c])
    except:
        data = 0
    return data


# Get liquid composition, return as pandas DataFrame
def query_liquid(scratchdir, index=1):
    prefix = scratchdir + f'out{index}/'  # path to data files
    n_header_lines = 1
    
    # Read results and return them if possible
    try:
        data = pd.read_csv(prefix + 'Liquid_comp_tbl.txt', delim_whitespace=True, header=n_header_lines)
        # Ensure columns are numeric
        for c in data.columns:
            if data[c].dtype!='float64':
                data[c] = np.genfromtxt(data[c])
    except:
        data = 0
    return data


# Read solid composition, return as pandas DataFrame
def query_solid(scratchdir, index=1):
    prefix = scratchdir + f'out{index}/' # path to data files
    n_header_lines = 1
    
    # Read results and return them if possible
    try:
        data = pd.read_csv(prefix + 'Solid_comp_tbl.txt', delim_whitespace=True, header=n_header_lines)
        # Ensure columns are numeric
        for c in data.columns:
            if data[c].dtype!='float64':
                data[c] = np.genfromtxt(data[c])
    except:
        data = 0
    return data

# Read system thermodynamic data, return as pandas DataFrame
def query_system(scratchdir, index=1):
    prefix = scratchdir + f'out{index}/' # path to data files
    n_header_lines = 1
    
    # Read results and return them if possible
    try:
        data = pd.read_csv(prefix + 'System_main_tbl.txt', delim_whitespace=True, header=n_header_lines)
        # Ensure columns are numeric
        for c in data.columns:
            if data[c].dtype!='float64':
                data[c] = np.genfromtxt(data[c])
    except:
        data = 0
    return data

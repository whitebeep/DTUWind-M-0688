#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 5 11:10:00 2023

@author: eduardantoligil
"""

#%% Import Module
def import_status():
    print('Module imported successfully: power/power_systems.py')
    
#%% Import
import numpy as np
from numpy.linalg import inv

#%% Main functions

def init_network_data():    
    global bus, branch, gen
    
    bus_info = ['Nr', 'Type', 'Pld', 'Qld', 'Gsh', 'Bsh', 'Vm', 'Va']
    branch_info = ['From', 'To', 'R', 'X', 'G', 'B', 'Tap', 'Phase', 'rating']
    gen_info = ['bus', 'P', 'Q', 'V', 'MVA', 'Pmax', 'Qmax', 'Pmin', 'Qmin']
    
    bus = dict.fromkeys(bus_info)
    branch = dict.fromkeys(branch_info)
    gen = dict.fromkeys(gen_info)
    
    return

def load_network_data(power_system, rm_line=None):
    
    # Load network
    ps_name = power_system(rm_line)
    
    # Assign bus indices
    global pq_index, pv_index, ref, pvpq_index
    pq_index, pv_index, ref, pvpq_index = bus_indices()
    
    # Initialze voltage and power data
    global V_0, S_bus
    V_0, S_bus = voltage_power_data()
    
    # Initialize branch-bus incidence matrix
    global C_ft
    C_ft = make_Cft()
    
    print('\nSuccessfully loaded: "{}"'.format(ps_name))
    
    return

def load_pf_data(pf_type='DC'):
    
    pf_types = 'DC', 'AC_NR', 'AC_FDLF'
    
    if pf_type in pf_types:
        
        if pf_type == 'DC':
            global Bdc_bus, Bdc_line # include Pbusinj, Pfinj in the future?
            Bdc_bus, Bdc_line = make_Bdc() #, Pbusinj, Pfinj
            
        elif pf_type == 'AC_NR':
            global Y_bus, Y_from, Y_to
            Y_bus, Y_from, Y_to = make_Ybus()
        
        elif pf_type == 'AC_FDLF':
            global B_p, B_pp
            B_p, B_pp = make_B()
        
        load_status = '\nSuccessfully loaded: {} PowerFlow data'.format(pf_type)
        
    else:
        load_status = '\nValueError: please specify a valid pf_type argument: "DC", "AC_NR" or "AC_FDLF"'
        
    print(load_status)
    
    return

def load_contingency_data():
    
    global ptdf, ptdf_noref, lodf
    ptdf = calculate_ptdf()
    ptdf_noref = calculate_ptdf(full=False)
    lodf = calculate_lodf()
    print('\nSuccessfully loaded: PTDF and LODF matrices')
    
    return

#%% Auxiliary functions

def bus_indices():
   
    buscode = bus['Type']
    pq_index = np.where(buscode == 1)[0] # Find indices for all PQ-busses
    pv_index = np.where(buscode == 2)[0] # Find indices for all PV-busses
    ref = np.where(buscode == 3)[0] # Find index for ref bus
    pvpq_index = np.concatenate((pv_index, pq_index))
    
    return pq_index, pv_index, ref, pvpq_index

def make_Cft():
    
    n_bus = len(bus['Nr'])
    n_br = len(branch['From'])
    
    # Define incidence matrix
    f, t = branch['From'], branch['To']
    Cft = np.zeros((n_br, n_bus))
    Cft[np.arange(n_br),f] = 1
    Cft[np.arange(n_br),t] = -1
    
    return Cft
    
def make_Bdc():
    
    # Line susceptance matrix
    tap = branch['Tap']
    b_line = 1/(branch['X'])/tap
    B_diag = np.diag(b_line)
    Bf = B_diag.dot(C_ft)
    
    # Bus susceptance matrix
    Bbus = (C_ft.T).dot(Bf)
    
    # # Line and bus injection
    # Pfinj = b_line * (-branch['Phase'] * np.pi/180)
    # Pbusinj = (C_ft.T).dot(Pfinj)

    return Bbus, Bf#, Pbusinj, Pfinj

def make_B():
   
    # Line admitance matrix
    Ys = 1/(1j*branch['X'])    # Yp = branch['G']+  1j*branch['B']
    Ys_diag = np.diag(Ys)
    Ybus = (C_ft.T).dot(Ys_diag.dot(C_ft))
    
    # Bp(P) and Bpp(Q) susceptance matrices
    Bp = np.imag(Ybus[pvpq_index, :][:, pvpq_index])
    Bpp = np.imag(Ybus[pq_index, :][:, pq_index])

    return Bp, Bpp

def make_Ybus():

    n_bus = len(bus['Nr'])
    n_br = len(branch['From'])
    
    # Calculate shunt, series and parallel admitances
    Ysh = bus['Gsh'] + 1j*bus['Bsh']
    Ys = 1/(branch['R'] + 1j*branch['X'])
    Yp = branch['G'] + 1j*branch['B']

    Yff = Ys + Yp/2
    Yft = - Ys
    Ytf = - Ys

    # Construct bus admitance matrix (Y_bus)
    Ybus = np.zeros((n_bus,n_bus),dtype=complex)
    f, t = branch['From'], branch['To']

    for i in range(n_bus):
        Ybus[i,i] += Ysh[i]
        
    for i in range(n_br):
        at = branch['Tap'][i]
        if at > 0:
            Ybus[f[i],f[i]] += Yff[i]/at+Yff[i]*(1-at)/(at**2)
            Ybus[f[i],t[i]] += Yft[i]/at
            Ybus[t[i],f[i]] += Ytf[i]/at
            Ybus[t[i],t[i]] += Yff[i]/at+Yff[i]*(at-1)/at
        else:
            Ybus[f[i],f[i]] += Yff[i]
            Ybus[f[i],t[i]] += Yft[i]
            Ybus[t[i],f[i]] += Ytf[i]
            Ybus[t[i],t[i]] += Yff[i]

    # Construct branch admitance matrices (Y_from and Y_to)
    Y_branch_diag = np.diag(Ys)
    Y_from = Y_branch_diag.dot(C_ft)
    Y_to = - Y_from

    return Ybus, Y_from, Y_to

def voltage_power_data():
    
    Sbus = bus['Pld'] + 1j*bus['Qld']
    V0 = bus['Vm'] * np.exp( 1j*np.deg2rad(bus['Va']) )

    return V0, Sbus

def calculate_ptdf(full=True):
    
    # Remove slack bus from Bdc_bus and Bdc_line 
    Bbus = Bdc_bus[pvpq_index, :][:, pvpq_index]
    Bline = Bdc_line[:, pvpq_index]
    
    # Calculate ptdf for the reduced system
    ptdf = Bline.dot(inv(Bbus))
    
    if full:
        # Include slack bus
        ptdf = np.insert(ptdf, ref, 0, axis=1)
        
    return ptdf

def calculate_lodf():
    
    nl = ptdf.shape[0]
    
    H = ptdf @ C_ft.T
    h = np.expand_dims(np.diag(H),1)
    lodf = H / (np.ones((nl, nl)) - np.ones((nl, 1)) @ h.T)
    np.fill_diagonal(lodf, -1)
    
    return lodf  
    
#%% List of Power Systems

def slides_3bus(rm_line=None):
    
    # Power system name ID
    psys_name = '3-bus system from slides'
    
    # Bus data
    n_bus = 3
    Sbase = 100
    
    bus_nr = np.arange(n_bus)
    bus_type = np.array([3, 2, 1])
    Pgen = np.array([0.0, 100.0, 0.0])/Sbase
    Qgen = np.array([0.0, 0.0, 0.0])/Sbase
    Pld = np.array([0.0, 0.0, 160.0])/Sbase
    Qld = np.array([0.0, 0.0, 10.0])/Sbase
    Gsh = np.zeros(n_bus)
    Bsh = np.zeros(n_bus)
    Vm = np.array([1.00, 1.00, 1.00])
    Va = np.zeros(n_bus)
    
    bus['Nr'] = bus_nr
    bus['Type'] = bus_type
    bus['Pld'] = Pgen-Pld
    bus['Qld'] = Qgen-Qld
    bus['Bsh'] = Bsh
    bus['Gsh'] = Gsh
    bus['Vm'] = Vm
    bus['Va'] = Va
    
    # Line data
    br_f = np.array([1, 1, 2]) -1
    br_t = np.array([2, 3, 3]) -1
    n_br = len(br_f)
    R_br = np.array([0.01, 0.02, 0.03])
    X_br = np.array([0.1, 0.25, 0.2])
    G_br = np.array([0.0, 0.0, 0.0])
    B_br = np.array([0.0, 0.0, 0.0])
    tap_br = np.ones(n_br)
    phase_br = np.zeros(n_br)
    rating_br = np.ones(n_br)

    if rm_line is not None: # Remove lines from the system
        br_f = np.delete(br_f, rm_line)
        br_t = np.delete(br_t, rm_line)
        R_br = np.delete(R_br, rm_line)
        X_br = np.delete(X_br, rm_line)
        G_br = np.delete(G_br, rm_line)
        B_br = np.delete(B_br, rm_line)
        tap_br = np.delete(tap_br, rm_line)
        phase_br = np.delete(phase_br, rm_line)
        rating_br = np.delete(rating_br, rm_line)
        
    branch['From'] = br_f
    branch['To'] = br_t
    branch['R'] = R_br
    branch['X'] = X_br
    branch['G'] = G_br
    branch['B'] = B_br
    branch['Tap'] = tap_br
    branch['Phase'] = phase_br
    branch['rating'] = rating_br # np.full(n_br,0.104/X_br) # TODO add propper line rating

    # Generator data
    genbus = np.flatnonzero(np.logical_or.reduce((Pgen != 0,Qgen !=0, bus_type==3))) # TODO multiple generators on the same bus
    n_gen = len(genbus)
    
    gen['bus'] = genbus
    gen['P'] = Pgen[genbus]
    gen['Q'] = Qgen[genbus]
    gen['V'] = Vm[genbus]
    gen['MVA'] = 100*np.ones(n_gen)
    gen['Pmax'] = np.ones(n_gen)
    gen['Qmax'] = np.ones(n_gen)
    gen['Pmin'] = np.zeros(n_gen)
    gen['Qmin'] = np.zeros(n_gen)

    return psys_name

def IEEE_5bus(rm_line=None):
    
    # Power system name ID
    psys_name = 'IEEE 5-bus system'
    
    # Bus data
    n_bus = 5
    Sbase = 100
    
    bus_nr = np.arange(n_bus)
    bus_type = np.array([3, 1, 1, 1, 1])
    Pgen = np.array([0.0, 40.0, 0.0, 0.0, 0.0])/Sbase
    Qgen = np.array([0.0, 30.0, 0.0, 0.0, 0.0])/Sbase
    Pld = np.array([0.0, 20.0, 45.0, 40.0, 60.0])/Sbase
    Qld = np.array([0.0, 10.0, 15.0, 5.0, 10.0])/Sbase
    Gsh = np.zeros(n_bus)
    Bsh = np.zeros(n_bus)
    Vm = np.array([1.06, 1.00, 1.00, 1.00, 1.00])
    Va = np.zeros(n_bus)
    
    bus['Nr'] = bus_nr
    bus['Type'] = bus_type
    bus['Pld'] = Pgen-Pld
    bus['Qld'] = Qgen-Qld
    bus['Bsh'] = Bsh
    bus['Gsh'] = Gsh
    bus['Vm'] = Vm
    bus['Va'] = Va
    
    # Line data
    br_f = np.array([1, 1, 2, 2, 2, 3, 4])-1
    br_t = np.array([2, 3, 3, 4, 5, 4, 5])-1
    n_br = len(br_f)
    R_br = np.array([0.02, 0.08, 0.06, 0.06, 0.04, 0.01, 0.08])
    X_br = np.array([0.06, 0.24, 0.25, 0.18, 0.12, 0.03, 0.24])
    G_br = np.array([0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00])
    B_br = np.array([0.03, 0.025, 0.02, 0.02, 0.015, 0.01, 0.025])*2
    tap_br = np.ones(n_br)
    phase_br = np.zeros(n_br)
    rating_br = np.ones(n_br)

    if rm_line is not None: # Remove lines from the system
        br_f = np.delete(br_f, rm_line)
        br_t = np.delete(br_t, rm_line)
        R_br = np.delete(R_br, rm_line)
        X_br = np.delete(X_br, rm_line)
        G_br = np.delete(G_br, rm_line)
        B_br = np.delete(B_br, rm_line)
        tap_br = np.delete(tap_br, rm_line)
        phase_br = np.delete(phase_br, rm_line)
        rating_br = np.delete(rating_br, rm_line)
        
    branch['From'] = br_f
    branch['To'] = br_t
    branch['R'] = R_br
    branch['X'] = X_br
    branch['G'] = G_br
    branch['B'] = B_br
    branch['Tap'] = tap_br
    branch['Phase'] = phase_br
    branch['rating'] = rating_br # np.full(n_br,0.104/X_br) # TODO add propper line rating

    # Generator data
    genbus = np.flatnonzero(np.logical_or.reduce((Pgen != 0,Qgen !=0, bus_type==3))) # TODO multiple generators on the same bus
    n_gen = len(genbus)
    
    gen['bus'] = genbus
    gen['P'] = Pgen[genbus]
    gen['Q'] = Qgen[genbus]
    gen['V'] = Vm[genbus]
    gen['MVA'] = 100*np.ones(n_gen)
    gen['Pmax'] = np.ones(n_gen)
    gen['Qmax'] = np.ones(n_gen)
    gen['Pmin'] = np.zeros(n_gen)
    gen['Qmin'] = np.zeros(n_gen)

    return psys_name

def IEEE_5bus_mp(rm_line=None):
    
    # Power system name ID
    psys_name = 'IEEE 5-bus MP system'
    
    # Bus data
    n_bus = 5
    Sbase = 100
    
    bus_nr = np.arange(n_bus)
    bus_type = np.array([2, 1, 2, 3, 2])
    Pgen = np.array([210.0, 0.0, 323.49, 0.0, 466.51])/Sbase
    Qgen = np.array([0.0, 0.0, 0.0, 0.0, 0.0])/Sbase
    Pld = np.array([0.0, 300.0, 300.0, 400.0, 0.0])/Sbase
    Qld = np.array([0.0, 98.61, 98.61, 131.47, 0.0])/Sbase
    Gsh = np.zeros(n_bus)
    Bsh = np.zeros(n_bus)
    Vm = np.array([1.00, 1.00, 1.00, 1.00, 1.00])
    Va = np.zeros(n_bus)
    
    bus['Nr'] = bus_nr
    bus['Type'] = bus_type
    bus['Pld'] = Pgen-Pld
    bus['Qld'] = Qgen-Qld
    bus['Bsh'] = Bsh
    bus['Gsh'] = Gsh
    bus['Vm'] = Vm
    bus['Va'] = Va
    
    # Line data
    br_f = np.array([1, 1, 1, 2, 3, 4])-1
    br_t = np.array([2, 4, 5, 3, 4, 5])-1
    n_br = len(br_f)
    R_br = np.array([0.00281, 0.00304, 0.00064, 0.00108, 0.00297, 0.00297])
    X_br = np.array([0.0281, 0.0304, 0.0064, 0.0108, 0.0297, 0.0297])
    G_br = np.array([0.00, 0.00, 0.00, 0.00, 0.00, 0.00])
    B_br = np.array([0.00712 ,0.00658, 0.03126, 0.01852, 0.00674, 0.00674])
    tap_br = np.ones(n_br)
    phase_br = np.zeros(n_br)
    rating_br = np.ones(n_br)

    if rm_line is not None: # Remove lines from the system
        br_f = np.delete(br_f, rm_line)
        br_t = np.delete(br_t, rm_line)
        R_br = np.delete(R_br, rm_line)
        X_br = np.delete(X_br, rm_line)
        G_br = np.delete(G_br, rm_line)
        B_br = np.delete(B_br, rm_line)
        tap_br = np.delete(tap_br, rm_line)
        phase_br = np.delete(phase_br, rm_line)
        rating_br = np.delete(rating_br, rm_line)
        
    branch['From'] = br_f
    branch['To'] = br_t
    branch['R'] = R_br
    branch['X'] = X_br
    branch['G'] = G_br
    branch['B'] = B_br
    branch['Tap'] = tap_br
    branch['Phase'] = phase_br
    branch['rating'] = rating_br # np.full(n_br,0.104/X_br) # TODO add propper line rating

    # Generator data
    genbus = np.flatnonzero(np.logical_or.reduce((Pgen != 0,Qgen !=0, bus_type==3))) # TODO multiple generators on the same bus
    n_gen = len(genbus)
    
    gen['bus'] = genbus
    gen['P'] = Pgen[genbus]
    gen['Q'] = Qgen[genbus]
    gen['V'] = Vm[genbus]
    gen['MVA'] = 100*np.ones(n_gen)
    gen['Pmax'] = np.ones(n_gen)
    gen['Qmax'] = np.ones(n_gen)
    gen['Pmin'] = np.zeros(n_gen)
    gen['Qmin'] = np.zeros(n_gen)

    return psys_name

def IEEE_5bus_mod(rm_line=None):
    
    # Power system name ID
    psys_name = 'IEEE 5-bus mod system'
    
    # Bus data
    n_bus = 5
    Sbase = 100
    
    bus_nr = np.arange(n_bus)
    bus_type = np.array([3, 1, 1, 1, 1])
    Pgen = np.array([0.0, 20.0, 0.0, 0.0, 0.0])/Sbase
    Qgen = np.array([0.0, 5.0, 0.0, 0.0, 0.0])/Sbase
    Pld = np.array([0.0, 10.0, 15.0, 10.0, 15.0])/Sbase
    Qld = np.array([0.0, 5.0, 7.5, 2.5, 5.0])/Sbase
    Gsh = np.zeros(n_bus)
    Bsh = np.zeros(n_bus)
    Vm = np.array([1.05, 1.00, 1.00, 1.00, 1.00])
    Va = np.zeros(n_bus)
    
    bus['Nr'] = bus_nr
    bus['Type'] = bus_type
    bus['Pld'] = Pgen-Pld
    bus['Qld'] = Qgen-Qld
    bus['Bsh'] = Bsh
    bus['Gsh'] = Gsh
    bus['Vm'] = Vm
    bus['Va'] = Va
    
    # Line data
    br_f = np.array([1, 1, 1, 1, 2, 3, 4])-1
    br_t = np.array([2, 3, 4, 5, 3, 4, 5])-1
    n_br = len(br_f)
    R_br = np.array([0.02, 0.08, 0.06, 0.06, 0.04, 0.01, 0.08])*0
    X_br = np.array([ 1.03424425,  0.52637178,  0.90166637,  1.17900178, 29.20287936, 15.59879165,  0.67128607])
    G_br = np.array([0.00, 0.00, 0.00, 0.00, 0.00, 0.00, 0.00])
    B_br = np.array([0.03, 0.025, 0.02, 0.02, 0.015, 0.01, 0.025])*2*0
    tap_br = np.ones(n_br)
    phase_br = np.zeros(n_br)
    rating_br = np.ones(n_br)

    if rm_line is not None: # Remove lines from the system
        br_f = np.delete(br_f, rm_line)
        br_t = np.delete(br_t, rm_line)
        R_br = np.delete(R_br, rm_line)
        X_br = np.delete(X_br, rm_line)
        G_br = np.delete(G_br, rm_line)
        B_br = np.delete(B_br, rm_line)
        tap_br = np.delete(tap_br, rm_line)
        phase_br = np.delete(phase_br, rm_line)
        rating_br = np.delete(rating_br, rm_line)
        
    branch['From'] = br_f
    branch['To'] = br_t
    branch['R'] = R_br
    branch['X'] = X_br
    branch['G'] = G_br
    branch['B'] = B_br
    branch['Tap'] = tap_br
    branch['Phase'] = phase_br
    branch['rating'] = rating_br # np.full(n_br,0.104/X_br) # TODO add propper line rating

    # Generator data
    genbus = np.flatnonzero(np.logical_or.reduce((Pgen != 0,Qgen !=0, bus_type==3))) # TODO multiple generators on the same bus
    n_gen = len(genbus)
    
    gen['bus'] = genbus
    gen['P'] = Pgen[genbus]
    gen['Q'] = Qgen[genbus]
    gen['V'] = Vm[genbus]
    gen['MVA'] = 100*np.ones(n_gen)
    gen['Pmax'] = np.ones(n_gen)
    gen['Qmax'] = np.ones(n_gen)
    gen['Pmin'] = np.zeros(n_gen)
    gen['Qmin'] = np.zeros(n_gen)

    return psys_name

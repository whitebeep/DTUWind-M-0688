#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Apr 14 20:46:12 2023

@author: eduardantoligil
"""

#%% Import Module
def import_status():
    print('Module imported successfully: postpro/postprocessing.py')
    
    
#%% Imports
import numpy as np
from datetime import datetime

# # Load ps.ptdf and ps.lodf matrices
# import power.power_systems as ps
# ps.load_contingency_data()
# PTDF, LODF = ps.ptdf, ps.lodf
    
#%% Post-processing functions

def postpro_hhl_dc(theta_hhl, bus, B, Bf, pvpq_index):
    
    # Assign theta from HHL DC Power Flow
    theta = bus['Va'].copy()
    theta[pvpq_index] = theta_hhl
    
    # Get real power injection on busses
    P0_bus = B.dot(theta)
    
    # Get real power flow on lines
    P0_line = Bf.dot(theta)
   
    # Get voltage on busses
    n_bus = len(theta)
    V = np.ones(n_bus) * np.exp(1j*theta)
    
    success = 1
    
    return V, success, P0_line, P0_bus, theta  

def lodf_dc(pf0, lodf):
    delta_pf = lodf * pf0
    pf0_aux = np.expand_dims(pf0, 1)
    pfc = pf0_aux + delta_pf
    return pfc, pf0

#%% Functions to export results

def branch_list(branch):
    branches = list(zip(branch['From']+1, branch['To']+1))
    return branches

def get_timestamp_str():
    now = datetime.now()
    timestamp_str = now.strftime("%y%m%d_%H%M%S")
    return timestamp_str

def generate_latex_table_pfc(name, pf0, pfc, branch, export=False):
    # Headers
    branches = list( '{'+ str(branch) + '}' for branch in branch_list(branch) )
    headers = ['line', '{base}'] + branches
    
    # Columns
    num_columns = len(pf0)
    col_definition = 'lS' + 'S' * num_columns
    
    # Create table
    table = r"\begin{table}[H]" + "\n"
    table += r"\centering" + "\n"
    table += r"\caption{pfc-change-this-caption}" + "\n"
    table += r"\label{tab:pfc-change-this-label}" + "\n"
    
    # Create tabular
    table += r"\begin{tabular}{@{}" + col_definition + r"@{}}" + "\n" + r"\toprule" + "\n"
    
    # Add headers
    table += " & ".join(headers) + r" \\" + "\n" + r"\midrule" + "\n"    
    # Add rows
    decimal_precision = 5
    for i, row in enumerate(pfc):
        branch_str = [ branches[i] ]
        row_str = list( str(np.round(val, decimal_precision)) for val in row)
        pf0_str = [ str(np.round(pf0[i], decimal_precision)) ]
        row = branch_str + pf0_str + row_str
        
        table += " & ".join(row) + r" \\" + "\n"
    
    # Close table
    table += r"\bottomrule" + "\n" + r"\end{tabular}" + r"\end{table}"
    # print(table)
    
    # Save table into a .tex file
    if export:
        timestamp = get_timestamp_str()
        with open('postpro/tables/tab_pfc_' + name + '_' + timestamp + '.tex', 'w') as f:
            f.write(table)
        f.close()
    
    return table
    
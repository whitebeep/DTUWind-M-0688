#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon May 29 18:47:20 2023

@author: eduardantoligil
"""

#%% Imports
# import numpy as np

import power.power_systems as ps
import power.power_flow as pf
import postpro.postprocessing as pp

from LODF.aux_func import prepare_lodf_matrix, prepare_pfc_matrix
from LODF.aux_func import save_quantum_circuit, save_quasiprob, save_statevector_as_latex_table
from LODF.aux_func import postpro_statevector, postpro_counts, pfc_from_quantum_circuit

from LODF.circ_func import vector_to_quantum_circuit, matrix_to_quantum_circuit
from LODF.circ_func import compose_quantum_circuit, simulate_quantum_circuit

import os
os.chdir("/Users/eduardantoligil/Documents/Estudis universitaris/DTU Technical University of Denmark/DTU Courses/Master Thesis/code/DTUWind-M-0688/N_1/")

#%% Controls

# General
switch_prints = False
switch_export = True

# Main
switch_dc = True
prints = False
togate = True

#%% Module import status
if switch_prints:
    ps.import_status()
    pf.import_status()
    pp.import_status()
    # hhl_dc.import_status()
    # ibm_cred.import_status()

#%% power_systems.py

# Select Network Data
# network = ps.IEEE_5bus_mod
network = ps.slides_3bus

# Initialize Network
ps.init_network_data() # ps.bus, ps.branch, ps.gen

if switch_prints:
    print('\nInitialize Network:')
    print('\tBus(es): ', ps.bus)
    print('\tBranch(s): ', ps.branch)
    print('\tGenerator(s): ', ps.gen)

# Load Network
ps.load_network_data(network)

if switch_prints:
    print('\nLoaded Network:')
    print('\tBus(es): ', ps.bus)
    print('\tBranch(s): ', ps.branch)
    print('\tGenerator(s): ', ps.gen)

# Load PowerFlow
# ps.load_pf_data('AC_NR')
# ps.load_pf_data('AC_FDLF')
ps.load_pf_data('DC')

# Load PTDF and LODF matrices
ps.load_contingency_data()

#%% power_flow.py
""" DC Power Flow: Classical Solution """

if switch_dc:
    # Run DC Power Flow
    pf.PowerFlowDC(ps.bus, ps.Bdc_bus, ps.Bdc_line, ps.S_bus, ps.pvpq_index)
    # Available variables: pf.V_dc, pf.success_dc, pf.P0_line_dc, pf.P0_bus_dc, pf.theta_dc
    if switch_prints:
        print('\nDC-PF variables: pf.V_dc, pf.success_dc, pf.P0_line_dc, pf.P0_bus_dc, pf.theta_dc')
        print('\nDC Power Flow Results:','\n V_dc: ', pf.V_dc, '\n success_dc: ', pf.success_dc, '\n Pline0_dc: ', pf.P0_line_dc, '\n Pbus0_dc: ', pf.P0_bus_dc, '\n Theta_dc: ', pf.theta_dc)
    
#%% post_processing.py
""" DC Power Flow: Classical Post-processing"""

if switch_dc:
    ##### From Classical DC Power Flow
    # Calculate post-contingency values
    PFC_dc, PF0_dc = pp.lodf_dc(pf.P0_line_dc, ps.lodf)
    # Generate .tex table
    table_pfc_dc = pp.generate_latex_table_pfc('dc', PF0_dc, PFC_dc, ps.branch, switch_export)

#%%  Input
print('\n=== Input ===')

Pinj_noref  = ps.bus['Pld'][ps.pvpq_index]

# Quantum Circuit for input
qc_input, vec_norm = vector_to_quantum_circuit(Pinj_noref, name='input', prints=prints)

# Scaling
scaling_input = vec_norm

# Simulate QC
sv_input, ct_input = simulate_quantum_circuit(qc_input)

# Save QC
save_quantum_circuit(qc_input, togate=togate, export=switch_export)
qc_input.draw('mpl')

#%% PTDF.py
print('\n=== PTDF ===')

# Quantum Circuit for PTDF
qc_mat_ptdf, mat_norm_ptdf = matrix_to_quantum_circuit(ps.ptdf_noref, name='PTDF\nno_ref', togate=togate, prints=prints)

# Compose Quantum Circuit
qc_ptdf = compose_quantum_circuit(qc_input, qc_mat_ptdf)
qc_ptdf.name = 'PTDF'

# Scaling
scaling_ptdf = vec_norm * mat_norm_ptdf

# Simulate QC
sv_ptdf, ct_ptdf = simulate_quantum_circuit(qc_ptdf)

# Extract results from QC
sv_ptdf_result = postpro_statevector(sv_ptdf, scaling=scaling_ptdf, name='PTDF: statevector scaled')
ct_ptdf_result= postpro_counts(ct_ptdf, scaling=scaling_ptdf)
# print('\nsv_ptdf_result {} = \n'.format(sv_ptdf_result.shape), sv_ptdf_result)
# print('\nct_ptdf_result {} = \n'.format(ct_ptdf_result.shape), ct_ptdf_result)

# Draw Quantum Circuit
save_statevector_as_latex_table(sv_ptdf, qc_ptdf, scaling=scaling_ptdf, export=switch_export)
save_quasiprob(qc_ptdf, ct_ptdf, export=switch_export)
save_quantum_circuit(qc_ptdf, togate=togate, export=switch_export)
qc_ptdf.draw('mpl')

#%% LODF.py - LODF matrix
print('\n=== LODF ===')

lodf_matrix = ps.lodf

# Initialize qubits with statevector from PTDF
sv_init = sv_ptdf.data
mat_lodf = prepare_lodf_matrix(lodf_matrix)

# Quantum circuit for input
qc_input_lodf, vec_norm_lodf = vector_to_quantum_circuit(sv_init, name='input-lodf', prints=prints)

# Quantum Circuit for LODF
qc_mat_lodf, mat_norm_lodf = matrix_to_quantum_circuit(mat_lodf, name='LODF\nmatrix', togate=togate, prints=prints) #qc_mat_lodf.draw('mpl')

# Compose circuit
qc_lodf_mat = compose_quantum_circuit(qc_input_lodf, qc_mat_lodf)
qc_lodf_mat.name = 'LODF'

# Scaling
scaling_lodf_mat = scaling_ptdf * mat_norm_lodf

# Simulate QC
sv_lodf_mat, ct_lodf_mat = simulate_quantum_circuit(qc_lodf_mat)

# Extract results from QC
sv_lodf_result_mat = postpro_statevector(sv_lodf_mat, scaling=scaling_lodf_mat, name='LODF: statevector scaled')
ct_lodf_result_mat = postpro_counts(ct_lodf_mat, scaling=scaling_lodf_mat)
# print('\nsv_lodf_result_mat {} = \n'.format(sv_lodf_result_mat.shape), sv_lodf_result_mat)
# print('\nct_lodf_result_mat {} = \n'.format(ct_lodf_result_mat.shape), ct_lodf_result_mat)

# Draw Quantum Circuit
save_statevector_as_latex_table(sv_lodf_mat, qc_lodf_mat, scaling=scaling_lodf_mat, export=switch_export)
save_quasiprob(qc_lodf_mat, ct_lodf_mat, export=switch_export)
save_quantum_circuit(qc_lodf_mat, togate=togate, export=switch_export)
qc_lodf_mat.draw('mpl')

#%% LODF.py - PFC matrix
print('\n=== PFC ===')

# Initialize qubits with statevector from PTDF
sv_init = sv_lodf_mat.data
mat_pfc = prepare_pfc_matrix(lodf_matrix)

# Quantum circuit for input
qc_input_pfc, vec_norm_pfc = vector_to_quantum_circuit(sv_init, name='input-pfc', prints=prints)

# Quantum Circuit for PFC
qc_mat_pfc, mat_norm_pfc = matrix_to_quantum_circuit(mat_pfc, name='PFC\nmatrix', togate=togate, prints=prints) #qc_mat_pfc.draw('mpl')

# Compose circuit
qc_pfc_mat = compose_quantum_circuit(qc_input_pfc, qc_mat_pfc)
qc_pfc_mat.name = 'PFC'

# Scaling
scaling_pfc_mat = scaling_ptdf * mat_norm_lodf * mat_norm_pfc

# Simulate QC
sv_pfc_mat, ct_pfc_mat = simulate_quantum_circuit(qc_pfc_mat)

# Extract results from QC
sv_pfc_result_mat = postpro_statevector(sv_pfc_mat, scaling=scaling_pfc_mat, name='PFC: statevector scaled')
ct_pfc_result_mat = postpro_counts(ct_pfc_mat, scaling=scaling_pfc_mat)
# print('\nsv_pfc_result_mat {} = \n'.format(sv_pfc_result_mat.shape), sv_pfc_result_mat)
# print('\nct_pfc_result_mat {} = \n'.format(ct_pfc_result_mat.shape), ct_pfc_result_mat)

# Postprocess results from QC
qc_pfc, qc_pf0 = pfc_from_quantum_circuit(sv_pfc_result_mat, lodf_matrix)

# Draw Quantum Circuit
save_statevector_as_latex_table(sv_pfc_mat, qc_pfc_mat, scaling=scaling_pfc_mat, export=switch_export)
save_quasiprob(qc_pfc_mat, ct_pfc_mat, export=switch_export)
save_quantum_circuit(qc_pfc_mat, togate=togate, export=switch_export)
qc_pfc_mat.draw('mpl')

#%% LODF.py : N-1 Contingency Assessment
print('\n=== N-1 ===')

# Inputs
Pinj_noref  = ps.bus['Pld'][ps.pvpq_index]
ptdf_noref = ps.ptdf_noref
lodf_matrix = ps.lodf
mat_lodf = prepare_lodf_matrix(lodf_matrix)
mat_pfc = prepare_pfc_matrix(lodf_matrix)

# Quantum circuit for input
qc_input, vec_norm = vector_to_quantum_circuit(Pinj_noref, name='input', prints=prints)

# Quantum circuit for PTDF
qc_mat_ptdf, mat_norm_ptdf = matrix_to_quantum_circuit(ptdf_noref, name='PTDF\nno_ref', togate=togate, prints=prints)

# Quantum Circuit for LODF
qc_mat_lodf, mat_norm_lodf = matrix_to_quantum_circuit(mat_lodf, name='LODF\nmatrix', togate=togate, prints=prints) #qc_mat_lodf.draw('mpl')

# Quantum Circuit for PFC
qc_mat_pfc, mat_norm_pfc = matrix_to_quantum_circuit(mat_pfc, name='PFC\nmatrix', togate=togate, prints=prints) #qc_mat_pfc.draw('mpl')

# Compose circuit
qc_N1 = compose_quantum_circuit(qc_input, qc_mat_ptdf, qc_mat_lodf, qc_mat_pfc)
qc_N1.name = 'N1'

# Scaling
scaling_N1 = vec_norm * mat_norm_ptdf * mat_norm_lodf * mat_norm_pfc

# Simulate QC
sv_N1, ct_N1 = simulate_quantum_circuit(qc_N1)

# Extract results from QC
sv_N1_result = postpro_statevector(sv_N1, scaling=scaling_N1, name='PTDF_LODFmat_PFC Statevector scaled')
ct_N1_result = postpro_counts(ct_N1, scaling=scaling_N1)

print('\nsv_N1_result {} = \n'.format(sv_N1_result.shape), sv_N1_result)
print('\nct_N1_result {} = \n'.format(ct_N1_result.shape), ct_N1_result)

# Postprocess results from QC
qc_pfc, qc_pf0 = pfc_from_quantum_circuit(sv_N1_result, lodf_matrix)

# Draw Quantum Circuit
save_statevector_as_latex_table(sv_N1, qc_N1, scaling=scaling_N1 ,export=switch_export)
save_quasiprob(qc_N1, ct_N1, export=switch_export)
save_quantum_circuit(qc_N1, togate=togate, export=switch_export)
qc_N1.draw('mpl')

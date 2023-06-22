#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat May 6 19:04:32 2023

@author: eduardantoligil
"""

#%% Control
# basis_gates=['id', 'rz', 'sx', 'x', 'cx']
prints = False
export = True
togate = True

#%% Imports 

import numpy as np
from aux_func import prepare_lodf_matrix, prepare_pfc_matrix, chop
from aux_func import save_quantum_circuit, save_quasiprob, save_statevector_as_latex_table
from aux_func import postpro_statevector, postpro_counts, pfc_from_quantum_circuit, pfc_from_quantum_circuit_old
from circ_func import vector_to_quantum_circuit, matrix_to_quantum_circuit
from circ_func import compose_quantum_circuit, simulate_quantum_circuit

import os
os.chdir("/Users/eduardantoligil/Documents/Estudis universitaris/DTU Technical University of Denmark/DTU Courses/Master Thesis/code/DTUWind-M-0688/N_1/")

#%% Test PTDF circuit
print('\n=== Test PTDF ===')

# Prepare inputs
Pinj_noref = np.array([1. , -1.6])
ptdf_noref = np.array([[-0.81818182, -0.45454545],
                        [-0.18181818, -0.54545455],
                        [ 0.18181818, -0.45454545]])
P0_line = ptdf_noref.dot(Pinj_noref)
print('\nptdf_noref {} = \n'.format(ptdf_noref.shape), ptdf_noref)
print('\nPinj_noref {} = \n'.format(Pinj_noref.shape), Pinj_noref)
print('\nP0_line {} = \n'.format(P0_line.shape), P0_line)
# It does not make any difference if the circuit is run with the following matrix:
# ptdf_noref = np.array([[-0.81818182, -0.45454545,  0.        ,  0.        ],
#                         [-0.18181818, -0.54545455,  0.        ,  0.        ],
#                         [ 0.18181818, -0.45454545,  0.        ,  0.        ],
#                         [ 0.        ,  0.        ,  0.        ,  0.        ]])

# Quantum Circuit for input
qc_input, vec_norm = vector_to_quantum_circuit(Pinj_noref, name='input', prints=prints)

# Quantum Circuit for PTDF
qc_mat_ptdf, mat_norm_ptdf = matrix_to_quantum_circuit(ptdf_noref, name='PTDF\nno_ref', togate=togate, prints=prints)

# Compose Quantum Circuit
qc_ptdf = compose_quantum_circuit(qc_input, qc_mat_ptdf)
qc_ptdf.name = 'PTDF'

# Scaling
scaling_ptdf = vec_norm * mat_norm_ptdf

# Simulate QC
sv_ptdf, ct_ptdf = simulate_quantum_circuit(qc_ptdf)

# Extract results from QC
sv_ptdf_result = postpro_statevector(sv_ptdf, scaling=scaling_ptdf, name='PTDF Statevector scaled')
ct_ptdf_result= postpro_counts(ct_ptdf, scaling=scaling_ptdf)

print('\nsv_ptdf_result {} = \n'.format(sv_ptdf_result.shape), sv_ptdf_result)
print('\nct_ptdf_result {} = \n'.format(ct_ptdf_result.shape), ct_ptdf_result)

# Draw Quantum Circuit
save_statevector_as_latex_table(sv_ptdf, qc_ptdf, scaling=scaling_ptdf, export=export)
save_quasiprob(qc_ptdf, ct_ptdf, export=export)
save_quantum_circuit(qc_ptdf, togate=togate, export=export)
qc_ptdf.draw('mpl')

#%% Test LODF row circuit
print('\n=== Test LODF(row) ===')

# Prepare inputs
sv_init = sv_ptdf.data
row_lodf = prepare_lodf_matrix(np.array([[-1.,  1., -1.]]))
print('\nsv_init {} = \n'.format(sv_init.shape), chop(sv_init))
print('\nrow_lodf {} = \n'.format(row_lodf.shape), row_lodf)

# Quantum Circuit for input
qc_input_lodf, vec_norm_lodf = vector_to_quantum_circuit(sv_init, name='input-lodf', prints=prints)

# Quantum Circuit for LODF row
qc_row_lodf, row_norm_lodf = matrix_to_quantum_circuit(row_lodf, name='LODF\nrow', togate=togate, prints=prints)

# Compose Quantum Circuit
qc_lodf_row = compose_quantum_circuit(qc_input_lodf, qc_row_lodf)
qc_lodf_row.name = 'LODFrow'

# Scaling
scaling_lodf_row = scaling_ptdf * row_norm_lodf

# Simulate QC
sv_lodf_row, ct_lodf_row = simulate_quantum_circuit(qc_lodf_row)

# Extract results from QC
sv_lodf_result_row = postpro_statevector(sv_lodf_row, scaling=scaling_lodf_row, name='LODFrow Statevector scaled')
ct_lodf_result_row = postpro_counts(ct_lodf_row, scaling=scaling_lodf_row)

print('\nsv_lodf_result_row {} = \n'.format(sv_lodf_result_row.shape), sv_lodf_result_row)
print('\nct_lodf_result_row {} = \n'.format(ct_lodf_result_row.shape), ct_lodf_result_row)

# Draw Quantum Circuit
save_statevector_as_latex_table(sv_lodf_row, qc_lodf_row, scaling=scaling_lodf_row, export=export)
save_quasiprob(qc_lodf_row, ct_lodf_row, export=export)
save_quantum_circuit(qc_lodf_row, togate=togate, export=export)
qc_lodf_row.draw('mpl')

#%% End-to-End test PTDF+LODF(row)
print('\n=== End-to-End: PTDF+LODF(row) ===')

# Prepare inputs
Pinj_noref = np.array([1. , -1.6])
ptdf_noref = np.array([[-0.81818182, -0.45454545],
                       [-0.18181818, -0.54545455],
                       [ 0.18181818, -0.45454545]])
row_lodf = prepare_lodf_matrix(np.array([[-1.,  1., -1.]]))
print('\nptdf_noref {} = \n'.format(ptdf_noref.shape), ptdf_noref)
print('\nPinj_noref {} = \n'.format(Pinj_noref.shape), Pinj_noref)
print('\nrow_lodf {} = \n'.format(row_lodf.shape), row_lodf)

# Quantum circuit for input
qc_input, vec_norm = vector_to_quantum_circuit(Pinj_noref, name='input', prints=prints)

# Quantum circuit for PTDF
qc_mat_ptdf, mat_norm_ptdf = matrix_to_quantum_circuit(ptdf_noref, name='PTDF\nno_ref', togate=togate, prints=prints)

# Quantum circuit for LODF row
qc_row_lodf, row_norm_lodf = matrix_to_quantum_circuit(row_lodf, name='LODF\nrow', togate=togate, prints=prints)

# Compose Quantum Circuit
qc_main = compose_quantum_circuit(qc_input, qc_mat_ptdf, qc_row_lodf)
qc_main.name = 'PTDF_LODFrow'

# Scaling
scaling = vec_norm * mat_norm_ptdf * row_norm_lodf

# Simulate QC
statevector, counts = simulate_quantum_circuit(qc_main)

# Extract results from QC
sv_result = postpro_statevector(statevector, scaling=scaling, name='PTDF_LODFrow State vector scaled')
ct_result= postpro_counts(counts, scaling=scaling)

print('\nsv_result {} = \n'.format(sv_result.shape), sv_result)
print('\ncounts_result {} = \n'.format(ct_result.shape), ct_result)

# Draw Quantum Circuit
save_statevector_as_latex_table(statevector, qc_main, scaling=scaling, export=export)
save_quasiprob(qc_main, counts, export=export)
save_quantum_circuit(qc_main, togate=togate, export=export)
qc_main.draw('mpl')

#%% Test LODF matrix circuit
print('\n=== Test LODF matrix ===')

lodf_matrix = np.array([[-1.,  1., -1.],
                        [ 1., -1.,  1.],
                        [-1.,  1., -1.]])

# Initialize qubits with statevector from PTDF
sv_init = sv_ptdf.data
mat_lodf = prepare_lodf_matrix(lodf_matrix)
print('\nsv_init {} = \n'.format(sv_init.shape), chop(sv_init))
print('\nmat_lodf {} = \n'.format(mat_lodf.shape), mat_lodf)

# Quantum circuit for input
qc_input_lodf, vec_norm_lodf = vector_to_quantum_circuit(sv_init, name='input-lodf', prints=prints)

# Quantum Circuit for LODF matrix
qc_mat_lodf, mat_norm_lodf = matrix_to_quantum_circuit(mat_lodf, name='LODF\nmatrix', togate=togate, prints=prints) #qc_mat_lodf.draw('mpl')

# Compose circuit
qc_lodf_mat = compose_quantum_circuit(qc_input_lodf, qc_mat_lodf)
qc_lodf_mat.name = 'LODFmat'

# Scaling
scaling_lodf_mat = scaling_ptdf * mat_norm_lodf

# Simulate QC
sv_lodf_mat, ct_lodf_mat = simulate_quantum_circuit(qc_lodf_mat)

# Extract results from QC
sv_lodf_result_mat = postpro_statevector(sv_lodf_mat, scaling=scaling_lodf_mat, name='LODFmat Statevector scaled')
ct_lodf_result_mat = postpro_counts(ct_lodf_mat, scaling=scaling_lodf_mat)

print('\nsv_lodf_result_mat {} = \n'.format(sv_lodf_result_mat.shape), sv_lodf_result_mat)
print('\nct_lodf_result_mat {} = \n'.format(ct_lodf_result_mat.shape), ct_lodf_result_mat)

# Postprocess results from QC
qc_pfc, qc_pf0 = pfc_from_quantum_circuit_old(sv_lodf_result_mat, lodf_matrix)

# Draw Quantum Circuit
save_statevector_as_latex_table(sv_lodf_mat, qc_lodf_mat, scaling=scaling_lodf_mat, export=export)
save_quasiprob(qc_lodf_mat, ct_lodf_mat, export=export)
save_quantum_circuit(qc_lodf_mat, togate=togate, export=export)
qc_lodf_mat.draw('mpl')

#%% Test PFC matrix circuit
print('\n=== Test PFC matrix ===')

lodf_matrix = np.array([[-1.,  1., -1.],
                        [ 1., -1.,  1.],
                        [-1.,  1., -1.]])

# Initialize qubits with statevector from PTDF
sv_init = sv_ptdf.data
mat_lodf = prepare_lodf_matrix(lodf_matrix)
mat_pfc = prepare_pfc_matrix(lodf_matrix)
print('\nsv_init {} = \n'.format(sv_init.shape), chop(sv_init))
print('\nmat_lodf {} = \n'.format(mat_lodf.shape), mat_lodf)

# Quantum circuit for input
qc_input_lodf, vec_norm_lodf = vector_to_quantum_circuit(sv_init, name='input-lodf', prints=prints)

# Quantum Circuit for LODF matrix
qc_mat_lodf, mat_norm_lodf = matrix_to_quantum_circuit(mat_lodf, name='LODF\nmatrix', togate=togate, prints=prints) #qc_mat_lodf.draw('mpl')

# Quantum Circuit for PFC
qc_mat_pfc, mat_norm_pfc = matrix_to_quantum_circuit(mat_pfc, name='PFC\nmatrix', togate=togate, prints=prints) #qc_mat_pfc.draw('mpl')

# Compose circuit
qc_lodf_mat = compose_quantum_circuit(qc_input_lodf, qc_mat_lodf, qc_mat_pfc)
qc_lodf_mat.name = 'LODFmat_PFC'

# Scaling
scaling_lodf_mat = scaling_ptdf * mat_norm_lodf * mat_norm_pfc

# Simulate QC
sv_lodf_mat, ct_lodf_mat = simulate_quantum_circuit(qc_lodf_mat)

# Extract results from QC
sv_lodf_result_mat = postpro_statevector(sv_lodf_mat, scaling=scaling_lodf_mat, name='LODFmat_PFC Statevector scaled')
ct_lodf_result_mat = postpro_counts(ct_lodf_mat, scaling=scaling_lodf_mat)

print('\nsv_lodf_result_mat {} = \n'.format(sv_lodf_result_mat.shape), sv_lodf_result_mat)
print('\nct_lodf_result_mat {} = \n'.format(ct_lodf_result_mat.shape), ct_lodf_result_mat)

# Postprocess results from QC
qc_pfc, qc_pf0 = pfc_from_quantum_circuit(sv_lodf_result_mat, lodf_matrix)

# Draw Quantum Circuit
save_statevector_as_latex_table(sv_lodf_mat, qc_lodf_mat, scaling=scaling_lodf_mat, export=export)
save_quasiprob(qc_lodf_mat, ct_lodf_mat, export=export)
save_quantum_circuit(qc_lodf_mat, togate=togate, export=export)
qc_lodf_mat.draw('mpl')

#%% End-to-End test PTDF+LODF+PFC
print('\n=== End-to-End: PTDF+LODF+PFC ===')

# Prepare inputs
Pinj_noref = np.array([1. , -1.6])
ptdf_noref = np.array([[-0.81818182, -0.45454545],
                       [-0.18181818, -0.54545455],
                       [ 0.18181818, -0.45454545]])
lodf_matrix = np.array([[-1.,  1., -1.],
                        [ 1., -1.,  1.],
                        [-1.,  1., -1.]])
mat_lodf = prepare_lodf_matrix(lodf_matrix)
mat_pfc = prepare_pfc_matrix(lodf_matrix)
print('\nptdf_noref {} = \n'.format(ptdf_noref.shape), ptdf_noref)
print('\nPinj_noref {} = \n'.format(Pinj_noref.shape), Pinj_noref)
print('\nlodf_matrix {} = \n'.format(lodf_matrix.shape), lodf_matrix)
print('\nmat_lodf {} = \n'.format(mat_lodf.shape), mat_lodf)
print('\nmat_pfc {} = \n'.format(mat_pfc.shape), mat_pfc)

# Quantum circuit for input
qc_input, vec_norm = vector_to_quantum_circuit(Pinj_noref, name='input', prints=prints)

# Quantum circuit for PTDF
qc_mat_ptdf, mat_norm_ptdf = matrix_to_quantum_circuit(ptdf_noref, name='PTDF\nno_ref', togate=togate, prints=prints)

# Quantum Circuit for LODF matrix
qc_mat_lodf, mat_norm_lodf = matrix_to_quantum_circuit(mat_lodf, name='LODF\nmatrix', togate=togate, prints=prints) #qc_mat_lodf.draw('mpl')

# Quantum Circuit for PFC
qc_mat_pfc, mat_norm_pfc = matrix_to_quantum_circuit(mat_pfc, name='PFC\nmatrix', togate=togate, prints=prints) #qc_mat_pfc.draw('mpl')

# Compose circuit
qc_N1 = compose_quantum_circuit(qc_input, qc_mat_ptdf, qc_mat_lodf, qc_mat_pfc)
qc_N1.name = 'PTDF_LODFmat_PFC'

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
save_statevector_as_latex_table(sv_N1, qc_N1, scaling=scaling_N1 ,export=export)
save_quasiprob(qc_N1, ct_N1, export=export)
save_quantum_circuit(qc_N1, togate=togate, export=export)
qc_N1.draw('mpl')

#%% Running circuit on a real backend (Outdated)

# # Loading IBMQ Account
# import sys
# sys.path.append('..')
# import hhl.IBM_credentials as ibm_cred
# ibm_cred.load_account()

# # Run on IBMQ backend
# from qiskit import IBMQ

# circuit_real = qc_ptdf.copy()
# circuit_real.measure_all()
# shots = int(2e4)
# switch_retrieve = True
# retrieve_job_id = 'chquib5a2bbvbud368s0' # Dummy job id
# retrieve_backend = 'ibmq_lima'

# if not switch_retrieve:
#     backend_real = ibm_cred.get_real_backend()
#     transpiled_circuit_real = transpile(circuit_real, backend=backend_real, optimization_level=1)
#     job_real = backend_real.run(transpiled_circuit_real, shots=shots)

# else:
#     provider = IBMQ.get_provider(hub='ibm-q')
#     backend_real = provider.get_backend(retrieve_backend)
#     transpiled_circuit_real = transpile(circuit_real, backend=backend_real, optimization_level=1)
#     job_real = backend_real.retrieve_job(retrieve_job_id)
    
# result_real = job_real.result()
# counts_real = result_real.get_counts()
# scaled_counts_real = {state: count/shots for state, count in counts_real.items()}
# P0_line_real = results_from_counts(scaled_counts_real, scaling)

# print('\nP0_line_real = \n', P0_line_real)

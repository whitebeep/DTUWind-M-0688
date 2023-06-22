#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 24 19:04:32 2023

@author: eduardantoligil
"""

#%% Control

# basis_gates=['id', 'rz', 'sx', 'x', 'cx']
prints = False
export = True

#%% Imports

import numpy as np
from aux_func import postpro_statevector, postpro_counts
from aux_func import save_quantum_circuit, save_quasiprob
from circ_func import vector_to_quantum_circuit, matrix_to_quantum_circuit
from circ_func import compose_quantum_circuit, simulate_quantum_circuit

import os
os.chdir("/Users/eduardantoligil/Documents/Estudis universitaris/DTU Technical University of Denmark/DTU Courses/Master Thesis/code/DTUWind-M-0688/N_1/")

#%% Step-by-step: PTDF Quantum Circuit
from aux_func import decompose_to_unitary, compute_control_matrix, pad_to_n_qubits
from qiskit import QuantumRegister, QuantumCircuit

# Prepare inputs
Pinj_noref = np.array([1. , -1.6])
ptdf_noref = np.array([[-0.81818182, -0.45454545],
                        [-0.18181818, -0.54545455],
                        [ 0.18181818, -0.45454545]])

# Quantum Circuit for input
qc_input, vec_norm = vector_to_quantum_circuit(Pinj_noref, name='input', prints=prints)

# Quantum Circuit for PTDF
matrix = ptdf_noref

# Matrix information
mat_norm = np.linalg.norm(matrix)
mat = matrix/mat_norm
mat_shape = mat.shape
n_row = int(np.ceil(np.log2(mat_shape[0])))
n_col = int(np.ceil(np.log2(mat_shape[1])))

## Matrix-to-QuantumCircuit using SVD+LCU
u, s, vh = np.linalg.svd(mat)     # SVD decomposition of mat
s_diag = np.diag(s)
s_diag_aux = np.insert(s_diag, len(s_diag), np.zeros((1,len(s_diag))), axis=0)

# Since s_diag is REAL, extract only UB and VB
UB, VB, UC, VC = decompose_to_unitary(s_diag)

# Pad UB/VB matrices 
UB_pad = pad_to_n_qubits(UB)
VB_pad = pad_to_n_qubits(VB)

# Add control qubit to UB/VB for the LCU
UB_ctrl = compute_control_matrix(UB_pad, 1, ctrl_state=0)
VB_ctrl = compute_control_matrix(VB_pad, 1, ctrl_state=1)

# Combine both matrices into a unitary matrix
UV = UB_ctrl@VB_ctrl

# Pad u/vh matrices
u_pad = pad_to_n_qubits(u)
vh_pad = pad_to_n_qubits(vh)

# Quantum Circuit (QC)
# Initialize Quantum Registers (QR)
qr_aux = qr_aux = QuantumRegister(1, 'aux')
qr_mat = QuantumRegister(max(n_row,n_col), 'mat')

num_op = 6
for op in range(num_op):
    print('\nIteration {}'.format(op))
    print('-------------')
    # Build QC
    qc_mat = QuantumCircuit(qr_mat, qr_aux, name='ptdf_iter_{}'.format(op))
    qc_mat.barrier()
    if op >= 1:
        qc_mat.unitary(vh_pad, qr_mat[:min(n_row, n_col)], label='VH')
    if op >= 2:
        qc_mat.h(qr_aux)
    if op >= 3:
        qc_mat.unitary(UV, qr_aux[:] + qr_mat[:min(n_row, n_col)], label='S_diag')
    if op >= 4:
        qc_mat.h(qr_aux)
    if op >= 5:
        qc_mat.unitary(u_pad, qr_mat[:n_row], label='U')
    
    if prints:
        print('\nmatrix {} = \n'.format(matrix.shape), matrix)
        # print('\nmat_shape = \n', mat_shape)
        # print('\n(n_row, n_col) = \n', (n_row, n_col))
        print('\nmat_norm = \n', np.round(mat_norm, 6))
        print('\nmat {} = \n'.format(mat.shape), mat)
        
        print('\nu {} = \n'.format(u.shape), u)
        print('\ns {} = \n'.format(s.shape), s)
        print('\nvh {} = \n'.format(vh.shape), vh)
        
        print('\ns_diag {} = \n'.format(s_diag.shape), s_diag)
        print('\nUB {} = \n'.format(UB.shape), UB)
        print('\nVB {} = \n'.format(VB.shape), VB)
        
        # Check if UB/UB_pad and VB/VB_pad are equal
        equal_UB = np.allclose(UB, UB_pad)
        equal_VB = np.allclose(VB, VB_pad)
        if not equal_UB:
            print('\nAre UB and UB_pad equal?\n', equal_UB)
            print('\nUB_pad {} = \n'.format(UB_pad.shape), UB_pad)
        if not equal_VB:
            print('\nAre VB and VB_pad equal?\n', equal_VB)
            print('\nVB_pad {} = \n'.format(VB_pad.shape), VB_pad)
            
        print('\nUB_ctrl {} = \n'.format(UB_ctrl.shape), UB_ctrl)
        print('\nVB_ctrl {} = \n'.format(VB_ctrl.shape), VB_ctrl)
        print('\nUV {} = \n'.format(UV.shape), UV)
            
        print('\nu_pad {} = \n'.format(u_pad.shape), u_pad)
        # print('\nUV {} = \n'.format(UV.shape), UV)
        print('\nvh_pad {} = \n'.format(vh_pad.shape), vh_pad)
        
        print('\nqc_mat = \n', qc_mat)
        
    qc_mat_ptdf, mat_norm_ptdf = qc_mat, mat_norm
        
    # Scaling
    scaling_ptdf = vec_norm * mat_norm_ptdf
    
    # Compose Quantum Circuit
    qc_ptdf = compose_quantum_circuit(qc_input, qc_mat_ptdf)
    print('\nqc_ptdf = \n', qc_ptdf)
    
    # Simulate QC
    sv_ptdf, ct_ptdf = simulate_quantum_circuit(qc_ptdf)
    
    # # Extract results from QC
    # sv_ptdf_result = postpro_statevector(sv_ptdf, scaling=scaling_ptdf, name='PTDF Statevector scaled')
    sv_ptdf_result = postpro_statevector(sv_ptdf, scaling=1, name='PTDF Statevector scaled')
    ct_ptdf_result= postpro_counts(ct_ptdf, scaling=scaling_ptdf)
    
    # print('\nsv_ptdf_result = \n', sv_ptdf_result)
    # print('\nct_ptdf_result = \n', ct_ptdf_result)
    
    # # Checks
    # print('\nnp.allclose(u.dot(s_diag_aux.dot(vh))*mat_norm, matrix) = \n', 
    #       np.allclose(u.dot(s_diag_aux.dot(vh))*mat_norm, matrix))
    
    # print('\nnp.allclose(s_diag, .5*(UB+VB)) = \n', 
    #       np.allclose(s_diag, .5*(UB+VB)))
    
    # Draw Quantum Circuit
    save_quantum_circuit(qc_ptdf, togate=False, export=export)
    qc_ptdf.draw('mpl')
    
#%% Test PTDF circuit
print('\nTest PTDF')

# Prepare inputs
Pinj_noref = np.array([1. , -1.6])
ptdf_noref = np.array([[-0.81818182, -0.45454545],
                        [-0.18181818, -0.54545455],
                        [ 0.18181818, -0.45454545]])
P0_line = ptdf_noref.dot(Pinj_noref)
print('\nptdf_noref = \n', ptdf_noref)
print('\nPinj_noref = \n', Pinj_noref)
print('\nP0_line = \n', P0_line)
# It does not make any difference if run the circuit with the following matrix:
# ptdf_noref = np.array([[-0.81818182, -0.45454545,  0.        ,  0.        ],
#                         [-0.18181818, -0.54545455,  0.        ,  0.        ],
#                         [ 0.18181818, -0.45454545,  0.        ,  0.        ],
#                         [ 0.        ,  0.        ,  0.        ,  0.        ]])

# Quantum Circuit for input
qc_input, vec_norm = vector_to_quantum_circuit(Pinj_noref, name='input', prints=prints)

# Quantum Circuit for PTDF
qc_mat_ptdf, mat_norm_ptdf = matrix_to_quantum_circuit(ptdf_noref, name='PTDF', togate=False, prints=prints)

# Compose Quantum Circuit
qc_ptdf = compose_quantum_circuit(qc_input, qc_mat_ptdf)
qc_ptdf.name = 'ptdf'

# Scaling
scaling_ptdf = vec_norm * mat_norm_ptdf

# Simulate QC
sv_ptdf, ct_ptdf = simulate_quantum_circuit(qc_ptdf)

# Extract results from QC
sv_ptdf_result = postpro_statevector(sv_ptdf, scaling=scaling_ptdf, name='PTDF Statevector scaled')
ct_ptdf_result= postpro_counts(ct_ptdf, scaling=scaling_ptdf)

print('\nsv_ptdf_result = \n', sv_ptdf_result)
print('\nct_ptdf_result = \n', ct_ptdf_result)

# Draw Quantum Circuit
save_quasiprob(qc_ptdf, ct_ptdf, export=export)
save_quantum_circuit(qc_ptdf, togate=False, export=export)
qc_ptdf.draw('mpl')    

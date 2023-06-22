#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 3 22:41:11 2023

@author: eduardantoligil
"""

import numpy as np
from aux_func import is_diagonal, decompose_to_unitary, pad_to_n_qubits, compute_control_matrix
from qiskit import QuantumCircuit, QuantumRegister#, ClassicalRegister
from qiskit import Aer, transpile
# from qiskit.visualization import plot_histogram

#%% Encode ptdf matrix into a quantum circuit

def construct_lcu(s_diag, prints=False):
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
    
    if prints:
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
            
    return UV

def vector_to_quantum_circuit(vector, name='input', prints=False):
    # Vector information
    vec_norm = np.linalg.norm(vector)
    vec = vector/vec_norm
    nb_vec = int(np.ceil(np.log2(len(vec))))
    vec_pad = np.pad(vec, (0, int(2**nb_vec)-len(vec)))

    # Quantum Circuit
    # Initialize Quantum Registers (QR)
    qr_vec = QuantumRegister(nb_vec, 'qubit')
    
    # Initialize QC
    qc_vec = QuantumCircuit(qr_vec, name=name) 
    
    # Build QC
    # qc_vec.initialize(vec_pad)
    qc_vec.isometry(vec_pad, qr_vec, None)
    
    if prints:
        print('\nvector {} = \n'.format(vector.shape), vector)
        print('\nvec_norm = \n', vec_norm)
        print('\nvec {} = \n'.format(vec.shape), vec)
        
        print('\nnb_vec = \n', nb_vec)
        #Check if vec/vec_pad are equal
        equal_vec = np.allclose(vec, vec_pad)
        if not equal_vec:
            print('\nvec_pad {} = \n'.format(vec_pad.shape), vec_pad)
            
        print('\nqr_vec = \n', qr_vec)
        print('\nqc_vec = \n', qc_vec)
        
    return qc_vec, vec_norm
    
def arbitrary_matrix_to_quantum_circuit(matrix, name='arb_mat', togate=False, prints=False):
    # Matrix information
    mat_norm = np.linalg.norm(matrix)
    mat = matrix/mat_norm
    mat_shape = mat.shape
    n_row = int(np.ceil(np.log2(mat_shape[0])))
    n_col = int(np.ceil(np.log2(mat_shape[1])))

    ## Matrix-to-QuantumCircuit using SVD+LCU
    u, s, vh = np.linalg.svd(mat)
    s_diag = np.diag(s)
    UV = construct_lcu(s_diag, prints=prints)
    
    # Pad u/vh matrices
    u_pad = pad_to_n_qubits(u)
    vh_pad = pad_to_n_qubits(vh)
    
    # Quantum Circuit (QC)
    # Initialize Quantum Registers (QR)
    qr_aux = qr_aux = QuantumRegister(1, 'aux')
    qr_mat = QuantumRegister(max(n_row,n_col), 'qubit')
    
    # Build QC
    qc_mat = QuantumCircuit(qr_mat, qr_aux, name=name)
    if not togate:
        qc_mat.barrier()
    qc_mat.unitary(vh_pad, qr_mat[:min(n_row, n_col)], label='VH')
    qc_mat.h(qr_aux)
    qc_mat.unitary(UV, qr_aux[:] + qr_mat[:min(n_row, n_col)], label='S_diag')
    qc_mat.h(qr_aux)
    qc_mat.unitary(u_pad, qr_mat[:n_row], label='U')
    
    if togate:
        # Build Gate
        qc_gate = qc_mat.to_gate(label=name)
        # Rebuild QC
        qc_mat = QuantumCircuit(qr_mat, qr_aux, name=name)
        qc_mat.barrier()
        qc_mat.append(qc_gate, qr_mat[:] + qr_aux[:])
    
    if prints:
        print('\nmatrix {} = \n'.format(matrix.shape), matrix)
        # print('\nmat_shape = \n', mat_shape)
        # print('\n(n_row, n_col) = \n', (n_row, n_col))
        print('\nmat_norm = \n', np.round(mat_norm, 6))
        print('\nmat {} = \n'.format(mat.shape), mat)
        
        print('\nu {} = \n'.format(u.shape), u)
        print('\ns {} = \n'.format(s.shape), s)
        print('\nvh {} = \n'.format(vh.shape), vh)
            
        print('\nu_pad {} = \n'.format(u_pad.shape), u_pad)
        # print('\nUV {} = \n'.format(UV.shape), UV)
        print('\nvh_pad {} = \n'.format(vh_pad.shape), vh_pad)
        
        print('\nqc_mat = \n', qc_mat)
    
    return qc_mat, mat_norm

def diagonal_matrix_to_quantum_circuit(matrix, name='diag_mat', togate=False, prints=False):
    # Matrix information
    mat_norm = np.linalg.norm(matrix)
    mat = matrix/mat_norm
    
    ## matrix-to-QuantumCircuit using LCU
    n_qb = int(np.ceil(np.log2(len(matrix)))) + 1 # (+1) from aux_qubit
    mat_pad = pad_to_n_qubits(mat, n=n_qb, pad_val=0)
    UV = construct_lcu(mat_pad, prints=prints) 
    
    # Quantum Circuit (QC)
    # Initialize Quantum Registers (QR)
    qr_aux = QuantumRegister(1, 'aux') 
    qr_mat = QuantumRegister(n_qb, 'qubit')
    
    # Build Gate
    qc_mat = QuantumCircuit(qr_mat, qr_aux, name=name)
    if not togate:
        qc_mat.barrier() 
    qc_mat.h(qr_aux)
    qc_mat.unitary(UV, qr_aux[:] + qr_mat[:], label='LODF')
    qc_mat.h(qr_aux)
    
    if togate:
        # Build Gate
        qc_gate = qc_mat.to_gate(label=name)
        # Rebuild QC
        qc_mat = QuantumCircuit(qr_mat, qr_aux, name=name)
        qc_mat.barrier()
        qc_mat.append(qc_gate, qr_mat[:] + qr_aux[:])
    
    if prints:        
        print('\nmatrix = \n', matrix)
        print('\nmat_norm = \n', np.round(mat_norm, 6))
        print('\nmat = \n', np.round(mat, 6))
        
        print('\nn_qb = \n', n_qb)
        print('\nmat_pad = \n', np.round(mat_pad,6))
        
        print('\nqc_mat = \n', qc_mat)

    return qc_mat, mat_norm

def matrix_to_quantum_circuit(matrix, name='ptdf', togate=False, prints=False):
    if is_diagonal(matrix):
        qc_mat, mat_norm = diagonal_matrix_to_quantum_circuit(matrix, name, togate, prints)
    else:
        qc_mat, mat_norm = arbitrary_matrix_to_quantum_circuit(matrix, name, togate, prints)
    return qc_mat, mat_norm

def compose_quantum_circuit(*qc_list):
    if len(qc_list) == 1:
        return qc_list[0].copy()
    
    qc_first = qc_list[0].copy()
    qc_rest = compose_quantum_circuit(*qc_list[1:])
    
    qc_out = qc_rest.compose(qc_first, front=True)
    return qc_out

def simulate_quantum_circuit(qc, backend=Aer.get_backend('statevector_simulator'), prints=False):
    circuit = qc.copy()
    # Run Quantum Circuit
    transpiled_circuit = transpile(circuit, backend=backend, basis_gates=['id', 'rz', 'sx', 'x', 'cx'])
    job = backend.run(transpiled_circuit)
    result = job.result()
    # Extract Statevector and counts
    sv = result.get_statevector()
    counts = result.get_counts()
    
    if prints:
        # print('\ntranspiled_circuit = \n', transpiled_circuit)
        print('\ntranspiled_circuit depth = \n', transpiled_circuit.depth())
    
    return sv, counts

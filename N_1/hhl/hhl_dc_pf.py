#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 19:12:12 2023

@author: eduardantoligil

Code inspired by Brynjar Sævarsson (@autor: brysa)
"""

#%% Imports

import numpy as np

#%% Import Module

def import_status():
    print('Module imported successfully: hhl/hhl_dc_pf.py')
    return

#%% Functions to solve HHL for DC-PF

# Imports for solve_hhl_sv()
from qiskit import Aer, execute
from qiskit.algorithms.linear_solvers.hhl import HHL

backend_sim = Aer.get_backend("statevector_simulator") # backend_sim = Aer.get_backend('aer_simulator')
hhl = HHL(1e-3) # The default is Statevector simulation , quantum_instance=backend_sim

def matrix_vector(Sbus, pvpq_index, B):
    pvpq = pvpq_index
    Pbus = np.real(Sbus)
    # Initialize A and b from Ax=b
    Pbus_noref = Pbus[pvpq] # b
    B_noref = B[:, pvpq][pvpq,:] # A
    # Assign matrix and vector
    matrix = B_noref
    vector = Pbus_noref
    return matrix, vector

def solve_hhl_sv(matrix, vector):
    bNorm = np.linalg.norm(vector)
    
    num_qubits = int(np.log2(vector.shape[0])) # number of variables (assuming n=2^k)
    
    hhl_circuit = hhl.construct_circuit(matrix, vector, neg_vals=False)

    total_qubits = hhl_circuit.num_qubits
    job = execute(hhl_circuit, backend = backend_sim)
    state = job.result().get_statevector().data
    # print('State: ', state)
    
    amplitudes = np.zeros(len(vector),dtype=complex)
    for i, a in enumerate(state):
    
        # get binary representation
        b = ("{0:0%sb}" % total_qubits).format(i)
        amplitude = a #np.abs(a)#*np.sign(np.real(-a))# ** 2
        
        # extract value of Z and corresponding probability
        i_normal = int(b[-num_qubits:], 2)

        if int(b[0], 2) == 1:

            amplitudes[i_normal] += bNorm*amplitude
            # print('Iteration: ', i)
            # print('b: ', b)
            # print('bNorm: ', bNorm)
            # print(amplitude)
            # print('bNorm*amplitude: ', bNorm*amplitude)
            # print('abs(bNorm*amplitude): ', np.abs(bNorm*amplitude))
            # print('')
        
        # nl_qubits = total_qubits - (num_qubits+1)
        # if int(b[0], 2) == 1 and b[-(num_qubits+nl_qubits):-num_qubits]=='0'*nl_qubits:
        #     amplitudes[i_normal] += bNorm*amplitude     
        
    # print('Result: ', np.abs(amplitudes)*np.sign(np.real(amplitudes)))
    return np.abs(amplitudes)*np.sign(np.real(amplitudes))

#%% HHL Test

# from qiskit.algorithms.linear_solvers import HHL
# from qiskit.quantum_info import Statevector
# # from linear_solvers import HHL, NumPyLinearSolver

# def get_solution_vector(hhl_solution):
#     """Extracts and normalizes simulated state vector
#     from LinearSolverResult."""
#     solution_vector = Statevector(hhl_solution.state).data[16:18].real
#     norm = hhl_solution.euclidean_norm
#     return norm * solution_vector / np.linalg.norm(solution_vector)

# matrix = np.array([ [1, -1/3], [-1/3, 1] ])
# vector = np.array([1, 0])

# print(matrix)
# print(vector)

# hhl_solution = HHL().solve(matrix, vector)
# classical_solution = np.linalg.solve(matrix, vector)

# print('full naive solution vector:', get_solution_vector(hhl_solution))
# print('classical state:', classical_solution)

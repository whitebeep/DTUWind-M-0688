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

# Imports for solve_hhl_tomo()
import time
import hhl.IBM_credentials as ibm
from qiskit.compiler import transpile
from qiskit.tools.monitor import job_monitor
from qiskit.ignis.verification.tomography import state_tomography_circuits, StateTomographyFitter

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
    
def solve_hhl_tomo(matrix, vector, sim=True, error_correction=False):
    
    if error_correction:
        Niter = 10
    else:
        Niter = 1
        
    astar = 0
    tol = 1e-6
    for n in range(Niter):
        bNorm = np.linalg.norm(vector)
    
        hhl_circuit = hhl.construct_circuit(matrix, vector, neg_vals=False)
        t = time.time()
    
        num_qubits = int(np.log2(vector.shape[0])) # number of variables (assuming n=2^k)
        total_qubits = hhl_circuit.num_qubits
        
        qubits = []
        for k in range(num_qubits):
            qubits.append(k)
        qubits.append(total_qubits-1)
        
        qst = state_tomography_circuits(hhl_circuit, qubits)
        Shots = 2**12
        if sim:
            circuitToRun = transpile(qst, basis_gates=['id', 'rz', 'sx', 'x', 'cx'], optimization_level=1)
            job = execute(circuitToRun, backend = Aer.get_backend('qasm_simulator'), shots = Shots)
        else:
            backend_real = ibm.get_real_backend()
            circuitToRun = transpile(qst, basis_gates=['id', 'rz', 'sx', 'x', 'cx'], optimization_level=1)
            job = execute(circuitToRun, backend = backend_real, shots = Shots)
            job_monitor(job)
           
        depth = circuitToRun[0].depth()
        ops = circuitToRun[0].count_ops()
        print('Time taken:', time.time() - t)
        
        tomo = StateTomographyFitter(job.result(), qst)
        
        rho = tomo.fit()
        
        amp = np.zeros(len(vector),dtype=complex)
        signs = np.zeros(len(vector),dtype=complex)
        for i, a in enumerate(rho):
    
            # get binary representation
            b = ("{0:0%sb}" % (num_qubits+1)).format(i)
             
            i_normal = int(b[-num_qubits:], 2)
            if int(b[0], 2) == 1:
    
                amp[i_normal] += bNorm*np.sqrt(a[i])
    
                signs[i_normal] += rho[2**(num_qubits),i]
                
                
        amplitudes = abs(amp)*np.sign(signs.real)*np.sign(vector[0])
        
        # Error correction
        astar = astar + amplitudes
        error = np.linalg.norm(matrix@amplitudes - vector)
        print('Error: ', error)
        if error <= tol: 
            break
        vector = vector-matrix@amplitudes 

    return astar, depth, ops

def solve_random(matrix, vector):
    
    Niter = 1000

    astar = 0
    tol = 1e-6
    for n in range(Niter):
        bNorm = np.linalg.norm(vector)

        amplitudes = np.random.random(len(vector))*vector
        
        # Error correction
        astar += amplitudes
        error = np.linalg.norm(matrix@amplitudes - vector)
        print('\nIteration nº: %d' % n)
        print('vector_0:\t ', vector)
        print('bNorm:\t ', bNorm)
        print('amplitudes:\t ', amplitudes)
        print('astar:\t ', astar)
        print('error:\t ', error)
        
        # print('Error: ', error)
        if error <= tol: 
            print(n)
            break
        vector = vector-matrix@amplitudes 
        print('vector_1:\t ', vector)

    return astar

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

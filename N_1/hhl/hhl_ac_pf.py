#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 19:13:02 2023

@author: eduardantoligil

Code inspired by Brynjar Sævarsson (@autor: brysa)
"""

#%% Imports

import numpy as np
from typing import Optional, Union, List

from qiskit import QuantumCircuit, QuantumRegister, ClassicalRegister
from qiskit import execute
from qiskit import Aer

from qiskit.algorithms.linear_solvers.matrices.numpy_matrix import NumPyMatrix
from qiskit.compiler import transpile

from qiskit.circuit.library import PhaseEstimation
from qiskit.circuit.library.arithmetic.exact_reciprocal import ExactReciprocal

from qiskit.ignis.verification.tomography import state_tomography_circuits, StateTomographyFitter

#%% HHL for AC FDLF functions

# Simplified version of the HHL from qiskit on "qiskit.algorithms.linear_solvers.hhl"
# The circuit is splitted into two functions so only the input part of the HHL circuit needs to be updated

def get_delta(n_l: int, lambda_min: float, lambda_max: float) -> float:
    """Calculates the scaling factor to represent exactly lambda_min on nl binary digits.
    Args:
        n_l: The number of qubits to represent the eigenvalues.
        lambda_min: the smallest eigenvalue.
        lambda_max: the largest eigenvalue.
    Returns:
        The value of the scaling factor.
    """
    formatstr = "#0" + str(n_l + 2) + "b"
    lambda_min_tilde = np.abs(lambda_min * (2**n_l - 1) / lambda_max)
    # floating point precision can cause problems
    if np.abs(lambda_min_tilde - 1) < 1e-7:
        lambda_min_tilde = 1
    binstr = format(int(lambda_min_tilde), formatstr)[2::]
    lamb_min_rep = 0
    for i, char in enumerate(binstr):
        lamb_min_rep += int(char) / (2 ** (i + 1))
    return lamb_min_rep

def construct_circuit(
    matrix: Union[List, np.ndarray, QuantumCircuit],
    vector: Union[List, np.ndarray, QuantumCircuit],
    neg_vals: Optional[bool] = True,
    tomography: Optional[bool] = True,
    ol: Optional[int] = 3,
) -> QuantumCircuit:
    
    nb = int(np.ceil(np.log2(len(vector))))
    vector_circuit = QuantumCircuit(nb,name='input')
    
    if all(vector==0):
        vector+=1e-6
    
    vector_circuit.isometry(vector / np.linalg.norm(vector), list(range(nb)), None)
    
    matrix_circuit = NumPyMatrix(matrix, evolution_time=2*np.pi)
    
    # Set the tolerance for the matrix approximation
    if hasattr(matrix_circuit, "tolerance"):
        matrix_circuit.tolerance = 1e-2 / 6

    # check if the matrix can calculate the condition number and store the upper bound
    if (
        hasattr(matrix_circuit, "condition_bounds")
        and matrix_circuit.condition_bounds() is not None
    ):
        kappa = matrix_circuit.condition_bounds()[1]
    else:
        kappa = 1
    # Update the number of qubits required to represent the eigenvalues
    # The +neg_vals is to register negative eigenvalues because
    # e^{-2 \pi i \lambda} = e^{2 \pi i (1 - \lambda)}
    nl = max(nb + 1, int(np.log2(kappa)) + 1) + neg_vals

    # check if the matrix can calculate bounds for the eigenvalues
    if hasattr(matrix_circuit, "eigs_bounds") and matrix_circuit.eigs_bounds() is not None:
        lambda_min, lambda_max = matrix_circuit.eigs_bounds()
        # Constant so that the minimum eigenvalue is represented exactly, since it contributes
        # the most to the solution of the system. -1 to take into account the sign qubit
        delta = get_delta(nl - neg_vals, lambda_min, lambda_max)
        # Update evolution time
        matrix_circuit.evolution_time = 2 * np.pi * delta / lambda_min / (2**neg_vals)
        # Update the scaling of the solution
        # scaling = lambda_min
    else:
        delta = 1 / (2**nl)
        print("The solution will be calculated up to a scaling factor.")

    reciprocal_circuit = ExactReciprocal(nl, delta, neg_vals=neg_vals)

    # Initialise the quantum registers
    qb = QuantumRegister(nb,'nb')  # right hand side and solution
    ql = QuantumRegister(nl,'nl')  # eigenvalue evaluation qubits

    qf = QuantumRegister(1,'na')  # flag qubits


    qc = QuantumCircuit(qb, ql, qf)
    
    ic = QuantumCircuit(qb, ql, qf)
    # State preparation
    ic.append(vector_circuit, qb[:])
    qc.barrier()
    # QPE
    phase_estimation = PhaseEstimation(nl, matrix_circuit)

    qc.append(phase_estimation, ql[:] + qb[:])

    # Conditioned rotation
    qc.append(reciprocal_circuit, ql[::-1] + [qf[0]])

    # QPE inverse
    qc.append(phase_estimation.inverse(), ql[:] + qb[:])
    
    if tomography:
        qubits = np.arange(nb+1).tolist()
        qubits[-1] = nb+nl # TODO it is probably not necessary to perform tomography on the ancilla qubit but the current code can't run without it
        qc = state_tomography_circuits(qc, qubits)

        qc = transpile(qc,basis_gates=['id', 'rz', 'sx', 'x', 'cx'], optimization_level=ol)
        ic.add_register(qc[0].cregs[0])
    else:
        qc = transpile(qc,basis_gates=['id', 'rz', 'sx', 'x', 'cx'], optimization_level=ol)

    return qc, ic

def update_circuit(
        qc: Union[List, np.ndarray, QuantumCircuit],
        ic: Union[List, np.ndarray, QuantumCircuit],
        vec: Union[List, np.ndarray, QuantumCircuit],
        ol: Optional[int] = 1,
        ) -> QuantumCircuit:
    """Update the input circuit and combine it with the constant part of the HHL circuit.
    Args:
        matrix: The matrix specifying the system, i.e. A in Ax=b.
        vector: The vector specifying the right hand side of the equation in Ax=b.
        neg_vals: States whether the matrix has negative eigenvalues. If False the
        computation becomes cheaper.
    Returns:
        The HHL circuit.
    Raises:
        ValueError: If the input is not in the correct format.
        ValueError: If the type of the input matrix is not supported.
    """
    
    # Remove the old vector_circuit
    if len(ic.data)>0:
        ic.data.pop(0)
    
    # Create a new vector_circuit and add it to the input 
    nb = int(np.log2(len(vec)))
    vector_circuit = QuantumCircuit(nb,name='input')
    vector_circuit.isometry(vec / np.linalg.norm(vec), list(range(nb)), None)
    
    ic.append(vector_circuit,range(nb))

    # Transpile the circuit into basis gates
    # ic = transpile(ic,basis_gates=['id', 'rz', 'sx', 'x', 'cx'], optimization_level=ol)
    
    if type(qc) == list: # including state tomography 
        circuit = []
        for circ in qc:
            # c = ic+circ # XXX this is deprecated.
            ind = [ic.qubits.index(q) for q in circ.qubits]
            c = ic.compose(circ,qubits = [ic.qubits[i] for i in ind])
            
            c.name = circ.name
            circuit.append(c)
    else:    
        circuit = ic.compose(qc)
        cr = ClassicalRegister(nb+1, 'c')
        circuit.add_register(cr)
        circuit.barrier()
        qubits = np.arange(nb+1).tolist()
        qubits[-1] = qc.num_qubits-1
        circuit.measure(qubits,np.arange(nb+1).tolist())
    
    return circuit



#%% Function for solving HHL

def solve_hhl(mat, vec, hhl_cir, ic, backend, tomography=True, neg_vals=False, Shots = 2**12):
    
    bNorm = np.linalg.norm(vec)

    true_solution = np.linalg.solve(mat, vec if neg_vals else -vec) # XXX when tomography is not used we get the sign of the result classically

    num_qubits = int(np.log2(vec.shape[0])) # number of variables (assuming n=2^k)
    
    # if not neg_vals and np.all(np.linalg.eigvals(mat)<0):
    #     vec = -vec

    hhl_circuit = update_circuit(hhl_cir, ic, vec)
    
    n_var = len(vec)
    approx_solution = np.zeros(n_var)
    
    if tomography:
        
        depth = hhl_circuit[0].depth()
        ops = hhl_circuit[0].count_ops()
        
        
        job = execute(hhl_circuit, backend = backend, shots = Shots)
        tomo = StateTomographyFitter(job.result(), hhl_circuit)
        
        rho = tomo.fit()
        
        amp = np.zeros(len(vec),dtype=complex)
        signs = np.zeros(len(vec),dtype=complex)
        for i, a in enumerate(rho):
    
            # get binary representation
            b = ("{0:0%sb}" % (num_qubits+1)).format(i)
             
            i_normal = int(b[-num_qubits:], 2)
            if int(b[0], 2) == 1:
    
                amp[i_normal] += bNorm*np.sqrt(a[i]).real
    
                signs[i_normal] += rho[2**(num_qubits),i]
            
            
        approx_solution = amp.real*np.sign(signs.real)*np.sign(vec[0])
        
        if np.all(np.sign(approx_solution)==np.sign(true_solution)):
            print('True')
        else:
            print('False')
            print('vec: ',vec)
            print('amp: ',approx_solution)
            print('true: ', true_solution)
            
    else:
        
        # circuitToRun = transpile(hhl_circuit,basis_gates=['id', 'rz', 'sx', 'x', 'cx'],optimization_level=3)
        circuitToRun = transpile(hhl_circuit,basis_gates=['id', 'rz', 'sx', 'x', 'cx'])
        
        depth = circuitToRun.depth()
        ops = circuitToRun.count_ops()
        result_indx = np.arange(2 ** (num_qubits), 2 ** (num_qubits+1))

        job = execute(circuitToRun, backend = backend, shots = Shots)
        
        # job_monitor(job)
    
        job_result = job.result()
        M = job_result.get_counts()
    
       
        for i in range(n_var):
            try:
                approx_solution[i] = bNorm*np.sqrt(M[str(np.binary_repr(result_indx[i]))]/Shots)*np.sign(true_solution[i])
            except:
                pass
            finally:
                pass
    
    return approx_solution, job_result, depth, ops


#%% Quantum Power Flow

def HHL_FDLF(Ybus, Bp, Bpp, Sbus, V0, pvpq_index, pq_index, backend=Aer.get_backend('qasm_simulator'), max_iter=100, err_tol=1e-5):
    
    # Iterations dictionary
    dict_keys = ['n', 'V', 'Vm', 'Va', 'normP', 'normQ']
    dict_iter = {key: [] for key in dict_keys}
    # jobs,depth_p, depth_q, times, gates_p, gates_q, Vm_iter, Th_iter, dp_iter, dq_iter
    
    # Flags and iteration counter
    neg_vals = False # Using negative eigenvalues requires more qubits and gates
    tomography = False # Tomography is used to extract the full quantum state but requires more runs. If False then we need to find the sign of the result classically
    Q_flag = False
    success = False # Check if the algorithm is successful
    n = 0
    dict_iter['n'].append(n)

    # Assign initial voltage guess (V0)
    V = V0.copy()
    Vm = abs(V)
    Va = np.angle(V)
    dict_iter['V'].append(V)
    dict_iter['Vm'].append(Vm)
    dict_iter['Va'].append(Vm)
    
    # Active/Reactive power indices
    dP, dQ = pvpq_index, pq_index
    
    # Calculate dpf and dfq power mismatch
    mis = (Sbus - V0 * np.conj(Ybus @ V0))/ Vm
    dfp = np.real(mis[dP])
    dfq = np.imag(mis[dQ])
    
    # Check stop condition
    normP = np.linalg.norm(dfp, np.inf)
    normQ = 0
    if len(pq_index) > 0:
        Q_flag = True
        normQ = np.linalg.norm(dfq, np.inf)
    dict_iter['normP'].append(normP)
    dict_iter['normQ'].append(normQ)
        
    if normP < err_tol and normQ < err_tol:
        success = True
        
    else:
        # Create the initial HHL circuit # Bp == Bpp
        if neg_vals:
            hhl_cir, ic = construct_circuit(Bp, dfp, neg_vals=neg_vals, tomography=tomography)
        else:
            hhl_cir, ic = construct_circuit(-Bp, -dfp, neg_vals=neg_vals, tomography=tomography)
    
    while not success and n < max_iter:
        # Iteration start
        n += 1
        dict_iter['n'].append(n)
        
        ## Active Power (Bp)
        # Assign b vector
        vec = dfp if neg_vals else -dfp
        
        # Compute delta voltage angles (dVa)
        approx_dVa, job_result, depth, ops = solve_hhl(Bp, vec, hhl_cir, ic, backend, tomography=tomography, neg_vals=neg_vals)
        dVa = -approx_dVa   
        # classic_dVa = np.linalg.solve(Bp, dfp)
        
        # Update voltage angles (Va)
        Va[dP]  += dVa
        dict_iter['Vm'].append(['-']*len(Va))
        dict_iter['Va'].append(Va)

        
        # Update voltage vector
        V = Vm * np.exp(1j * Va)
        dict_iter['V'].append(V)
    
        # Calculate dfp and dfq mismatch
        mis = (Sbus - V * np.conj(Ybus @ V)) / Vm
        dfp = np.real(mis[dP])
        dfq = np.imag(mis[dQ])
        
        # Check stop condition    
        normP = np.linalg.norm(dfp, np.inf)
        if Q_flag:
            normQ = np.linalg.norm(dfq, np.inf)
        dict_iter['normP'].append(normP)
        dict_iter['normQ'].append(normQ)
            
        if normP < err_tol and normQ < err_tol:
            success = True
            break
        
        ## Reactive Power (Bpp)
        if Q_flag:
            # Same iteration
            dict_iter['n'].append(n)
            
            # Assign b vector
            vec = dfq if neg_vals else -dfq
            
            # Compute delta voltage magnitudes (dVm)
            approx_dVm, job_result, depth, ops = solve_hhl(Bpp, vec, hhl_cir, ic, backend, tomography=tomography, neg_vals=neg_vals)
            dVm = -approx_dVm
            # classic_dVm = np.linalg.solve(Bpp, dfq)
            
            # Update voltage magnitudes (Vm)
            Vm[dQ] += dVm
            dict_iter['Vm'].append(Vm)
            dict_iter['Va'].append(['-']*len(Vm))
            
            # Update voltage vector
            V = Vm * np.exp(1j * Va)
        
            # Compute dfp and dfq mismatch
            mis = (Sbus - V * np.conj(Ybus @ V)) / Vm
            dfp = np.real(mis[dP])
            dfq = np.imag(mis[dQ])
        
            # Check stop condition
            normP = np.linalg.norm(dfp, np.inf)  
            normQ = np.linalg.norm(dfq, np.inf)
            dict_iter['normP'].append(normP)
            dict_iter['normQ'].append(normQ)
                
            if normP < err_tol and normQ < err_tol:
                success = True
                break
        else:
            pass

    return V, success, n, dict_iter


#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed May 3 22:42:31 2023

@author: eduardantoligil
"""

import numpy as np
import scipy as sp
from qiskit.visualization import plot_histogram

#%% aux_func

def pad_to_square_matrix(array, nqb=None):
    if nqb is None:
        nqb = int(np.ceil(np.log2(max(array.shape))))
    ndim = 2**nqb
    nr, nc = array.shape
    pad = ((0,ndim-nr),(0,ndim-nc))
    array_pad = np.pad(array, pad)
    return array_pad

def prepare_lodf_matrix(array):
    # Create a diagonal matrix for each lodf row
    # diag_matrices = [np.diag(row) for row in array]  
    diag_matrices = [np.eye(len(array[0]))]+[np.diag(row) for row in array]
    # Stack diag matrices verticalls
    mat = np.vstack(diag_matrices)
    mat = pad_to_square_matrix(mat)
    if len(array) != 1:
        # Pad to a square matrix
        mat = pad_to_square_matrix(mat)
    return mat

def prepare_pfc_matrix(array):
    n_p0 = array.shape[0]
    n_cas = array.size
    n_tot = n_p0 + n_cas
    nqb = int(np.ceil(np.log2(n_tot))) + 1 # (+1) from aux_qubit
    mat = np.eye(n_tot)
    for i,n in enumerate(np.arange(n_p0, n_tot, n_p0)):
        mat[n:n+n_p0, i] = 1
    mat = pad_to_square_matrix(mat, nqb=nqb)
    return mat

def postpro_statevector(sv, scaling=1, name='State vector', decimals=6, prints=True):
    array = scaling * chop(sv.data)
    n_qubits = int(np.ceil(np.log2(array.size)))
    if prints:
        print('\n', '{}:'.format(name))
        for i, a in enumerate(array):
            str_bin = ("{:0%sb}" % n_qubits).format(i)
            print('\t{}: {}'.format(str_bin, round(a, decimals)))      
    res = array[:2**(n_qubits-1)]
    return res.real

def postpro_counts(counts, scaling=1):
    bit_str = list(counts.keys())
    res_qubits = len(bit_str[0])-1
    res = np.zeros(2**res_qubits)
    
    for i, a in counts.items():
        if i[0]=='0':
            res_idx = int(i[1:], 2)
            res[res_idx] += np.sqrt(a)
    
    res *= scaling
    return res

def pfc_from_quantum_circuit(qc_res, lodf_mat, prints=True): # Currently implemented for qc_res==statevector
    # Extract lodf size and shape
    lodf_size = lodf_mat.size
    lodf_shape = lodf_mat.shape

    # Extract PF0 and PFC from qc_res
    qc_PF0 = qc_res[:lodf_shape[0]]
    qc_PFC_row = qc_res[lodf_shape[0]:lodf_shape[0]+lodf_size]
    qc_PFC = qc_PFC_row.reshape(lodf_shape)
    
    if prints:
        print('\nqc_PF0 = \n', qc_PF0)
        print('\nqc_PFC = \n', qc_PFC)
    return qc_PFC, qc_PF0

def pfc_from_quantum_circuit_old(qc_res, lodf_mat, prints=True): # Currently implemented for qc_res==statevector
    # Extract lodf size and shape
    lodf_size = lodf_mat.size
    lodf_shape = lodf_mat.shape

    # Extract PF0 and deltaPF from qc_res
    qc_PF0 = qc_res[:lodf_shape[0]]
    qc_deltaPF_row = qc_res[lodf_shape[0]:lodf_shape[0]+lodf_size]
    
    # Prepare PF0 and deltaPF to obtain PFC
    qc_PF0_aux = np.expand_dims(qc_PF0, 1)
    qc_deltaPF = qc_deltaPF_row.reshape(lodf_shape)
    
    # Calculate PFC
    qc_PFC = chop(qc_PF0_aux + qc_deltaPF)
    
    if prints:
        print('\nqc_PF0 = \n', qc_PF0)
        print('\nqc_deltaPF = \n', qc_deltaPF)
        print('\nqc_PFC = \n', qc_PFC)
    return qc_PFC, qc_PF0

def save_quantum_circuit(qc, output='mpl', togate=False, export=False):
    if togate:
        # filename ='draw/togate_'+qc.name+'.png'
        filename ='LODF/draw/togate_'+qc.name+'.png'
    else:
        # filename = 'draw/'+qc.name+'.png'
        filename ='LODF/draw/'+qc.name+'.png'
    options = {'output': output,
               'filename': filename,
               'plot_barriers': True,
               'initial_state': False
               }
    if export:
        qc.draw(**options)
        print('\nSaved Quantum Circuit: {}'.format(qc.name))
    return

def save_quasiprob(qc, qc_count, export=False):
    options = {'data': qc_count,
               'color': '#1F3DFF',
               # 'filename': 'hist/quasiprob_'+qc.name+'.png'
               'filename': 'LODF/hist/quasiprob_'+qc.name+'.png'
               }
    if export:
        plot_histogram(**options)
        print('\nSaved Quasi-probabilities: {}'.format(qc.name))
    return

def save_statevector_as_latex_table(statevector, qc, scaling=1, export=False):
    n_qb = statevector.num_qubits
    n_val = len(statevector.data)
    n_col = 2  # Number of columns per row (change to 4)
    n_row = (n_val + n_col - 1) // n_col  # Calculate number of rows
    headers = n_col * ['Bitstring & {Amplitude}']
    caption = 'Statevector Table'
    label = 'sv_for_' + qc.name

    # Create table
    table = r"\begin{table}[H]" + "\n"
    table += r"\centering" + "\n"
    table += r"\caption{" + caption + "}" + "\n"
    table += r"\label{tab:" + label + "}" + "\n"

    # Create tabular
    table += r"\begin{tabular}{@{}" + "rl" * n_col + r"@{}}" + "\n" + r"\toprule" + "\n"

    # Add headers
    table += " & ".join(headers) + r" \\" + "\n" + r"\midrule" + "\n"

    for i in range(n_row):
        s_idx = i * n_col
        e_idx = min((i + 1) * n_col, n_val)
        if scaling != 1:
            row_val = scaling * chop(statevector.data[s_idx:e_idx])
        else: 
            row_val = chop(statevector.data[s_idx:e_idx])
        row_bit = [format(j, '0{}b'.format(n_qb)) for j in range(s_idx, e_idx)]
        row_val_rounded = [0 if (v.real==0 and v.imag==0) else np.round(v, decimals=4) for v in row_val]

        table += " & ".join(["{} & {}".format(b, '{-}' if v==0 else str(v)) for b, v in zip(row_bit, row_val_rounded)])
        table += r" \\" + "\n"

    # Close table
    table += r"\bottomrule" + "\n" + r"\end{tabular}" + r"\end{table}"

    if export:
        # filename = 'sv/sv_' + qc.name + '.tex'
        filename = 'LODF/sv/sv_' + qc.name + '.tex'
        with open(filename, 'w') as f:
            f.write(table)
        f.close()
        print('\nSaved Statevector: {}'.format(qc.name))
        
    return table


def chop(array, epsilon=1e-12):
    """
    Function copied from qiskit docs 0.19.6 (currently deprecated)
    Truncate small values of a complex array.

    Args:
        array (array_like): array to truncate small values.
        epsilon (float): threshold.

    Returns:
        np.array: A new operator with small values set to zero.
    """
    ret = np.array(array)

    if np.isrealobj(ret):
        ret[abs(ret) < epsilon] = 0.0
    else:
        ret.real[abs(ret.real) < epsilon] = 0.0
        ret.imag[abs(ret.imag) < epsilon] = 0.0
    return ret

def is_diagonal(*mat):
    """
    Check if a matrix is diagonal.

    Parameters
    ----------
    *mat (numpy.ndarray): The matrix to check. One or more arrays

    Returns
    -------
    bool: True if the matrix is diagonal, False otherwise.
    
    """
    # Check if all off-diagonal elements are zero
    if len(mat) > 1:
        return [np.allclose(m, np.diag(np.diagonal(m))) if issquare(m) else False for m in mat]
    elif issquare(mat[0]):
        return np.allclose(mat[0], np.diag(np.diagonal(mat[0])))
    else:
        return False
    
#%% Functions from "a_variousfunctions.py" by @autor: brysa

def ishermitian(*mat):
    '''
    Check if a matrix is hermitian

    Parameters
    ----------
    *mat : np.array
        One or more arrays

    Returns
    -------
    bool or list of bools
        True if the matrix is hermitian.

    '''
    if len(mat) > 1:
        return [np.allclose(m, m.T.conj()) if issquare(m) else False for m in mat]
    elif issquare(mat[0]):
        return np.allclose(mat[0], mat[0].T.conj())
    else:
        return False
    
def issquare(mat):
    '''
    Check if a matrix is square

    Parameters
    ----------
    mat : np.array

    Returns
    -------
    bool
        True if the matrix is square.

    '''
    return mat.shape[0] == mat.shape[1] and len(mat.shape) == 2

def decompose_to_unitary(A):
    '''
    Decompose a hermitian matrix with a norm=1 to four unitary matrices
    UB, VB, UC, VC so that A = (UB+VB)/2 + 1j*(UC+VC)/2

    Parameters
    ----------
    A : np.array
        Hermitian matrix with norm=1

    Returns
    -------
    UB: np.array
        Unitary matrix.
    VB: np.array
        Unitary matrix.
    UC: np.array
        Unitary matrix.
    VC: np.array
        Unitary matrix.

    '''
    
    if not ishermitian(A) or not round(np.linalg.norm(A),6)==1:
        print('Matrix must be Hermitian and with norm=1')
        return 0, 0, 0, 0
    # Decompose a complex matrix into unitary matrices
    B = 1/2*(A+A.conj()) # A.real
    C = 1/(2j)*(A-A.conj()) # A.imag

    UB = B +1j*sp.linalg.sqrtm(np.eye(B.shape[0])-B@B)
    VB = B -1j*sp.linalg.sqrtm(np.eye(B.shape[0])-B@B)
    
    UB[np.where(abs(UB)<1e-9)] = 0
    VB[np.where(abs(VB)<1e-9)] = 0

    UC = C + 1j*sp.linalg.sqrtm(np.eye(C.shape[0])-C@C)
    VC = C - 1j*sp.linalg.sqrtm(np.eye(C.shape[0])-C@C)
    
    UC[np.where(abs(UC)<1e-9)] = 0
    VC[np.where(abs(VC)<1e-9)] = 0
    
    return UB, VB, UC, VC

def compute_control_matrix(base_mat, num_ctrl_qubits, ctrl_state=None):
    '''
    Add a control qubit to a unitary matrix. This function is taken from
    Qiskit source code.

    Parameters
    ----------
    base_mat : np.array
        matrix to be controlled.
    num_ctrl_qubits : int
        number of control qubits.
    ctrl_state : str, optional
        '0' or '1' The value of the ctrl qubit to apply the control action. The default is None.

    Returns
    -------
    full_mat : np.array
        matrix with control action applied.

    '''
    num_target = int(np.log2(base_mat.shape[0]))
    ctrl_dim = 2**num_ctrl_qubits
    ctrl_grnd = np.repeat([[1], [0]], [1, ctrl_dim - 1])
    if ctrl_state is None:
        ctrl_state = ctrl_dim - 1
    elif isinstance(ctrl_state, str):
        ctrl_state = int(ctrl_state, 2)
    
    ctrl_proj = np.diag(np.roll(ctrl_grnd, ctrl_state))
    full_mat = np.kron(np.eye(2**num_target), np.eye(ctrl_dim) - ctrl_proj) + np.kron(
        base_mat, ctrl_proj
    )
    return full_mat

def pad_to_n_qubits(A, n=None, pad_val = 1):
    '''
    Pad a matrix to size 2**n with ones on the diagonal

    Parameters
    ----------
    A : np.array
        input matrix.
    n : int, optional
        the desired size of the matrix as 2**n. If None then it pads to the 
        nearest value of 2**n for the given matrix

    Returns
    -------
    Apad : np.array
        a padded array of size 2**n.

    '''
    size_mat = len(A)
    
    if n is None:
        n = int(np.ceil(np.log2(size_mat)))

    Apad = np.pad(A,(0,int(2**n)-size_mat))
    Apad[np.arange(size_mat,int(2**n)),np.arange(size_mat,int(2**n))] = pad_val
    
    return Apad

# def isunitary(*mat):
#     '''
#     Check if a matrix is unitary

#     Parameters
#     ----------
#     *mat : np.array
#         One or more arrays

#     Returns
#     -------
#     bool or list of bools
#         True if the matrix is unitary.

#     '''
#     if len(mat) > 1:
#         return [np.allclose(np.eye(len(m)), m.dot(m.T.conj())) if issquare(m) else False for m in mat]
#     elif issquare(mat[0]):
#         return np.allclose(np.eye(len(mat[0])), mat[0].dot(mat[0].T.conj()))
#     else:
#         return False
    
# def isidentity(*mat):
#     '''
#     Check if a matrix is an identity matrix

#     Parameters
#     ----------
#     *mat : np.array
#         One or more arrays

#     Returns
#     -------
#     bool or list of bools
#         True if the matrix is an identity matrix.

#     '''
#     if len(mat) > 1:
#         return [np.allclose(np.eye(len(m)), m) if issquare(m) else False for m in mat]
#     elif issquare(mat[0]):
#         return np.allclose(np.eye(len(mat[0])), mat[0])
#     else:
#         return False 

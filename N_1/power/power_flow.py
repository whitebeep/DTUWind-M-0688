#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sun Mar 5 11:12:24 2023

@author: eduardantoligil
"""

#%% Import Module
def import_status():
    print('Module imported successfully: power/power_flow.py')
    
#%% Imports
import numpy as np
from numpy.linalg import inv

#%% DC Power Flow (DC-PF)

def PowerFlowDC(bus, B, Bf, Sbus, pvpq_index):
    """
    This PowerFlow DC function is constructed as simple as possible. It does
    not contain elements such as:
        - Gsh
        - Pbusinj
        - Pfinj
    Therefore, this could be included properly in the future.
    """

    pvpq = pvpq_index
    n_bus = B.shape[0]
    
    # Initialize voltage angle (Va_0) vector
    theta = bus['Va'].copy()
    
    # Initialize real power injection on busses
    Pbus = np.real(Sbus)
    
    # Get voltage angles (Va) on busses
    theta[pvpq] = inv(B[pvpq, :][:, pvpq]).dot(Pbus[pvpq])    
    
    # Get real power injection on busses
    P0_bus = B.dot(theta)
    
    # Get real power flow on lines
    P0_line = Bf.dot(theta)
   
    # Get voltage on busses
    V = np.ones(n_bus) * np.exp(1j*theta)
    
    success = 1
    
    # Comment if not needed
    global V_dc, success_dc, P0_line_dc, P0_bus_dc, theta_dc
    V_dc, success_dc, P0_line_dc, P0_bus_dc, theta_dc = V, success, P0_line, P0_bus, theta
    
    return V, success, P0_line, P0_bus, theta

#%% AC Fast Decoupled Load Flow (FDLF)
# There are two functions to solve the FDLF. The first one (PowerFlowFD) is more accurate than the second implementation (PowerFlowFD_V2)
# The first implementation uses Bp,Bpp and the second the imaginary part of Ybus.

def PowerFlowFD(Ybus, Bp, Bpp, Sbus, V0, pvpq_index, pq_index, max_iter = 100, err_tol = 1e-4):

    # Set iteration counter and success flag
    Q_flag = False
    success = False
    n = 0
    
    # Assign initial voltage guess (V0)
    V = V0.copy()
    Va = np.angle(V)
    Vm = abs(V)

    # Active/Reactive power indices
    dP, dQ = pvpq_index, pq_index
    
    # Calculate P and Q power mismatch
    mis = (Sbus - V * np.conj(Ybus @ V)) / Vm
    P = np.real(mis[dP])
    Q = np.imag(mis[dQ])

    # Check stop condition
    normP = np.linalg.norm(P, np.inf)
    normQ = 0
    
    if len(pq_index) > 0:
        Q_flag = True
        normQ = np.linalg.norm(Q, np.inf)
        
    if normP < err_tol and normQ < err_tol:
        success = True
    
    else: 
        # Store Bp and Bpp inverse
        Bp_inv = inv(Bp)
        if Q_flag:
            Bpp_inv = inv(Bpp)
        
    while not success and n < max_iter:
        # Iteration start
        n += 1
        
        ## Active Power (Bp)
        # Compute delta voltage angles (dVa)
        dVa = -Bp_inv @ P

        # Update voltage angles (Va)
        Va[dP] += dVa
        V = Vm * np.exp(1j * Va)
        
        # Calculate P and Q mismatch
        mis = (Sbus - V * np.conj(Ybus @ V)) / Vm
        P = np.real(mis[dP])
        Q = np.imag(mis[dQ])

        # Check stop condition
        normP = np.linalg.norm(P, np.inf)
        if Q_flag:
            normQ = np.linalg.norm(Q, np.inf)

        if normP < err_tol and normQ < err_tol:
            success = True
            break
        
        ## Reactive Power (Bpp)
        # Compute delta voltage magnitudes (dVm)
        if Q_flag:
            dVm = -Bpp_inv @ Q
            
            Vm[dQ] += dVm
            V = Vm * np.exp(1j * Va)

            mis = (Sbus - V * np.conj(Ybus @ V)) / Vm
            P = np.real(mis[dP])
            Q = np.imag(mis[dQ])

            # Check stop condition
            normP = np.linalg.norm(P, np.inf)
            normQ = np.linalg.norm(Q, np.inf)

            if normP < err_tol and normQ < err_tol:
                success = True
                break   
        else:
            pass
        
    # Comment if not needed
    global V_fdlf, success_fdlf, n_fdlf
    V_fdlf, success_fdlf, n_fdlf = V, success, n

    return V, success, n

def PowerFlowFD_V2(Ybus, Sbus, V0, pvpq_index, pq_index, max_iter = 100, err_tol = 1e-4):

    # Set iteration counter and success flag
    Q_flag = False
    success = False
    n = 0
    
    # Assign initial voltage guess (V0)
    V = V0.copy()
    Va = np.angle(V)
    Vm = abs(V)

    # Active/Reactive power indices
    dP, dQ = pvpq_index, pq_index
    
    # Calculate P and Q power mismatch
    mis = (Sbus - V * np.conj(Ybus @ V)) / Vm
    P = np.real(mis[dP])
    Q = np.imag(mis[dQ])

    # Check stop condition
    normP = np.linalg.norm(P, np.inf)
    normQ = 0
    if len(pq_index) > 0:
        Q_flag = True
        normQ = np.linalg.norm(Q, np.inf)
        
    if normP < err_tol and normQ < err_tol:
        success = True
    
    else: 
        # Compute Bp and Bpp
        Bp = np.imag(Ybus[dP, :][:, dP])
        Bpp = np.imag(Ybus[dQ, :][:, dQ])

        # Store Bp and Bpp inverse
        Bp_inv = inv(Bp)
        if Q_flag:
            Bpp_inv = inv(Bpp)
        
    while not success and n < max_iter:
        # Iteration start
        n += 1
        
        ## Active Power (Bp)
        # Compute delta voltage angles (dVa)
        dVa = -Bp_inv @ P

        # Update voltage angles (Va)
        Va[dP] += dVa
        
        # Update voltag vector
        V = Vm * np.exp(1j * Va)
        
        # Calculate P and Q mismatch
        mis = (Sbus - V * np.conj(Ybus @ V)) / Vm
        P = np.real(mis[dP])
        Q = np.imag(mis[dQ])

        # Check stop condition
        normP = np.linalg.norm(P, np.inf)
        if Q_flag:
            normQ = np.linalg.norm(Q, np.inf)

        if normP < err_tol and normQ < err_tol:
            success = True
            break
        
        ## Reactive Power (Bpp)
        if Q_flag:
            # Compute voltage magnitudes (dVm)
            dVm = -Bpp_inv @ Q
            
            # Update voltage magnitudes (Vm)
            Vm[dQ] += dVm
            
            # Update voltage vector
            V = Vm * np.exp(1j * Va)
            
            # Calculate P and Q mismatch
            mis = (Sbus - V * np.conj(Ybus @ V)) / Vm
            P = np.real(mis[dP])
            Q = np.imag(mis[dQ])

            # Check stop condition
            normP = np.linalg.norm(P, np.inf)
            normQ = np.linalg.norm(Q, np.inf)

            if normP < err_tol and normQ < err_tol:
                success = True
                break   
        else:
            pass
        
    # Comment if not needed
    global V_fdlf, success_fdlf, n_fdlf
    V_fdlf, success_fdlf, n_fdlf = V, success, n

    return V, success, n

#%% AC Newton-Raphson Power Flow

import scipy.linalg as la
# V1, success1, n1 = PowerFlowNewton(Ybus, Sbus, V0, pv_index, pq_index, max_iter=100, err_tol=1e-4)

def PowerFlowNewton(Ybus, Sbus, V0, pv_index, pq_index, max_iter=100, err_tol=1e-4):
    success = False #Initialization of status flag and iteration counter
    n = 0
    V = V0
    # print(' iteration maximum P & Q mismatch (pu)')
    # print(' --------- ---------------------------')
    # Determine mismatch between initial guess and and specified value for P and Q
    F = calculate_F(Ybus, Sbus, V, pv_index, pq_index)
    # Check if the desired tolerance is reached
    success, normF = CheckTolerance(F, err_tol)
    # Start the Newton iteration loop

    while (not success) and (n < max_iter) :
        # print(n,normF)
        n += 1 # Update counter
        # Compute derivatives and generate the Jacobian matrix
        J_dS_dVm , J_dS_dTheta = generate_Derivatives(Ybus, V)
        J = generate_Jacobian(J_dS_dVm, J_dS_dTheta, pv_index, pq_index)
        # Compute the update step
        try:
            dx = la.solve(J, F)
            # Update voltages and check if tolerance is now reached
            V = Update_Voltages(dx, V, pv_index, pq_index)
            F = calculate_F(Ybus, Sbus, V, pv_index, pq_index)
            success, normF = CheckTolerance(F, err_tol)
        except la.LinAlgWarning:
            print('Ill conditioned matrix')
            break
        except:
            print('Singular matrix')
            break
        finally:
            pass
    # print(n,normF)
    # if not success: #print out message concerning wether the power flow converged or not
    #     print('No Convergence !!!\n Stopped after %d iterations without solution...' %(n, ))

    if any(np.isnan(V)):
        print('No Convergence !!!\n NAN in voltage vector')
        success = False
    # else :
    #     print('The Newton Rapson Power Flow Converged in %d iterations!' %(n, ))
    return V, success, n

def calculate_F(Ybus, Sbus, V, pv_index, pq_index):
    Delta_S = Sbus-V*(Ybus.dot(V)).conj() # This function calculates the mismatch between the specified values of P and Q (In term of S)

    # We only use the above function for PQ and PV buses.
    Delta_P = np.real(Delta_S)
    Delta_Q = np.imag(Delta_S)

    F = np.concatenate((Delta_P[pv_index], Delta_P[pq_index], Delta_Q[pq_index]), axis = 0)

    return F

def CheckTolerance(F, err_tol):
    normF = np.linalg.norm(F,np.inf)

    if normF > err_tol:
        success = False
        # print('Not Success')
    else:
        success = True
    #     print('Success')
    # print('Highest error %.3f' %normF)
    return success, normF

def generate_Derivatives0(Ybus, V):
    J_dS_dVm = np.diag(V/np.absolute(V)).dot(np.diag((Ybus.dot(V)).conj())) + \
    np.diag(V).dot(Ybus.dot(np.diag(V/np.absolute(V))).conj())
    J_dS_dTheta = 1j*np.diag(V).dot((np.diag(Ybus.dot(V))-Ybus.dot(np.diag(V))).conj())
    return J_dS_dVm, J_dS_dTheta

def generate_Derivatives(Ybus, V):
    V = V.reshape(-1, 1)

    J_dS_dVm = (V.conj()*Ybus*(V/np.abs(V)).T).conj()
    np.fill_diagonal(J_dS_dVm, J_dS_dVm.diagonal() + (V/np.abs(V)*(Ybus@V).conj()).squeeze())
    J_dS_dTheta = -1j*V*(Ybus*V.T).conj()
    np.fill_diagonal(J_dS_dTheta, J_dS_dTheta.diagonal() + (1j*V*(Ybus@V).conj()).squeeze())

    return J_dS_dVm, J_dS_dTheta

def generate_Jacobian(J_dS_dVm, J_dS_dTheta, pv_index, pq_index):
    pvpq_ind=np.append(pv_index, pq_index)

    J_11 = np.real(J_dS_dTheta[np.ix_(pvpq_ind, pvpq_ind)])
    J_12 = np.real(J_dS_dVm[np.ix_(pvpq_ind, pq_index)])
    J_21 = np.imag(J_dS_dTheta[np.ix_(pq_index, pvpq_ind)])
    J_22 = np.imag(J_dS_dVm[np.ix_(pq_index, pq_index)])

    J = np.block([[J_11,J_12],[J_21,J_22]])

    return J

def Update_Voltages(dx, V, pv_index, pq_index):
    N1 = 0; N2 = len(pv_index) # dx[N1:N2]-ang. on the pv buses
    N3 = N2; N4 = N3 + len(pq_index) # dx[N3:N4]-ang. on the pq buses
    N5 = N4; N6 = N5 + len(pq_index) # dx[N5:N6]-mag. on the pq buses
    Theta = np.angle(V); Vm = np.absolute(V)
    if len(pv_index)>0:
        Theta[pv_index] += dx[N1:N2]
    if len(pq_index)>0:
        Theta[pq_index] += dx[N3:N4]
        Vm[pq_index] += dx[N5:N6]
    V = Vm * np.exp(1j*Theta)

    return V

#%% Post-processing functions

def ac_apparent_power(V, Ybus, Y_from, Y_to, branch):
    # Get branches indices
    br_f, br_t = branch['From'], branch['To']
    
    # Calculate S
    S_to = V[br_t]*(Y_to.dot(V)).conj()
    S_from = V[br_f]*(Y_from.dot(V)).conj()
    S_inj = V*(Ybus.dot(V)).conj()
    
    return S_to, S_from, S_inj, br_f, br_t

def DisplayResults(V, Ybus, Y_from, Y_to, branch):
    
    S_to, S_from, S_inj, br_f, br_t = ac_apparent_power(V, Ybus, Y_from, Y_to, branch)     

    dash = '=' * 60
    print()
    print(dash)
    print('|{:^58s}|'.format('Bus results'))
    print(dash)
    print('{:^6s} {:^17s} {:^17s} {:^17s}'.format('Bus', 'Voltage', 'Generation','Load'))
    print('{:^6s} {:^8s} {:^8s} {:^8s} {:^8s} {:^8s} {:^8s}'.format('#', 'Mag(pu)', 'Ang(deg)','P(pu)', 'Q(pu)','P(pu)', 'Q(pu)'))
    print('{:^6s} {:^8s} {:^8s} {:^8s} {:^8s} {:^8s} {:^8s}'.format('-'*6,'-'*8, '-'*8,'-'*8, '-'*8,'-'*8, '-'*8))

    for i in range(0 ,len(V)):
        if np.real(S_inj[i]) > 0:
            print('{:^6d} {:^8.3f} {:^8.3f} {:^8.3f} {:^8.3f} {:^8s} {:^9s}'.format(i+1, np.abs(V[i]), np.rad2deg(np.angle(V[i])), np.real(S_inj[i]), np.imag(S_inj[i]),'-','-'))
        else:
            print('{:^6d} {:^8.3f} {:^8.3f} {:^8s} {:^8s} {:^8.3f} {:^9.3f}'.format(i+1, np.abs(V[i]), np.rad2deg(np.angle(V[i])),'-','-', -np.real(S_inj[i]), -np.imag(S_inj[i])))

    print(dash)
    print('|{:^58s}|'.format('Branch Flow'))
    print(dash)
    print('{:^6s} {:<6s} {:<6s} {:^19s} {:^19s}'.format('Branch', 'From','To','From bus Injection', 'To bus Injection'))
    print('{:^6s} {:<6s} {:<6s} {:^9s} {:^9s} {:^9s} {:^9s}'.format('#','Bus','Bus','P(pu)', 'Q(pu)','P(pu)', 'Q(pu)'))


    print('{:^5s} {:^5s} {:^5s} {:^8s} {:^8s} {:^8s} {:^8s}'.format('-'*6,'-'*6, '-'*6,'-'*9, '-'*9,'-'*9, '-'*9))

    for i in range(0 ,len(br_f)):
        print('{:^6d} {:^6d} {:^6d} {:^9.3f} {:^9.3f} {:^9.3f} {:^9.3f}'.format(i+1, br_f[i]+1, br_t[i]+1, -np.real(S_from[i]), -np.imag(S_from[i]), -np.real(S_to[i]), -np.imag(S_to[i])))

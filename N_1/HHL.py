#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Apr 10 19:04:12 2023

@author: eduardantoligil
"""

#%% Controls

# General
switch_import_status = False
switch_network_init_show = False
switch_network_load_show = False
switch_ibm_account = False

# Main
switch_dc = True
switch_ac = False
switch_hhl = False

# Visualization
switch_export_table = True

#%% Imports
import numpy as np

import power.power_systems as ps
import power.power_flow as pf
# import LODF.LODF as lodf_module
import postpro.postprocessing as pp
import hhl.hhl_dc_pf as hhl_dc
import hhl.hhl_ac_pf as hhl_ac
import hhl.IBM_credentials as ibm_cred

#%% IBMQ Credentials
if switch_ibm_account:
    ibm_cred.load_account()

#%% Module import status
if switch_import_status:
    ps.import_status()
    pf.import_status()
    pp.import_status()
    hhl_dc.import_status()
    ibm_cred.import_status()

#%% power_systems.py

# Select Network Data
# network = ps.IEEE_5bus_mod
network = ps.slides_3bus

# Initialize Network
ps.init_network_data() # ps.bus, ps.branch, ps.gen

if switch_network_init_show:
    print('\nInitialize Network:')
    print('\tBus(es): ', ps.bus)
    print('\tBranch(s): ', ps.branch)
    print('\tGenerator(s): ', ps.gen)

# Load Network
ps.load_network_data(network)

if switch_network_load_show:
    print('\nLoaded Network:')
    print('\tBus(es): ', ps.bus)
    print('\tBranch(s): ', ps.branch)
    print('\tGenerator(s): ', ps.gen)

# Load PowerFlow
ps.load_pf_data('AC_NR')
ps.load_pf_data('AC_FDLF')
ps.load_pf_data('DC')

# Load PTDF and LODF matrices
ps.load_contingency_data()

#%% power_flow.py

""" DC Power Flow: Classical Solution """
if switch_dc:
    # Run DC Power Flow
    pf.PowerFlowDC(ps.bus, ps.Bdc_bus, ps.Bdc_line, ps.S_bus, ps.pvpq_index)
    # Available variables: pf.V_dc, pf.success_dc, pf.P0_line_dc, pf.P0_bus_dc, pf.theta_dc
    # print('\nDC-PF variables: pf.V_dc, pf.success_dc, pf.P0_line_dc, pf.P0_bus_dc, pf.theta_dc')
    # print('\nDC Power Flow Results:','\n V_dc: ', pf.V_dc, '\n success_dc: ', pf.success_dc, '\n Pline0_dc: ', pf.P0_line_dc, '\n Pbus0_dc: ', pf.P0_bus_dc, '\n Theta_dc: ', pf.theta_dc)

""" AC FLDF: Classical Solution """
if switch_ac:
    # Run AC Fast Decoupled Load Flow (FDLF)
    pf.PowerFlowFD(ps.Y_bus, ps.B_p, ps.B_pp, ps.S_bus, ps.V_0, ps.pvpq_index, ps.pq_index)
    # Available variables: pf.V_fdlf, success_fdlf, n_fdlf
    # print('\nAC-FDLF Classical Variables: pf.V_fdlf, success_fdlf, n_fdlf')
    # print('\nAC-FDLF Classical Results:','\nV_fdlf: ', pf.V_fdlf, '\nsuccess_fdlf: ', pf.success_fdlf, '\nn_fdlf: ', pf.n_fdlf)

#%% hhl_dc_pf.py

""" DC Power Flow: HHL Solution (Quantum) """
if switch_hhl and switch_dc:    
    # Prepare A and b (Ax=b)
    matrix_hhldc, vector_hhldc = hhl_dc.matrix_vector(ps.S_bus, ps.pvpq_index, ps.Bdc_bus)
    
    # Classical solution for comparison
    theta_cla_dc = np.linalg.solve(matrix_hhldc, vector_hhldc)
    
    # Run HHL-DC using statevector simulation
    theta_hhldc_noref = hhl_dc.solve_hhl_sv(matrix_hhldc, vector_hhldc)
    
    # print('\nQC from HHL DC-PF: ',theta_hhldc_noref)
    # print('CC1 from np.linalg.solve: ', theta_cla_dc)
    # print('CC2 from PowerFlowDC: ', pf.theta_dc[ps.pvpq_index])
    # print('Theta HHL/Classic equal?: ', np.allclose(theta_hhldc_noref, theta_cla_dc, rtol=1e-06, atol=1e-08))
    
#%% hhl_ac_pf.py

""" AC FDLF: HHL Solution (Quantum) """
if switch_hhl and switch_ac:
    V_hhlac, success_hhlac, n_hhlac, dict_hhlac = hhl_ac.HHL_FDLF(ps.Y_bus, ps.B_p, ps.B_pp, ps.S_bus, ps.V_0, ps.pvpq_index, ps.pq_index)
    # print('\nEigenvalues of Bp: ', np.linalg.eigvals(ps.B_p))
    # print('\nV_hhlac: ', V_hhlac)
    # print('\nsuccess_hhlac: ', success_hhlac)
    # print('\nn_hhlac: ', n_hhlac)
    # print('\ndict_hhlac: ', dict_hhlac)
    
    # print('\nTODO: Calculate HHL FDLF Solution')
    # print('Equal results:',np.allclose(V, V2.T, rtol=err_tol, atol=1e-08))

# %% Post-processing: LODF

""" DC Power Flow: Post-processing"""
if switch_dc:
    ##### From Classical DC Power Flow
    # Calculate post-contingency values
    PFC_dc, PF0_dc = pp.lodf_dc(pf.P0_line_dc, ps.lodf)
    # Generate .tex table
    table_pfc_dc = pp.generate_latex_table_pfc('dc', PF0_dc, PFC_dc, ps.branch, switch_export_table)
    # print('\nPF0:\n', PF0)
    # print('\nLODF:\n', LODF)
    # print('\nLODF*PF0:\n', LODF * PF0)
    # print('\nPFC:\n', PFC)
    # print('\n'+table_pfc_dc)
    # print('\nBranches:\n', branches)
    # branches = pp.branch_list(ps.branch)

if switch_hhl and switch_dc:
    ##### From HHL DC Power Flow
    # Calculate other variables
    V_hhldc, success_hhldc, P0_line_hhldc, P0_bus_hhldc, theta_hhldc = pp.postpro_hhl_dc(theta_hhldc_noref, ps.bus, ps.Bdc_bus, ps.Bdc_line, ps.pvpq_index)
    # Calculate post-contingency values
    PFC_hhldc, PF0_hhldc = pp.lodf_dc(P0_line_hhldc, ps.lodf)
    # Generate .tex table
    table_pfc_hhldc = pp.generate_latex_table_pfc('dchhl', PF0_hhldc, PFC_hhldc, ps.branch, switch_export_table)
    
    # print('\nTODO: Contingency analysis for HHL-DC')
    # print('QC',V_hhldc)
    # print('CC',pf.V_dc)
    # print('Equal results:',np.allclose(V_hhldc, pf.V_dc, rtol=1e-06, atol=1e-08))

""" AC FDLF: Post-processing """
if switch_ac:
    ##### From Classical FDLF (TODO)
    Sto_fdlf, Sfrom_fdlf, Sinj_fdlf, br_f, br_t = pf.ac_apparent_power(pf.V_fdlf, ps.Y_bus, ps.Y_from, ps.Y_to, ps.branch)
    pf.DisplayResults(pf.V_fdlf, ps.Y_bus, ps.Y_from, ps.Y_to, ps.branch)
    
if switch_hhl and switch_ac:
    ##### From HHL AC FDLF (TODO)
    print('\nTODO: Contingency analysis for HHL-FDLF')
    Sto_hhlac, Sfrom_hhlac, Sinj_hhlac, br_f, br_t = pf.ac_apparent_power(V_hhlac, ps.Y_bus, ps.Y_from, ps.Y_to, ps.branch)
    pf.DisplayResults(V_hhlac, ps.Y_bus, ps.Y_from, ps.Y_to, ps.branch)

# %% TESTS

# # power_systems.py tests
# # PS_TEST_1: START ...
# import power.old.p_sys as old

# # Old variables
# bus_old, branch_old, gen_old = old.slides_3bus()
# Y_bus_old, V_0_old, S_bus_old, pq_index_old, pv_index_old, ref_old, Y_from_old, Y_to_old = old.make_Ybus(bus_old, branch_old)
# B_p_old, B_pp_old = old.make_Bmat(bus_old, branch_old)
# Bdc_bus_old, Bdc_line_old, Pbus_inj, Pline_inj = old.make_Bdc(bus_old, branch_old)
# ptdf_old = old.ptdf(bus_old, branch_old)
# lodf_old = old.lodf(branch_old, ptdf_old)

# list_old = [bus_old, branch_old, gen_old,
#             pq_index_old, pv_index_old, ref_old,
#             V_0_old, S_bus_old,
#             Bdc_bus_old, Bdc_line_old,
#             Y_bus_old, Y_from_old, Y_to_old,
#             B_p_old, B_pp_old,
#             ptdf_old, lodf_old
#             ]

# # New variables
# bus_new = ps.bus
# branch_new = ps.branch
# gen_new = ps.gen

# pq_index_new = ps.pq_index
# pv_index_new = ps.pv_index
# ref_new = ps.ref
# pvpq_index_new = ps.pvpq_index

# V_0_new = ps.V_0
# S_bus_new = ps.S_bus
# C_ft_new = ps.C_ft

# Bdc_bus_new = ps.Bdc_bus
# Bdc_line_new = ps.Bdc_line

# Y_bus_new = ps.Y_bus
# Y_from_new = ps.Y_from
# Y_to_new = ps.Y_to

# B_p_new = ps.B_p
# B_pp_new = ps.B_pp

# ptdf_new = ps.ptdf
# lodf_new = ps.lodf

# list_new = [bus_new, branch_new, gen_new, 
#             pq_index_new, pv_index_new, ref_new,
#             V_0_new, S_bus_new,
#             Bdc_bus_new, Bdc_line_new,
#             Y_bus_new, Y_from_new, Y_to_new,
#             B_p_new, B_pp_new,
#             ptdf_new, lodf_new
#             ]

# # Comparison
# print(pq_index_new == pq_index_old,
#         pv_index_new == pv_index_old,
#         ref_new == ref_old,          
#         V_0_new == V_0_old,
#         S_bus_new == S_bus_old,
#         Bdc_bus_new == Bdc_bus_old,
#         Bdc_line_new == Bdc_line_old,
#         Y_bus_new == Y_bus_old,
#         Y_from_new == Y_from_old,
#         Y_to_new == Y_to_old,
#         B_p_new == B_p_old,
#         B_pp_new == B_pp_old,
#         ptdf_new == ptdf_old,
#         lodf_new == lodf_new
#         )
# for i in range(len(list_old)):        
#     print('\n\nOld variable: \n', list_old[i])
#     print('New variable: \n', list_old[i])
# # ... PS_TEST_1: END

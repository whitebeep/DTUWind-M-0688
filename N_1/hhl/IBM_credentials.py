#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Apr 13 18:55:42 2023

@author: eduardantoligil
"""

#%% Import Module
def import_status():
    print('Module imported successfully: hhl/IBM_credentials.py')

#%% IBMQ Credentials

from qiskit import IBMQ
from qiskit.providers.ibmq import least_busy

your_token = '993242sdfadsfadsfasdfe1223432mbm32n4b2hg34h32bjh4v23bjh4bj23hbvvnmbmnnvhv345787565456478d687s4af3ds675f38sd5g6dsg4583fg74sadg678' # Dummy token

def load_account(save=False, token=your_token):
    if save:
        IBMQ.save_account(token)
    IBMQ.load_account()
    print('\nIBMQ: account Loaded')
    
def get_real_backend(hub='ibm-q', filters = None):
    provider = IBMQ.get_provider(hub)
    # Forced filters
    
    if filters == None:
        filters = lambda x: x.configuration().n_qubits >= 3 and not x.configuration().simulator and x.status().operational == True
                                          
    available_devices = provider.backends(filters = filters)
    device = least_busy(available_devices)
    print(f'\nReal backend: running on {device.name()}')
    return device

# load_account()
# get_real_backend()

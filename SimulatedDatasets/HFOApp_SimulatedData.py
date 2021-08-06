#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 23 15:52:41 2021

@author: gz
"""

import numpy as np
from scipy.io import savemat

from mne import create_info
from mne.io import RawArray
from mne_hfo import (RMSDetector, compute_chs_hfo_rates, events_to_annotations)
from mne_hfo.simulate import simulate_hfo


# Sampling rate
srate = 2000

# length of data, in minutes
data_len = 10

# HFO frequencies
f = np.array(range(100, 222, 40))

# Number of event for each frequency
nb_evs = 20

# Number of cycles
nb_cycles = list(range(3, 10))


# simulate the testing dataset
freqs = [2.5, 6.0, 10.0, 16.0, 32.5, 67.5, 165.0,
         250.0, 425.0, 500.0, 800.0, 1500.0]

freq_ind = np.tile( np.array([0, 1, 2, 3]), (1, nb_evs))
nb_freqs = len(f)
total_evs = nb_evs * nb_freqs

freq = f[ freq_ind[0][ np.random.permutation( total_evs)] ]


for idx in range(10):
    if idx < 10:
        filename = 't00' + str(idx)
    elif idx < 100:
        filename = 't0' + str(idx)
    else:
        filename = 't' + str(idx)

    savefile = '/Users/gz/Downloads/HFOAppSimuData/' + filename + '.mat'

    # random number of cycles for each event
    c = np.random.choice( nb_cycles, total_evs, replace=True)
    max_dur = int(np.ceil(srate * max(c / freq) * 5))
    N = data_len * 60 * srate
    start_choice = np.array(range( max_dur, N - max_dur));
    
    to_continue = True
    while to_continue:
        sloc = np.random.choice( start_choice, total_evs, replace=False)
        sloc.sort()
        if min(sloc) < 3 * max_dur:
            sloc = np.random.choice( start_choice, total_evs, replace=False)
            sloc.sort()
        else:
            to_continue = False
    
    data = np.zeros(N)
    
    # generate some sinusoidal data at specified frequencies
    x = np.arange(N)
    for tmp_freq in freqs:
        # freq_amp = basic_amp / freq
        y = np.sin(2 * np.pi * tmp_freq * x / srate)
        data += y
    
    # HFO
    eloc = np.array([])
    for i in range(0, nb_evs*nb_freqs):
        sim = simulate_hfo(srate, freq[i], c[i])[0]
        ev_start = sloc[i]
        data[ev_start: ev_start + len(sim)] += sim * 10
        eloc = np.append(eloc, ev_start + len(sim))   
    
    ev_loc = np.transpose(np.array((sloc, eloc)))
    savemat(savefile, {'mat': data, 'srate':srate, 'freq':freq, 'nbcycles':c, 'ev_loc':ev_loc})
    
    print('Done ', str(idx), '/100. Good.')

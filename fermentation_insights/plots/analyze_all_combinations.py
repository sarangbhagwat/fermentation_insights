#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# !!! Package name and short description goes here.
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.
import numpy as np
import itertools
import os 
import dill

#%%
os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')

product_IDs = ['HP', 'HP_neutral', 'HP_hexanol', 'HP_neutral_hexanol',
               'TAL', 'TAL_SA',
               'succinic', 'succinic_neutral']
feedstock_IDs = ['glucose', 'sugarcane', 'corn', 'cornstover']

all_filenames = []
refinery = {i: {} for i in feedstock_IDs}
for p,f in list(itertools.product(product_IDs, feedstock_IDs)):
    all_filenames.append(p+'_'+f)
    refinery[f][p] = np.load(all_filenames[-1]+'_coefficients'+'.npy')
    
for filename in [f'TAL_SA_{str(i)}x_blp_glucose' for i in [0.1, 0.5, 2, 10]]:
    all_filenames.append(filename)
    # refinery[f][p] = np.load(all_filenames[-1])
    
coeff = {i: np.load(i+'_coefficients'+'.npy') for i in all_filenames}

ri_mpsp_f = {}
ri_mpsps = {}
for filename in all_filenames:
    with open(filename+'_RI_f_MPSP.pkl', 'rb') as f:
        ri_mpsp_f[filename] = dill.load(f)
    ri_mpsps[filename] = np.load(f'{filename}_RI_MPSP.npy')

# a*g, b*g, c*g, d*g, g, Rsq for indicators_eval

#%% Fraction of the feasible theoretical fermentation spaces where RI_MPSP>0
ri_positive_fractions = []
for filename in all_filenames:
    ris = np.log(ri_mpsps[filename])
    non_nan_ris = ris[np.where(~np.isnan(ris))]
    ris_n = len(non_nan_ris)
    ris_pos_n = len(np.where(non_nan_ris>0)[0])
    ris_neg_n = len(np.where(non_nan_ris>0)[0])
    # print(filename, ris_pos_n/ris_n)
    ri_positive_fractions.append(ris_pos_n/ris_n)
ri_positive_fractions = np.array(ri_positive_fractions)
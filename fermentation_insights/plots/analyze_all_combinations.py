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
import pandas as pd

#%%
os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')

product_IDs = [
               'TAL', 'TAL_SA',
                'HP', 'HP_neutral', 'HP_hexanol', 'HP_neutral_hexanol',
               'succinic', 'succinic_neutral']
feedstock_IDs = ['glucose', 'sugarcane', 'corn', 'cornstover']


all_pf_combinations = list(itertools.product(product_IDs, feedstock_IDs))

all_filenames = []
refinery = {i: {} for i in feedstock_IDs}
for p,f in all_pf_combinations:
    all_filenames.append(p+'_'+f)
    refinery[f][p] = np.load(all_filenames[-1]+'_coefficients'+'.npy')
    
# for filename in [f'TAL_SA_{str(i)}x_blp_glucose' for i in [0.1, 0.5, 2, 10]]:
#     all_filenames.append(filename)
#     # refinery[f][p] = np.load(all_filenames[-1])
    
coeff = {i: np.load(i+'_coefficients'+'.npy') for i in all_filenames}

ri_mpsp_f = {}
ri_mpsps = {}
for filename in all_filenames:
    with open(filename+'_RI_f_MPSP.pkl', 'rb') as f:
        ri_mpsp_f[filename] = dill.load(f)
    ri_mpsps[filename] = np.load(f'{filename}_RI_MPSP.npy')

# a*g, b*g, c*g, d*g, g, Rsq for indicators_eval

#%% Inflection yields
inflection_yields = {}
for filename in all_filenames:
    inflection_yields[filename] = np.load(f'{filename}_inflection_product_yields.npy')
    
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

#%% Contributions of d' term to MPSP across TY for all biorefineries
dcon_all = []
for filename in all_filenames:
    dcon = list(np.load(f'{filename}_dcon.npy'))
    if len(dcon[0][0])==50:
        dcon_all += dcon
    elif len(dcon[0][0])==49:
        dcon = list(dcon)
        dcon_to_add = []
        for i in dcon: 
            i = list(i)
            i_to_add = []
            for j in i:
                j_to_add = list(j)
                j_to_add.append(np.median(j_to_add))
                i_to_add.append(np.array(j_to_add))
            dcon_to_add.append(np.array(i_to_add))
        dcon_all += dcon_to_add
    else:
        breakpoint()
        
dcon_all = np.array(dcon_all)
dcon_all_non_nan = dcon_all[np.where(~np.isnan(dcon_all))]

# %% Dataframe with coeffs at baseline, 0.2x baseline, and 5.0x baseline productivity

df_coeffs_dict = dict()


for prod in ['0.2bp', '1.0bp', '5.0bp']:
    avals, bvals, cvals, dvals, gvals, Rsqs = [], [], [], [], [], []
    feeds, products = [], []
    df_coeffs_dict[prod] = dict()
    
    for (p,f) in all_pf_combinations:
        filename = p+'_'+f if prod=='1.0bp' else p+'_'+prod+'_'+f
        a, b, c, d, g, Rsq = np.load(filename+'_coefficients'+'.npy')
        avals.append(a)
        bvals.append(b)
        cvals.append(c)
        dvals.append(d)
        gvals.append(g)
        Rsqs.append(Rsq)
        feeds.append(f)
        products.append(p)
        
    df_coeffs_dict[prod]['feedstock'] = feeds
    df_coeffs_dict[prod]['product'] = products
    df_coeffs_dict[prod]["a'"] = avals
    df_coeffs_dict[prod]["b'"] = bvals
    df_coeffs_dict[prod]["c'"] = cvals
    df_coeffs_dict[prod]["d'"] = dvals
    df_coeffs_dict[prod]["recovery"] = gvals
    df_coeffs_dict[prod]['Rsq'] = Rsqs

for k, v in df_coeffs_dict.items():
    df_coeffs = pd.DataFrame.from_dict(v)
    df_coeffs.to_excel(k+'_coeffs.xlsx')

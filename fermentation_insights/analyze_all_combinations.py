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
os.chdir('C://Users//saran//Documents//Academia//repository_clones//fermentation_insights//fermentation_insights//TRY_results')

product_IDs = ['HP', 'HP_neutral', 'HP_hexanol', 'HP_neutral_hexanol',
               'TAL', 'TAL_SA',
               'succinic', 'succinic_neutral']
feedstock_IDs = ['glucose', 'sugarcane', 'corn', 'cornstover']

all_filenames = []
refinery = {i: {} for i in feedstock_IDs}
for p,f in list(itertools.product(product_IDs, feedstock_IDs)):
    all_filenames.append(p+'_'+f+'_coefficients'+'.npy')
    refinery[f][p] = np.load(all_filenames[-1])
coeff = {i: np.load(i) for i in all_filenames}

# a*g, b*g, c*g, d*g, g, Rsq for indicators_eval
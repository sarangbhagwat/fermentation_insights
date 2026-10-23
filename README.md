# fermentation_insights: Computational methods to generate generalizable, systems-scale insights into fermentation development

This module contains code to evaluate 32 biomanufacturing configurations across titer-rate-yield (TRY) combinations
and under uncertainty as well as to perform regression analyses to generate widely applicable insights; discussed in [[1]](#1).

![Example regression analyses across TRY for biorefineries producing acrylic acid via 3-HP](https://github.com/sarangbhagwat/fermentation_insights/blob/main/fermentation_insights/results/figures/HP_all_TRY_fit.png)

Getting Started
---------------

To reproduce the results, follow the installation instructions below, then directly run the script of interest, and the results will be saved as Excel files, .npy files, and .png figures in the "results" folder.

To reproduce results reported in [[1]](#1), simply do the following:

(i) Clone the fermentation_insights repository (this commit) and add the 'fermentation_insights' folder to your PYTHONPATH.

(ii) Clone the Bioindustrial-Park repository [commit 49f2e618a81b2fc1a6aae14fdd588e838b1e777f](https://github.com/BioSTEAMDevelopmentGroup/Bioindustrial-Park/tree/49f2e618a81b2fc1a6aae14fdd588e838b1e777f) and add the 'bioindustrial-park' folder to your PYTHONPATH.

(iii) Clone the BioSTEAM repository [commit e2d3942dd1076a4516efc91ae194f9e558428551](https://github.com/BioSTEAMDevelopmentGroup/biosteam/tree/e2d3942dd1076a4516efc91ae194f9e558428551)
and add the 'biosteam' folder to your PYTHONPATH.

(iv) Navigate to your local filepath bioindustrial-park/biorefineries/HP in your command prompt, and use:

```python
pip install -r requirements.txt
```

Coefficient predictor (physics-informed deep learning)
------------------------------------------------------

[`fermentation_insights/coeff_predictor`](fermentation_insights/coeff_predictor/README.md) predicts the coefficients
a', b', c', d' and the recovery g of the unifying equation, MPSP = (a' + b'/y + c'/t + d'/(y·t))/g, directly from a
specification of the feedstock, fermentation product, final product and fermentation regime, without building a process
model. It is a physics-informed deep ensemble (cost-channel decomposition, price-explicit intensities, equation layer,
physical constraints, uncertainty intervals and an applicability domain) trained on BioSTEAM simulations of 7
Bioindustrial-Park biorefineries (41 configurations). Its intervals are split-conformal, calibrated on the
leave-one-configuration-out cross-validation of the development data (so they target new configurations; coverage is
marginal over the samples of new configurations, not guaranteed for each configuration, and the quantiles are pooled
over the families: for a new configuration the MPSP intervals cover only 0.77 of the ethanol and 0.69 of the lactic
samples, the recovery intervals under-cover every TAL configuration, known or new, and new lactic ones, while the
intervals of HP and succinic acid are too wide (MPSP typically ±15-30 %); `predict()` names the under-covering keys
of each row in `intervals_undercovered`). In the cross-validation of the shipped (v4)
model, the leave-one-configuration-out MPSP error is 8.8 % (mean over the folds, each configuration counting once;
19.5 % pooled over the held-out samples; 6.9 % over the 35 folds whose product has another configuration; by family,
ranges of the folds, 1.5-3.7 % for HP and succinic acid, 8-13 % for TAL, 7-17 % for lactic acid, 17 % for the lactic
configuration without an input-identical twin in training, and 16-47 % for ethanol), but for a new feedstock (36.8 %)
or a new product family (47.0 %) the model is not better than a nearest-configuration lookup at the fold level; for
an unseen feedstock or product the intervals are not calibrated (MPSP coverage 0.76 and 0.70 of a nominal 0.90), and
such predictions are flagged. It reproduces only about two thirds of the simulated economy of scale, so scale-up
effects should not be read from it, and its response to the reagent prices, the boiler and turbogenerator efficiencies
and the income tax has the wrong sign for some configurations, so sweeps of these parameters should not be run with it
either (both checked against BioSTEAM). The numbers are those of its
[validation report](fermentation_insights/coeff_predictor/results/REPORT.md) (v4 run, trained on the data with the
D16 / D17 simulation corrections; the v1, v2 and v3 reports are kept as history). It is independent of the
publication-reproduction steps above; see its [README](fermentation_insights/coeff_predictor/README.md) for
installation, a quickstart, validation and the reproduction of its training data, and [DESIGN.md](fermentation_insights/coeff_predictor/DESIGN.md) for the design.

Documentation & Tutorials
-------------------------

- [BioSTEAM & BioIndustrial-Park](https://biosteam.readthedocs.io/en/latest/)
- [GitHub](https://docs.github.com/en)
- [Python](https://docs.python.org/3/)

## References

<a id="1">[1]</a> 
    Bhagwat, S. S.; Rao, C.V.; Zhao, H.; Singh, V.; Guest, J. S. A unifying equation for fermentation sustainability across the titer-rate-yield landscape. *Nature Communications* 2026. https://doi.org/10.1038/s41467-026-75285-1

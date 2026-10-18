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
marginal over the samples of new configurations, not guaranteed for each configuration; the recovery intervals
under-cover for some product families). In the cross-validation of the shipped (v3) model, the leave-one-configuration-out
MPSP error is 8.5 % (6.8 % over the 35 folds whose product has another configuration; by family from about 2-3 % for
HP and succinic acid to 10-12 % for TAL and lactic acid and 25 % for ethanol), but for a new
feedstock (37.0 %) or a new product family (44.8 %) the model is not better than a nearest-configuration lookup at the
fold level; for an unseen feedstock or product the intervals are not calibrated (0.71 MPSP coverage of a nominal 0.90
for both), and such predictions are flagged. The numbers are those of its
[validation report](fermentation_insights/coeff_predictor/results/REPORT.md) (v3 run; the v1 and v2 reports are kept
as history). It is independent of the publication-reproduction steps above; see
its [README](fermentation_insights/coeff_predictor/README.md) for installation, a quickstart, validation and the
reproduction of its training data, and [DESIGN.md](fermentation_insights/coeff_predictor/DESIGN.md) for the design.

Documentation & Tutorials
-------------------------

- [BioSTEAM & BioIndustrial-Park](https://biosteam.readthedocs.io/en/latest/)
- [GitHub](https://docs.github.com/en)
- [Python](https://docs.python.org/3/)

## References

<a id="1">[1]</a> 
    Bhagwat, S. S.; Rao, C.V.; Zhao, H.; Singh, V.; Guest, J. S. A unifying equation for fermentation sustainability across the titer-rate-yield landscape. *Nature Communications* 2026. https://doi.org/10.1038/s41467-026-75285-1

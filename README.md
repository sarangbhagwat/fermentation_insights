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

## References
<a id="1">[1]</a> 
    Bhagwat, S. S.; Rao, C.V.; Zhao, H.; Singh, V.; Guest, J. S. A Unifying Equation for Fermentation Sustainability across the Titer-Rate-Yield Landscape. Pre-print available online:  [https://doi.org/10.26434/chemrxiv-2025-rd5lg](https://doi.org/10.26434/chemrxiv-2025-rd5lg).

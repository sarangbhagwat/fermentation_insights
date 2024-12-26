# -*- coding: utf-8 -*-
"""
Created on Wed Dec 25 16:50:32 2024

@author: sarangbhagwat
"""

from biorefineries import HP
from biorefineries.HP.systems.sugarcane.system_sc_hexanol import HP_tea, HP_lca, R302, spec, AA, simulate_and_print, get_AA_MPSP, unit_groups_dict
from biorefineries.HP.models.sugarcane import models_sc_hexanol as models
import biosteam as bst
import os
import numpy as np

#%%
f = bst.main_flowsheet
chdir = os.chdir
ig = np.seterr(invalid='ignore')
product = AA

#%% Filepaths
HP_filepath = HP.__file__.replace('\\__init__.py', '')
HP_results_filepath = HP_filepath + '\\analyses\\results\\'

#%% Load baseline

spec.reactor.neutralization = True # !!! set neutralization here

model = models.HP_model
system = HP_sys = models.HP_sys

get_AA_MPSP()

feedstock_tag = 'sugarcane'
product_tag = 'Acrylic'

mode = '300L_FGI'

dist_filename = f'parameter-distributions_{feedstock_tag}_{product_tag}_' + mode + '.xlsx'

product_folder = 'acrylic_acid_product' if product_tag=='Acrylic' else 'HP_salt_product'

parameter_distributions_filename = HP_filepath+\
    f'\\analyses\\full\\parameter_distributions\\{product_folder}\\'+dist_filename


print(f'\n\nLoading parameter distributions ({mode}) ...')
model.parameters = ()
model.load_parameter_distributions(parameter_distributions_filename, models.namespace_dict)

# load_additional_params()
print(f'\nLoaded parameter distributions ({mode}).')

parameters = model.get_parameters()

print('\n\nLoading samples ...')
samples = model.sample(N=2000, rule='L')
model.load_samples(samples)
print('\nLoaded samples.')

# ## Change working directory to biorefineries\\HP\\analyses\\results
# chdir(HP.__file__.replace('\\__init__.py', '')+'\\analyses\\results')
# ##

model.exception_hook = 'warn'
print('\n\nSimulating baseline ...')
baseline_initial = model.metrics_at_baseline()

print(get_AA_MPSP())
get_AA_MPSP()

#%% Utils
def streams_between(sys1, sys2):
    sys1_to_sys2 = []
    for i in sys1.outs:
        for j in sys2.ins:
            if i==j: sys1_to_sys2.append(i)
    
    sys2_to_sys1 = []
    for i in sys2.outs:
        for j in sys1.ins:
            if i==j: sys2_to_sys1.append(i)
    
    return (sys1_to_sys2, sys2_to_sys1)

def flows_between(sys1, sys2):
    (sys1_to_sys2, sys2_to_sys1) = streams_between(sys1, sys2)
    flows_sys1_to_sys2, flows_sys2_to_sys1 = [], []
    for i in sys1_to_sys2:
        flows_sys1_to_sys2
        
def get_mass_flow_by_area(flow_nodes, value_factor=0.1):
    source=[]
    target=[]
    value=[]
    label=[]
    
    # label.append('feedstock production')
    
    label.append('process water center')
    
    for fn in flow_nodes:
        label.append(fn.ID.replace('_', ' '))
        
    for fn1 in flow_nodes:
        for fn2 in flow_nodes:
            fn1_to_fn2, fn2_to_fn1 = streams_between(fn1, fn2)
            for i in fn1_to_fn2:
                if i.F_mass>0.:
                    source.append(label.index(fn1.ID.replace('_', ' ')))
                    target.append(label.index(fn2.ID.replace('_', ' ')))
                    value.append(i.F_mass*value_factor)
                
            # no need to check fn2_to_fn1 since covering all permutations
    
    
    # # sugarcane = feedstock_acquisition.feeds[0]
    # sugarcane = f.sugarcane # might be more future-resistant
    # source.append(label.index('feedstock production'))
    # target.append(label.index('feedstock acquisition'))
    # value.append(sugarcane.F_mass*value_factor)
    
    
    source.append(label.index('process water center'))
    target.append(label.index('fermentation'))
    value.append(f.M304.ins[1].F_mass*value_factor)
    
    return source, target, value, label

def plot_mass_flow(flow_nodes, value_factor=0.1,
                   width=600, height=400,
                   x=None, y=None,
                   nodecolors=None,
                   linkcolors=None,
                   font='Arial',
                   font_size=18,
                   font_color='black',
                   linkalphas=0.4):
    import plotly.graph_objs as go#create links
    from plotly.offline import init_notebook_mode,  plot
    init_notebook_mode()
    source, target, value, label = get_mass_flow_by_area(flow_nodes, value_factor)
    
    if not isinstance(linkalphas, list):
        linkalphas = [linkalphas]*len(source)
    if not isinstance(linkcolors, list):
        linkcolors = [linkcolors]*len(source)
        
    if not isinstance(nodecolors, list):
        nodecolors = [nodecolors]*len(label)
    
    link = dict(source=source, target=target, value=value, 
                color=[
                    'rgba'+str(tuple(hex_to_rgba(i)+[j])) 
                    for i,j in zip(linkcolors, linkalphas)]
                )
    # color=['rgba(0,0,0, 0.5)'] * len(source))#create nodes
    node = dict(label=label, 
                pad=0, 
                thickness=8, 
                color=nodecolors,
                x=x, y=y
               )#create a sankey object
    chart = go.Sankey(link=link, node=node, 
                      arrangement="freeform")#build a figure
    fig = go.Figure(chart)
    fig.update_layout(width=width, 
                      height=height,
                      font_family=font,
                      font_color=font_color,
                      font_size=font_size,)
    plot(fig)

def hex_to_rgba(h):
    return list(int(h.lstrip('#')[i:i+2], 16) for i in (0, 2, 4))

#%% Create new systems from units and unit groups

feedstock_acquisition = bst.System.from_units('feedstock_acquisition', unit_groups_dict['feedstock acquisition'].units)
feedstock_juicing = bst.System.from_units('feedstock_juicing', unit_groups_dict['feedstock juicing'].units)

fermentation_units = unit_groups_dict['fermentation'].units
juice_evaporation_units = [f.U302, f.S301, f.F301, f.F301_P, f.M304_P]
for i in juice_evaporation_units: fermentation_units.remove(i)
fermentation = bst.System.from_units('fermentation', fermentation_units)
juice_evaporation = bst.System.from_units('juice_evaporation', juice_evaporation_units)

separation = bst.System.from_units('separation', unit_groups_dict['separation'].units)
catalytic_upgrading = bst.System.from_units('catalytic_upgrading', unit_groups_dict['upgrading'].units)

WWRR_units = unit_groups_dict['wastewater'].units
WWRR_units.remove(f.M510)
WWRR = bst.System.from_units('WWRR', WWRR_units)

CHP_units = unit_groups_dict['boiler & turbogenerator'].units
CHP_units.append(f.M510)
CHP = bst.System.from_units('CHP', CHP_units)
storage = bst.System.from_units('product_storage', [f.T620, f.T620_P])

#%% Sankey diagrams
x = np.array([0.35, 0.15, 0.28, 0.4, 0.5, 0.6, 0.68, 0.8, 0.9, 0.9])
y = [0.8, 0.15, 0.15, 0.05, 0.8, 0.4, 0.8, 0.2, 0.2, 0.8]

#%%
default_nodecolor = '#4e4e4e' # GG dark gray
default_linkcolor = '#90918e' # GG gray

yield_influenced_linkcolor = '#79bf82' # GG green
titer_influenced_linkcolor = '#60c1cf' # GG blue
yield_titer_influenced_linkcolor = '#a280b9' # GG purple

# yield_influenced_nodecolor = '#4d7e53' # GG dark green
# titer_influenced_nodecolor = '#35767f' # GG dark blue
# yield_titer_influenced_nodecolor = '#4c385a' # GG dark purple


yield_influenced_nodecolor = '#79bf82' # GG green
titer_influenced_nodecolor = '#60c1cf' # GG blue
yield_titer_influenced_nodecolor = '#a280b9' # GG purple

width = 1600

print('\n\nPlotting Sankey diagrams ...\n\n')

#%%
flow_nodes = [feedstock_acquisition,
              feedstock_juicing,
              juice_evaporation,
              fermentation,
              separation,
              catalytic_upgrading,
              WWRR,
              CHP,
              storage]

#%%
s, t, v, l = get_mass_flow_by_area(flow_nodes=flow_nodes)

#%%
yield_influenced_links = [
    ('process water center', 'fermentation'),
    ('fermentation', 'separation'),
    # ('fermentation', 'WWRR'), # sugar solution evaporation top product
    ('separation', 'catalytic upgrading'),
    ('separation', 'WWRR'),
    ('catalytic upgrading', 'WWRR'),
    ('catalytic upgrading', 'product storage'),
    ('WWRR', 'CHP'),
    ]

titer_influenced_links = [
    ('process water center', 'fermentation'),
    ('fermentation', 'separation'),
    # ('fermentation', 'WWRR'), # sugar solution evaporation top product
    ('separation', 'WWRR'),
    ]

yield_titer_influenced_links = list(set(yield_influenced_links+titer_influenced_links))

yield_titer_ratio_same_influenced_links = list(set(yield_titer_influenced_links) - set(titer_influenced_links))



#%% baseline
spec.load_specifications(spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity)
get_AA_MPSP()

plot_mass_flow(flow_nodes=flow_nodes,
                           width=width, height=800/2,
                           x=x, y=y,
                           nodecolors=default_nodecolor,
                           linkcolors=default_linkcolor)

#%% 0.5x yield
spec.load_specifications(0.5*spec.baseline_yield, spec.baseline_titer, spec.baseline_productivity)
get_AA_MPSP()

linkcolors = []
nodecolors = [default_nodecolor]*len(l)
for si, ti in zip(s,t):
    if not (l[si], l[ti]) in yield_influenced_links:
        linkcolors.append(default_linkcolor)
        
    else:
        linkcolors.append(yield_influenced_linkcolor)
        nodecolors[l.index(l[si])] = yield_influenced_nodecolor
        nodecolors[l.index(l[ti])] = yield_influenced_nodecolor
    
plot_mass_flow(flow_nodes=flow_nodes,
                           width=width, height=555/1.8,
                           x=x, y=y,
                           nodecolors=nodecolors,
                           linkcolors=linkcolors)

#%% 0.5x titer
spec.load_specifications(spec.baseline_yield, 0.5*spec.baseline_titer, spec.baseline_productivity)
get_AA_MPSP()

linkcolors = []
nodecolors = [default_nodecolor]*len(l)
for si, ti in zip(s,t):
    if not (l[si], l[ti]) in titer_influenced_links:
        linkcolors.append(default_linkcolor)
        
    else:
        linkcolors.append(titer_influenced_linkcolor)
        nodecolors[l.index(l[si])] = titer_influenced_nodecolor
        nodecolors[l.index(l[ti])] = titer_influenced_nodecolor
    
plot_mass_flow(flow_nodes=flow_nodes,
                           width=width, height=1250/2.2,
                           x=x, y=y,
                           nodecolors=nodecolors,
                           linkcolors=linkcolors)

#%% 0.5x yield and 0.5x titer and show both yield and titer influences
spec.load_specifications(0.5*spec.baseline_yield, 0.5*spec.baseline_titer, spec.baseline_productivity)
get_AA_MPSP()

linkcolors = []
nodecolors = [default_nodecolor]*len(l)

titer_influenced_nodes = [lsj for lsj, ltj in titer_influenced_links] + [ltj for lsj, ltj in titer_influenced_links]

for si, ti in zip(s,t):
    if not (l[si], l[ti]) in yield_titer_influenced_links:
        linkcolors.append(default_linkcolor)
        
    else:
        if (l[si], l[ti]) in titer_influenced_links: # influenced by titer
            linkcolors.append(yield_titer_influenced_linkcolor)
            nodecolors[l.index(l[si])] = yield_titer_influenced_nodecolor
            nodecolors[l.index(l[ti])] = yield_titer_influenced_nodecolor
        
        elif (l[si], l[ti]) in yield_influenced_links: # not influenced by titer but influenced by yield
            linkcolors.append(yield_influenced_nodecolor)
            if not l[si] in titer_influenced_nodes:
                nodecolors[l.index(l[si])] = yield_influenced_nodecolor
            if not l[ti] in titer_influenced_nodes:
                nodecolors[l.index(l[ti])] = yield_influenced_nodecolor
        
plot_mass_flow(flow_nodes=flow_nodes,
                           width=width, 
                           height=790/1.99, # 0.5, 0.5
                           # height=790/2.52, # 0.25, 0.5
                           # height=790/1.64, # 0.75, 0.5
                           x=x, y=y,
                           nodecolors=nodecolors,
                           linkcolors=linkcolors)

#%% 0.75x yield and 0.5x titer and show both yield and titer influences
# spec.load_specifications(0.75*spec.baseline_yield, 0.5*spec.baseline_titer, spec.baseline_productivity)
# get_AA_MPSP()

# linkcolors = []
# nodecolors = [default_nodecolor]*len(l)

# titer_influenced_nodes = [lsj for lsj, ltj in titer_influenced_links] + [ltj for lsj, ltj in titer_influenced_links]

# for si, ti in zip(s,t):
#     if not (l[si], l[ti]) in yield_titer_influenced_links:
#         linkcolors.append(default_linkcolor)
        
#     else:
#         if (l[si], l[ti]) in titer_influenced_links: # influenced by titer
#             linkcolors.append(yield_titer_influenced_linkcolor)
#             nodecolors[l.index(l[si])] = yield_titer_influenced_nodecolor
#             nodecolors[l.index(l[ti])] = yield_titer_influenced_nodecolor
        
#         elif (l[si], l[ti]) in yield_influenced_links: # not influenced by titer but influenced by yield
#             linkcolors.append(yield_influenced_nodecolor)
#             if not l[si] in titer_influenced_nodes:
#                 nodecolors[l.index(l[si])] = yield_influenced_nodecolor
#             if not l[ti] in titer_influenced_nodes:
#                 nodecolors[l.index(l[ti])] = yield_influenced_nodecolor
        
# plot_mass_flow(flow_nodes=flow_nodes,
#                            width=width, 

                             # height=790/1.99, # 0.5, 0.5
#                            # height=790/2.52, # 0.25, 0.5
#                            height=790/1.64, # 0.75, 0.5
#                            x=x, y=y,
#                            nodecolors=nodecolors,
#                            linkcolors=linkcolors)
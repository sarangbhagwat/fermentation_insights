#!/usr/bin/env python3
# -*- coding: utf-8 -*-
# !!! Package name and short description goes here.
# Copyright (C) 2024-, Sarang Bhagwat <sarangb2@illinois.edu>
# 
# This module is under the UIUC open-source license. See 
# github.com/BioSTEAMDevelopmentGroup/biosteam/blob/master/LICENSE.txt
# for license details.


from graphviz import Digraph
import itertools
from collections import defaultdict

#%% Utils
color = {
    'green': '#7BBD84', 
    'orange':       '#E58835', 
    'yellow':      '#F7C652', 
    'blue':        '#63C6CE', 
    'pink':        '#F8858A', 
    'gray':        '#94948C', 
    'purple':        '#734A8C', 
    'lilac':        '#D1C0E1', 
    'darker blue':        '#648496', 
    'darker green': '#607429',
    'magenta': '#A100A1',
    'lighter magenta': '#FEC1FE'
    }

light_color = {
    'yellow': '#FEEABE',
    'blue': '#B3DFE1',
    'gray': '#D9D9D4',
    'purple': '#DAC6E7', 
    'pink': '#EDC8C8',
    'green': '#BCE3C2',
    }

dark_color = {
    'pink': '#D4272D',
    'blue': '#057882',
    'gray': '#787872',
    }

#%% Settings

# Font
overall_fontname = 'Arial'
indep_fontname = 'Arial bold'
dep_fontname = 'Arial'
cluster_fontname = 'Arial bold'

overall_fontsize = 14

# Nodes
overall_nodeshape = 'box'
indep_nodeshape = 'oval'

# Edges
positive_edge_attributes = {'color': dark_color['pink']}
negative_edge_attributes = {'color': dark_color['blue']}

ambiguous_direction_edge_attributes = {'arrowhead':'None', 'color': dark_color['gray'], 'style':'dotted'}
ambiguous_color_edge_attributes = {'color': dark_color['gray'], 'style':'dashed'}

# Clusters
show_cluster_labels = True
# fermentation_fillcolor = light_color['yellow']
# separation_fillcolor = light_color['blue']
# wwrr_fillcolor = light_color['gray']
# chp_fillcolor = light_color['purple']
# systemlevel_fillcolor = 'transparent'

# fermentation_linecolor = 'transparent'
# separation_linecolor = 'transparent'
# wwrr_linecolor = 'transparent'
# chp_linecolor = 'transparent'
# systemlevel_linecolor = 'transparent'

feedstock_fillcolor = 'transparent'
fermentation_fillcolor = 'transparent'
separation_fillcolor = 'transparent'
wwrr_fillcolor = 'transparent'
chp_fillcolor = 'transparent'
systemlevel_fillcolor = 'transparent'

feedstock_linecolor = dark_color['gray']
fermentation_linecolor = dark_color['gray']
separation_linecolor = dark_color['gray']
wwrr_linecolor = dark_color['gray']
chp_linecolor = dark_color['gray']
systemlevel_linecolor = 'transparent'

cluster_label = defaultdict(lambda: None)
if show_cluster_labels:
    cluster_label.update({
        'feedstock': 'feedstock acquisition, pretreatment, saccharification',
        'fermentation': 'fermentation',
        'separation': 'separation',
        'wwrr': 'WWRR',
        'chp': 'CHP',
        'system-level': None})
# Graph
splines = 'ortho'
nodesep = 0.6

# Desirability-based box and/or text colors
default_boxcolor = 'white'
default_textcolor = 'black'
desirability_based_box_color = True
desirability_based_text_color = True

indeps = [
         'product yield',
         'product titer',
         'productivity',
         'coproduct yield',
         'cell mass yield',
         ]
desirables = [
            'recoverable energy in waste streams',
            'biogas produced',
            'heating, cooling, and power utilities produced',
            'production rate',
            ]
undesirables = [
            'dilution water',
            'time',
            'total volume of vessels',
            'water to remove',
            'organic impurities to remove',
            'product recovery',
            'water in waste streams',
            'system cost & environmental impacts',
            'sustainability indicators (e.g., MPSP, CI, FEC)',
            'f-m', 'f-e',
            's-m', 's-e',
            'w-m', 'w-e',
            'c-m', 'c-e',
            ]

box_color = defaultdict(lambda: default_boxcolor)
if desirability_based_box_color:
    box_color.update({i: 'white' for i in indeps})
    box_color.update({i: dark_color['gray'] for i in desirables})
    box_color.update({i: 'white' for i in undesirables})


text_color = defaultdict(lambda: default_textcolor)
if desirability_based_text_color:
    text_color.update({i: 'black' for i in indeps})
    text_color.update({i: 'white' for i in desirables})
    text_color.update({i: 'black' for i in undesirables})

    
#%% Graph
G = Digraph(graph_attr={
    'splines': splines,
    'nodesep': str(nodesep),
    },
            node_attr={'shape': overall_nodeshape,
                       'fontname': overall_fontname,
                       'fontsize': str(overall_fontsize),
                       'fillcolor':'transparent', 
                       'style':'filled',
                       # 'newrank': 'True',
                       })

#%% Nodes
# Feedstock nodes
with G.subgraph(name='cluster feedstock', 
                node_attr={'shape': overall_nodeshape, 
                           },
                graph_attr={'fillcolor': feedstock_fillcolor, 
                            'color': feedstock_linecolor,
                            'style': 'filled',
                            'label': cluster_label['feedstock'],
                            'fontname': cluster_fontname,
                            # 'rank': 'same',
                            }) as c:
    indep_vars = ['composition',
             'capacity',
             'price',]
    
    dep_vars = []
    
    for i in indep_vars: c.node(i, fontname=indep_fontname, fillcolor=box_color[i], style='filled', fontcolor=text_color[i], shape=indep_nodeshape)
    for i in dep_vars: c.node(i, fontname=dep_fontname, fillcolor=box_color[i], style='filled', fontcolor=text_color[i],)
    
    c.node('fs-m', 'material & energy costs', fontname=dep_fontname, fillcolor=box_color['f-m'], style='filled', fontcolor=text_color['f-m'],)
    c.node('fs-e', 'equipment purchase costs', fontname=dep_fontname, fillcolor=box_color['f-e'], style='filled', fontcolor=text_color['f-e'],)
    
    
# Fermentation nodes
with G.subgraph(name='cluster fermentation', 
                node_attr={'shape': overall_nodeshape, 
                           },
                graph_attr={'fillcolor': fermentation_fillcolor, 
                            'color': fermentation_linecolor,
                            'style': 'filled',
                            'label': cluster_label['fermentation'],
                            'fontname': cluster_fontname,
                            # 'rank': 'same',
                            }) as c:
    indep_vars = ['product yield',
             'product titer',
             'productivity',
             'coproduct yield',
             'cell mass yield',]
    dep_vars = ['dilution water',
                'time',
                'total volume of vessels',]
    
    for i in indep_vars: c.node(i, fontname=indep_fontname, fillcolor=box_color[i], style='filled', fontcolor=text_color[i], shape=indep_nodeshape)
    for i in dep_vars: c.node(i, fontname=dep_fontname, fillcolor=box_color[i], style='filled', fontcolor=text_color[i],)
    
    c.node('f-m', 'material & energy costs', fontname=dep_fontname, fillcolor=box_color['f-m'], style='filled', fontcolor=text_color['f-m'],)
    c.node('f-e', 'equipment purchase costs', fontname=dep_fontname, fillcolor=box_color['f-e'], style='filled', fontcolor=text_color['f-e'],)
    
# Separation nodes
with G.subgraph(name='cluster separation', 
                node_attr={'shape': overall_nodeshape},
                graph_attr={'fillcolor': separation_fillcolor, 
                            'color': separation_linecolor,
                            'style': 'filled',
                            'label': cluster_label['separation'],
                            'fontname': cluster_fontname,
                            # 'rank': 'same',
                            }) as c:
    indep_vars = []
    dep_vars = ['water to remove',
                'organic impurities to remove',
                'product recovery',]
    
    for i in indep_vars: c.node(i, fontname=indep_fontname, fillcolor=box_color[i], style='filled', fontcolor=text_color[i],)
    for i in dep_vars: c.node(i, fontname=dep_fontname, fillcolor=box_color[i], style='filled', fontcolor=text_color[i],)
    
    c.node('s-m', 'material & energy costs', fontname=dep_fontname, fillcolor=box_color['s-m'], style='filled', fontcolor=text_color['s-m'],)
    c.node('s-e', 'equipment purchase costs', fontname=dep_fontname, fillcolor=box_color['s-e'], style='filled', fontcolor=text_color['s-e'],)
    
# WWRR nodes
with G.subgraph(name='cluster wastewater resource recovery', 
                node_attr={'shape': overall_nodeshape},
                graph_attr={'fillcolor': wwrr_fillcolor, 
                            'color': wwrr_linecolor,
                            'style': 'filled',
                            'label': cluster_label['wwrr'],
                            'fontname': cluster_fontname,
                            # 'rank': 'same',
                            }) as c:
    indep_vars = []
    dep_vars = ['water in waste streams',
                'recoverable energy in waste streams',
                'biogas produced',]
    
    for i in indep_vars: c.node(i, fontname=indep_fontname, fillcolor=box_color[i], style='filled', fontcolor=text_color[i],)
    for i in dep_vars: c.node(i, fontname=dep_fontname, fillcolor=box_color[i], style='filled', fontcolor=text_color[i],)

    c.node('w-m', 'material & energy costs', fontname=dep_fontname, fillcolor=box_color['w-m'], style='filled', fontcolor=text_color['w-m'],)
    c.node('w-e', 'equipment purchase costs', fontname=dep_fontname, fillcolor=box_color['w-e'], style='filled', fontcolor=text_color['w-e'],)
    
# CHP nodes
with G.subgraph(name='cluster combined heat & power generation', 
                node_attr={'shape': overall_nodeshape},
                graph_attr={'fillcolor': chp_fillcolor, 
                            'color': chp_linecolor,
                            'style': 'filled',
                            'label': cluster_label['chp'],
                            'fontname': cluster_fontname,
                            # 'rank': 'same',
                            }) as c:
    indep_vars = []
    dep_vars = ['heating, cooling, and power utilities produced',]
    
    for i in indep_vars: c.node(i, fontname=indep_fontname, fillcolor=box_color[i], style='filled', fontcolor=text_color[i],)
    for i in dep_vars: c.node(i, fontname=dep_fontname, fillcolor=box_color[i], style='filled', fontcolor=text_color[i],)
    
    c.node('c-m', 'material & energy costs', fontname=dep_fontname, fillcolor=box_color['c-m'], style='filled', fontcolor=text_color['c-m'],)
    c.node('c-e', 'equipment purchase costs', fontname=dep_fontname, fillcolor=box_color['c-e'], style='filled', fontcolor=text_color['c-e'],)
    
# System-level nodes
with G.subgraph(name='cluster systemlevel', 
                node_attr={'shape': overall_nodeshape},
                graph_attr={'fillcolor': systemlevel_fillcolor,
                            'color': systemlevel_linecolor,
                            'style': 'filled',
                            'label': cluster_label['system-level'],
                            'fontname': cluster_fontname,
                            # 'rank': 'same',
                            },
                ) as c:
    indep_vars = []
    dep_vars = ['system cost & environmental impacts',
                'production rate',]
    
    for i in indep_vars: c.node(i, fontname=indep_fontname, fillcolor=box_color[i], style='filled', fontcolor=text_color[i],)
    for i in dep_vars: c.node(i, fontname=dep_fontname, fillcolor=box_color[i], style='filled', fontcolor=text_color[i],)
    c.node('sustainability indicators (e.g., MPSP, CI, FEC)', fontname='Arial bold', 
           fillcolor=box_color['sustainability indicators (e.g., MPSP, CI, FEC)'], 
           fontcolor=text_color['sustainability indicators (e.g., MPSP, CI, FEC)'],
           style='filled')

#%% Edges
positive_edges, negative_edges = [], []

# Feedstock origin edges
positive_edges.append(('price', 'fs-m'))
positive_edges.append(('capacity', 'fs-m'))
positive_edges.append(('capacity', 'fs-e'))

# Fermentation origin edges
positive_edges.append(('product yield', 'dilution water'))
positive_edges.append(('product yield', 'production rate'))
positive_edges.append(('dilution water', 'total volume of vessels'))
negative_edges.append(('product titer', 'dilution water'))
positive_edges.append(('product titer', 'time'))
negative_edges.append(('productivity', 'time'))
negative_edges.append(('time', 'total volume of vessels'))

positive_edges.append(('total volume of vessels', 'f-m'))
positive_edges.append(('total volume of vessels', 'f-e'))
# positive_edges.append(('dilution water', 'f-m'))
# positive_edges.append(('dilution water', 'f-e'))

positive_edges.append(('dilution water', 'water to remove'))
positive_edges.append(('coproduct yield', 'organic impurities to remove'))
positive_edges.append(('cell mass yield', 'organic impurities to remove'))

# Separation origin edges
positive_edges.append(('water to remove', 'water in waste streams'))
positive_edges.append(('water to remove', 's-m'))
positive_edges.append(('water to remove', 's-e'))
positive_edges.append(('organic impurities to remove', 's-m'))
positive_edges.append(('organic impurities to remove', 's-e'))

positive_edges.append(('organic impurities to remove', 'recoverable energy in waste streams'))
negative_edges.append(('product recovery', 'recoverable energy in waste streams'))
positive_edges.append(('product recovery', 'production rate'))

# WWRR origin edges
positive_edges.append(('recoverable energy in waste streams', 'biogas produced'))

positive_edges.append(('water in waste streams', 'w-m'))
positive_edges.append(('water in waste streams', 'w-e'))

positive_edges.append(('biogas produced', 'heating, cooling, and power utilities produced'))

# CHP origin edges
positive_edges.append(('heating, cooling, and power utilities produced', 'c-m'))
positive_edges.append(('heating, cooling, and power utilities produced', 'c-e'))

negative_edges.append(('heating, cooling, and power utilities produced', 'system cost & environmental impacts'))

# System-level origin edges
positive_edges.append(('system cost & environmental impacts', 'sustainability indicators (e.g., MPSP, CI, FEC)'))
negative_edges.append(('production rate', 'sustainability indicators (e.g., MPSP, CI, FEC)'))

# G.edge('product yield', 'dilution water', color='red')
# G.edge('product titer', 'dilution water', color='blue')

# Connect edges
for i in positive_edges: G.edge(i[0], i[1], _attributes=positive_edge_attributes)
for i in negative_edges: G.edge(i[0], i[1], _attributes=negative_edge_attributes)

# Connect (Process-level inputs & eq costs) -> (System-level inputs & eq costs) edges
for i in list(itertools.product(['fs-', 'f-', 's-', 'w-', 'c-'], ['m','e'])):
    G.edge(i[0]+i[1], 'system cost & environmental impacts', _attributes=positive_edge_attributes)
    
# Connect ambiguous direction and color edges
ambiguous_direction_edges = [
    
    ('product yield', 'coproduct yield'),
    ('product yield', 'cell mass yield'),
    ('coproduct yield', 'cell mass yield'),
    ]
ambiguous_color_edges = [
    ('composition', 'fs-m'),
    ('composition', 'fs-e'),
    ('composition', 'recoverable energy in waste streams'),
    ('composition', 'dilution water'),
    ('capacity', 'recoverable energy in waste streams'),
    ('composition', 'organic impurities to remove'),
    ('capacity', 'organic impurities to remove'),
    ('product titer', 'product recovery'),
    ('product yield', 'product recovery'),
    ('coproduct yield', 'product recovery'),
    ('cell mass yield', 'product recovery'),
    ]

for i in ambiguous_direction_edges: G.edge(i[0], i[1], _attributes=ambiguous_direction_edge_attributes)
for i in ambiguous_color_edges: G.edge(i[0], i[1], _attributes=ambiguous_color_edge_attributes)

#%% Plot
print(G.source)
G.render(view=True)
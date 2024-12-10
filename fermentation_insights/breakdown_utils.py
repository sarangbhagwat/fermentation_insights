# -*- coding: utf-8 -*-
"""
Created on Fri Dec  6 17:59:05 2024

@author: sarangbhagwat
"""
import numpy as np
from biosteam import System, Unit, Stream, BoilerTurbogenerator
from biosteam.units.decorators import cost
import biosteam as bst

#%% Create a system for all contributors to costs directly proportional to y for FGI equation coeff a

def get_a(ID,
            system, 
            feedstock,
            tea_class,
            broth_stream, 
            fermentation_product_ID,
            tea_defaults='CellulosicEthanolTEA',
            units_basis='spec disruption',
            boiler=None,
            spec=None,
            exclude_boiler_ic=True, # only include for coefficient 'b'
            include_ic_of_mc_units=False,
            include_ic_of_hu_units=False,):
    units, boiler = None, None
    spec_1, spec_2, spec_3 = spec.spec_1, spec.spec_2, spec.spec_3
    
    if units_basis == 'spec disruption':
        ic = dict() # installed costs
        uc = dict() # utility costs
        mc = dict() # material costs
        
        hu = dict() # heat utilities
        
        for y, t in zip([spec.baseline_yield, 0.2*spec.baseline_yield],
                        [1.8*spec.baseline_titer, 0.7*spec.baseline_titer]):
            print(y,t)
            spec.load_specifications(y, t, spec.baseline_productivity)
            for i in range(3): system.simulate()
            ic[y] = dict()
            uc[y] = dict()
            mc[y] = dict()
            hu[y] = dict()
            for i in system.units:
                if i.installed_cost is not None: ic[y][i] = i.installed_cost
                if i.utility_cost is not None: uc[y][i] = i.utility_cost
                mc[y][i] = sum(i.cost for i in i.ins)
                try:
                    duties = np.array([j.duty for j in i.heat_utilities])
                    # print(duties)
                    if np.any(np.abs(duties)>0):
                        # print(i.heat_utilities)
                        hu[y][i] = duties
                except:
                    pass
                
                
        # !!! Update IC check for 'b' coeff. Flows often too small for IC to be directly prop to y (usually lower)
        # Potentially check flows instead
        print('\nIC:')
        tot_ic = 0.
        for k in ic[spec.baseline_yield].keys():
            # print(k, ic[spec.baseline_yield][k]/1e6, ic[0.2*spec.baseline_yield][k]/1e6)
            if (not ic[0.2*spec.baseline_yield][k]==0) and\
                np.round(ic[spec.baseline_yield][k]/ic[0.2*spec.baseline_yield][k],0)==5.:
                tot_ic += ic[0.2*spec.baseline_yield][k]
                print(k, ic[0.2*spec.baseline_yield][k])
        print(tot_ic)
        
        print('\nUC:')
        tot_ucmc = 0.
        for k in uc[spec.baseline_yield].keys():
            if (not uc[0.2*spec.baseline_yield][k]==0) and\
                np.round(uc[spec.baseline_yield][k]/uc[0.2*spec.baseline_yield][k],0)==5.:
                tot_ucmc += uc[0.2*spec.baseline_yield][k]
                print(k)
                
        print('\nMC:')
        for k in mc[spec.baseline_yield].keys():
            if (not mc[0.2*spec.baseline_yield][k]==0) and\
                np.round(mc[spec.baseline_yield][k]/mc[0.2*spec.baseline_yield][k],0)==5.:
                tot_ucmc += mc[0.2*spec.baseline_yield][k]
                print(k)
                if include_ic_of_mc_units:
                    tot_ic += ic[0.2*spec.baseline_yield][k]
        print(tot_ucmc)
        
        print('\nHU:')
        all_heat_utilities = []
        for k in hu[spec.baseline_yield].keys():
            # print(hu[0.2*spec.baseline_yield])
            try:
                arr1, arr2 = hu[spec.baseline_yield][k], hu[0.2*spec.baseline_yield][k]
                if arr1.shape==arr2.shape and np.allclose(arr1, 5.*arr2, rtol=1e-1):
                    all_heat_utilities += k.heat_utilities
                    print(k)
                    if include_ic_of_hu_units:
                        tot_ic += ic[0.2*spec.baseline_yield][k]
            except:
                breakpoint()
        print(all_heat_utilities)   
        units = []
        @cost('Volume', 'Reactor',
              CE=bst.CE, cost=tot_ic, n=1., kW=0., BM=1.,)
        class NewICUnit(Unit):
            _ins_size_is_fixed = False
            _N_ins = 2
            _units = {'Volume': 'm^3'}
            def _design(self):
                self.design_results['Volume'] = 1.
            def _run(self): pass
        
        new_ic_unit = NewICUnit(ID+'_ic_unit', ins=('', ''))
        new_ic_unit.simulate()
        new_ic_unit.ins[0].imass['Water'] = 1.
        new_ic_unit.ins[0].price = tot_ucmc
        new_ic_unit.heat_utilities = [i.copy() for i in all_heat_utilities]
        units.append(new_ic_unit)
        
        
    boiler_units = [boiler] if boiler is not None else [i for i in system.facilities if isinstance(i, BoilerTurbogenerator)]
    try: assert len(boiler_units)==1.
    except: breakpoint()
    boiler = boiler_units[0]
    new_boiler = BoilerTurbogenerator(ID=ID+'_BT', ins=('', '', '', '', '', '', ''))
    # new_boiler._ins = boiler._ins
    new_boiler.ins[0].imass[fermentation_product_ID] = boiler.ins[0].imass[fermentation_product_ID]
    #!!! TODO: add final product sent to boiler, not just fermentation product
    #!!! TODO: add biogas from fermentation and final products sent to WWT
    ng = Stream(ID+'_ng', price = new_boiler.natural_gas_price)
    # if units_basis == 'spec disruption':
    #     ng-1-new_ic_unit
    #     # new_ic_unit.ins[0].price = tot_ucmc
    # else:
    ng-3-new_boiler
    new_boiler._other_units = set([new_ic_unit])
    boiler_cost_f = new_boiler._cost
    if exclude_boiler_ic:
        def no_ic_boiler_cost_f():
            boiler_cost_f()
            new_boiler.baseline_purchase_costs = dict()
            new_boiler.installed_costs = dict()
        new_boiler._cost = no_ic_boiler_cost_f 
    units += [new_boiler]
    
    yc_sys = System.from_units(ID=ID, units=units)
    yc_tea = None
    
    get_dry_flow_tpd = lambda: (feedstock.F_mass - feedstock.imass['Water'])*24./907.185
    
    if tea_defaults in ('CellulosicEthanolTEA',):
        yc_tea = tea_class(system=yc_sys, IRR=0.10, duration=(2019, 2049),
                    depreciation='MACRS7', income_tax=0.21, operating_days=180,
                    lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
                    startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
                    startup_VOCfrac=0.75, WC_over_FCI=0.05,
                    finance_interest=0.08, finance_years=10, finance_fraction=0.4,
                    OSBL_units=[new_boiler],
                    warehouse=0.04, site_development=0.09, additional_piping=0.045,
                    proratable_costs=0.10, field_expenses=0.10, construction=0.20,
                    contingency=0.10, other_indirect_costs=0.10, 
                    # labor_cost=3212962*get_dry_flow_tpd()/2205,
                    labor_cost=0.,
                    labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
                    steam_power_depreciation='MACRS20', boiler_turbogenerator=new_boiler)
    else:
        raise ValueError(f"tea_defaults must be one of {('CellulosicEthanolTEA',)}, not {tea_defaults}.")
    
    yc_sys._TEA = yc_tea
    
    new_boiler.simulate()
    # print(broth_stream.F_mass, broth_stream.imass[fermentation_product_ID].sum())
    a = yc_tea.solve_price(broth_stream)*broth_stream.F_mass/broth_stream.imass[fermentation_product_ID].sum()
    
    spec.load_specifications(spec_1, spec_2, spec_3)
    for i in range(3): system.simulate()
    
    return a, yc_sys


#%% Create a system for all contributors to fixed costs for FGI equation coeff b

#!!! NOTE: 'b' is not just the fixed equipment costs, but also the 
# fixed portion of equipment whose costs change linearly with y or t 
# (i.e., the costs of those equipment are in i*y +j)
def get_b(ID,
            system, 
            feedstock,
            tea_class,
            sugar_solution_stream,
            sugars_splitter=None,
            tea_defaults='CellulosicEthanolTEA',
            units_basis='sugars splitter',
            boiler=None,
            spec=None,
            exclude_boiler_ic=False,
            include_ic_of_hu_units=False):
    units, boiler = None, None
    spec_1, spec_2, spec_3 = spec.spec_1, spec.spec_2, spec.spec_3
    
    # if units_basis == 'sugars splitter':
        
    #     S302 = sugars_splitter
    #     units = list(S302.get_upstream_units()) + [S302]
        
    if units_basis == 'spec disruption':
        ic = dict() # installed costs
        uc = dict() # utility costs
        mc = dict() # material costs
        
        hu = dict() # heat utilities
        
        for y,t in zip([spec.baseline_yield, 0.2*spec.baseline_yield], 
                       [spec.baseline_titer, 0.5*spec.baseline_titer]):
            spec.load_specifications(y, t, spec.baseline_productivity)
            for i in range(3): system.simulate()
            ic[t] = dict()
            uc[t] = dict()
            mc[t] = dict()
            hu[t] = dict()
            for i in system.units:
                if i.installed_cost is not None: ic[t][i] = i.installed_cost
                if i.utility_cost is not None: uc[t][i] = i.utility_cost
                mc[t][i] = sum(i.cost for i in i.ins)
                try:
                    duties = np.array([j.duty for j in i.heat_utilities])
                    # print(duties)
                    if np.any(np.abs(duties)>0):
                        # print(i.heat_utilities)
                        hu[t][i] = duties
                except:
                    pass
        print('\nIC:')
        tot_ic = 0.
        for k in ic[spec.baseline_titer].keys():
            if (not ic[0.5*spec.baseline_titer][k]==0) and\
                np.round(ic[spec.baseline_titer][k]/ic[0.5*spec.baseline_titer][k],1)==1.0:
                tot_ic += ic[0.5*spec.baseline_titer][k]
                print(k)
        # print(tot_ic)
        
        print('\nUC:')
        tot_ucmc = 0.
        for k in uc[spec.baseline_titer].keys():
            if (not uc[0.5*spec.baseline_titer][k]==0) and\
                np.round(uc[spec.baseline_titer][k]/uc[0.5*spec.baseline_titer][k],1)==1.0:
                tot_ucmc += uc[0.5*spec.baseline_titer][k]
                print(k)
        
        print('\nMC:')
        for k in mc[spec.baseline_titer].keys():
            if (not mc[0.5*spec.baseline_titer][k]==0) and\
                np.round(mc[spec.baseline_titer][k]/mc[0.5*spec.baseline_titer][k],1)==1.0:
                tot_ucmc += mc[0.5*spec.baseline_titer][k]
                print(k)
        # print(tot_ucmc)
        
        print('\nHU:')
        all_heat_utilities = []
        for k in hu[spec.baseline_titer].keys():
            # print(hu[0.5*spec.baseline_titer])
            try:
                arr1, arr2 = hu[spec.baseline_titer][k], hu[0.5*spec.baseline_titer][k]
                if arr1.shape==arr2.shape and np.allclose(arr1, arr2, rtol=1e-1):
                    all_heat_utilities += k.heat_utilities
                    print(k)
                    if include_ic_of_hu_units:
                        tot_ic += ic[0.5*spec.baseline_titer][k]
            except:
                breakpoint()
        print(all_heat_utilities)   
        units = []
        @cost('Volume', 'Reactor',
              CE=bst.CE, cost=tot_ic, n=1., kW=0., BM=1.,)
        class NewICUnit(Unit):
            _ins_size_is_fixed = False
            _N_ins = 2
            _units = {'Volume': 'm^3'}
            def _design(self):
                self.design_results['Volume'] = 1.
            def _run(self): pass
        
        new_ic_unit = NewICUnit(ID+'_ic_unit', ins=('', ''))
        new_ic_unit.simulate()
        new_ic_unit.ins[0].imass['Water'] = 1.
        new_ic_unit.ins[0].price = tot_ucmc
        new_ic_unit.heat_utilities = [i.copy() for i in all_heat_utilities]
        units.append(new_ic_unit)
        
        # @cost('Volume', 'Reactor',
        #       CE=bst.CE, cost=tot_ic, n=1., kW=0., BM=1.,)
        # class NewUCMCHUUnit(Unit):
        #     _units = {'Volume': 'm^3'}
        #     def _design(self):
        #         self.design_results['Volume'] = 1.
            
        #     # @property
        #     # def _inlet_cost(self):
        #     #     return sum(i.cost for i in self.ins)
            
        # new_ucmchu_unit = NewUCMCHUUnit(ID+'_ucmchu_unit', ins='')
        # new_ucmchu_unit.ins[0].imass['Water'] = 1.
        # new_ucmchu_unit.ins[0].price = tot_ucmc
        # new_ucmchu_unit.heat_utilities = all_heat_utilities
        # units.append(new_ucmchu_unit)
        
        
    boiler_units = [boiler] if boiler is not None else [i for i in system.facilities if isinstance(i, BoilerTurbogenerator)]
    try: assert len(boiler_units)==1.
    except: breakpoint()
    boiler = boiler_units[0]
    new_boiler = BoilerTurbogenerator(ID=ID+'_BT',)
    new_boiler._ins = boiler._ins
    ng = Stream(ID+'_ng', price = new_boiler.natural_gas_price)
    # if units_basis == 'spec disruption':
    #     ng-1-new_ic_unit
    #     # new_ic_unit.ins[0].price = tot_ucmc
    # else:
    ng-3-new_boiler
    new_boiler._other_units = set([new_ic_unit])
    boiler_cost_f = new_boiler._cost
    if exclude_boiler_ic:
        def no_ic_boiler_cost_f():
            boiler_cost_f()
            new_boiler.baseline_purchase_costs = dict()
            new_boiler.installed_costs = dict()
        new_boiler._cost = no_ic_boiler_cost_f 
    units += [new_boiler]
    
    fc_sys = System.from_units(ID=ID, units=units)
    fc_tea = None
    
    get_dry_flow_tpd = lambda: (feedstock.F_mass - feedstock.imass['Water'])*24./907.185
    
    if tea_defaults in ('CellulosicEthanolTEA',):
        fc_tea = tea_class(system=fc_sys, IRR=0.10, duration=(2019, 2049),
                    depreciation='MACRS7', income_tax=0.21, operating_days=180,
                    lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
                    startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
                    startup_VOCfrac=0.75, WC_over_FCI=0.05,
                    finance_interest=0.08, finance_years=10, finance_fraction=0.4,
                    OSBL_units=[new_boiler],
                    warehouse=0.04, site_development=0.09, additional_piping=0.045,
                    proratable_costs=0.10, field_expenses=0.10, construction=0.20,
                    contingency=0.10, other_indirect_costs=0.10, 
                    labor_cost=3212962*get_dry_flow_tpd()/2205,
                    labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
                    steam_power_depreciation='MACRS20', boiler_turbogenerator=new_boiler)
    else:
        raise ValueError(f"tea_defaults must be one of {('CellulosicEthanolTEA',)}, not {tea_defaults}.")
    
    fc_sys._TEA = fc_tea
    
    new_boiler.simulate()
    
    b = fc_sys.TEA.solve_price(sugar_solution_stream)*sugar_solution_stream.F_mass/sugar_solution_stream.imass['Glucose', 'Sucrose', 'Xylose'].sum()
    
    spec.load_specifications(spec_1, spec_2, spec_3)
    for i in range(3): system.simulate()
    return b, fc_sys


#%% Create a system for all contributors to fixed costs for FGI equation coeff b

def get_c(ID,
            system, 
            feedstock,
            tea_class,
            broth_stream, 
            sugars_splitter=None,
            tea_defaults='CellulosicEthanolTEA',
            units_basis='sugars splitter',
            boiler=None,
            spec=None,
            exclude_boiler_ic=False,
            include_ic_of_hu_units=False):
    units, boiler = None, None
    spec_1, spec_2, spec_3 = spec.spec_1, spec.spec_2, spec.spec_3
    
    # if units_basis == 'sugars splitter':
        
    #     S302 = sugars_splitter
    #     units = list(S302.get_upstream_units()) + [S302]
        
    if units_basis == 'spec disruption':
        ic = dict() # installed costs
        uc = dict() # utility costs
        mc = dict() # material costs
        
        hu = dict() # heat utilities
        
        for y,t in zip([spec.baseline_yield, 0.1*spec.baseline_yield], 
                       [spec.baseline_titer, 0.5*spec.baseline_titer]):
            spec.load_specifications(y, t, spec.baseline_productivity)
            for i in range(3): system.simulate()
            ic[t] = dict()
            uc[t] = dict()
            mc[t] = dict()
            hu[t] = dict()
            for i in system.units:
                if i.installed_cost is not None: ic[t][i] = i.installed_cost
                if i.utility_cost is not None: uc[t][i] = i.utility_cost
                mc[t][i] = sum(i.cost for i in i.ins)
                try:
                    duties = np.array([j.duty for j in i.heat_utilities])
                    # print(duties)
                    if np.any(np.abs(duties)>0):
                        # print(i.heat_utilities)
                        hu[t][i] = duties
                except:
                    pass
        print('\nIC:')
        tot_ic = 0.
        for k in ic[spec.baseline_titer].keys():
            if (not ic[0.5*spec.baseline_titer][k]==0) and\
                np.round(ic[spec.baseline_titer][k]/ic[0.5*spec.baseline_titer][k],1)==5.0:
                tot_ic += ic[0.5*spec.baseline_titer][k]
                print(k)
        # print(tot_ic)
        
        print('\nUC:')
        tot_ucmc = 0.
        for k in uc[spec.baseline_titer].keys():
            if (not uc[0.5*spec.baseline_titer][k]==0) and\
                np.round(uc[spec.baseline_titer][k]/uc[0.5*spec.baseline_titer][k],1)==5.0:
                tot_ucmc += uc[0.5*spec.baseline_titer][k]
                print(k)
        
        print('\nMC:')
        for k in mc[spec.baseline_titer].keys():
            if (not mc[0.5*spec.baseline_titer][k]==0) and\
                np.round(mc[spec.baseline_titer][k]/mc[0.5*spec.baseline_titer][k],1)==5.0:
                tot_ucmc += mc[0.5*spec.baseline_titer][k]
                print(k)
        # print(tot_ucmc)
        
        print('\nHU:')
        all_heat_utilities = []
        for k in hu[spec.baseline_titer].keys():
            # print(hu[0.5*spec.baseline_titer])
            try:
                arr1, arr2 = hu[spec.baseline_titer][k], hu[0.5*spec.baseline_titer][k]
                if arr1.shape==arr2.shape and np.allclose(arr1, 5.*arr2, rtol=1e-1):
                    all_heat_utilities += k.heat_utilities
                    print(k)
                    if include_ic_of_hu_units:
                        tot_ic += ic[0.5*spec.baseline_titer][k]
            except:
                breakpoint()
        print(all_heat_utilities)   
        units = []
        @cost('Volume', 'Reactor',
              CE=bst.CE, cost=tot_ic, n=1., kW=0., BM=1.,)
        class NewICUnit(Unit):
            _ins_size_is_fixed = False
            _N_ins = 2
            _units = {'Volume': 'm^3'}
            def _design(self):
                self.design_results['Volume'] = 1.
            def _run(self): pass
        
        new_ic_unit = NewICUnit(ID+'_ic_unit', ins=('', ''))
        new_ic_unit.simulate()
        new_ic_unit.ins[0].imass['Water'] = 1.
        new_ic_unit.ins[0].price = tot_ucmc
        new_ic_unit.heat_utilities = [i.copy() for i in all_heat_utilities]
        units.append(new_ic_unit)
        
        # @cost('Volume', 'Reactor',
        #       CE=bst.CE, cost=tot_ic, n=1., kW=0., BM=1.,)
        # class NewUCMCHUUnit(Unit):
        #     _units = {'Volume': 'm^3'}
        #     def _design(self):
        #         self.design_results['Volume'] = 1.
            
        #     # @property
        #     # def _inlet_cost(self):
        #     #     return sum(i.cost for i in self.ins)
            
        # new_ucmchu_unit = NewUCMCHUUnit(ID+'_ucmchu_unit', ins='')
        # new_ucmchu_unit.ins[0].imass['Water'] = 1.
        # new_ucmchu_unit.ins[0].price = tot_ucmc
        # new_ucmchu_unit.heat_utilities = all_heat_utilities
        # units.append(new_ucmchu_unit)
        
        
    boiler_units = [boiler] if boiler is not None else [i for i in system.facilities if isinstance(i, BoilerTurbogenerator)]
    try: assert len(boiler_units)==1.
    except: breakpoint()
    boiler = boiler_units[0]
    new_boiler = BoilerTurbogenerator(ID=ID+'_BT', ins=('', '', '', '', '', '', ''))
    # new_boiler._ins = boiler._ins
    # new_boiler.ins[0].imass[fermentation_product_ID] = boiler.ins[0].imass[fermentation_product_ID]
    #!!! TODO: add final product sent to boiler, not just fermentation product
    #!!! TODO: add biogas from fermentation and final products sent to WWT
    ng = Stream(ID+'_ng', price = new_boiler.natural_gas_price)
    # if units_basis == 'spec disruption':
    #     ng-1-new_ic_unit
    #     # new_ic_unit.ins[0].price = tot_ucmc
    # else:
    ng-3-new_boiler
    new_boiler._other_units = set([new_ic_unit])
    boiler_cost_f = new_boiler._cost
    if exclude_boiler_ic:
        def no_ic_boiler_cost_f():
            boiler_cost_f()
            new_boiler.baseline_purchase_costs = dict()
            new_boiler.installed_costs = dict()
        new_boiler._cost = no_ic_boiler_cost_f 
    units += [new_boiler]
    
    y_by_t_c_sys = System.from_units(ID=ID, units=units)
    y_by_t_c_tea = None
    
    get_dry_flow_tpd = lambda: (feedstock.F_mass - feedstock.imass['Water'])*24./907.185
    
    if tea_defaults in ('CellulosicEthanolTEA',):
        y_by_t_c_tea = tea_class(system=y_by_t_c_sys, IRR=0.10, duration=(2019, 2049),
                    depreciation='MACRS7', income_tax=0.21, operating_days=180,
                    lang_factor=None, construction_schedule=(0.08, 0.60, 0.32),
                    startup_months=3, startup_FOCfrac=1, startup_salesfrac=0.5,
                    startup_VOCfrac=0.75, WC_over_FCI=0.05,
                    finance_interest=0.08, finance_years=10, finance_fraction=0.4,
                    OSBL_units=[new_boiler],
                    warehouse=0.04, site_development=0.09, additional_piping=0.045,
                    proratable_costs=0.10, field_expenses=0.10, construction=0.20,
                    contingency=0.10, other_indirect_costs=0.10, 
                    labor_cost=3212962*get_dry_flow_tpd()/2205,
                    labor_burden=0.90, property_insurance=0.007, maintenance=0.03,
                    steam_power_depreciation='MACRS20', boiler_turbogenerator=new_boiler)
    else:
        raise ValueError(f"tea_defaults must be one of {('CellulosicEthanolTEA',)}, not {tea_defaults}.")
    
    y_by_t_c_sys._TEA = y_by_t_c_tea
    
    new_boiler.simulate()
    
    c = y_by_t_c_sys.TEA.solve_price(broth_stream)*broth_stream.F_mass/broth_stream.F_vol
    
    spec.load_specifications(spec_1, spec_2, spec_3)
    for i in range(3): system.simulate()
    return c, y_by_t_c_sys

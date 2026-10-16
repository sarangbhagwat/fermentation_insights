# -*- coding: utf-8 -*-
"""
Physics layer used to build training targets (numpy only).

Unifying equation (Bhagwat et al., Nat. Commun. 2026)::

    MPSP(y, t) = (a' + b'/y + c'/t + d'/(y*t)) / g

``A, B, C, D`` (MPSP-level coefficients) are fitted by ordinary least squares
to the simulated MPSP over the feasible TY points of a sample (unweighted, as
``fermentation_insights.utils.fit_shifted_rect_hyperbola_two_param`` with
``scipy.optimize.curve_fit``, which for this linear model converges to the OLS
solution); the primed coefficients are ``a' = g*A`` etc. with ``g`` the median
recovery over the fitted points (as in ``MPSP_TRY_fit/MPSP_YT_fit_only.py``).

Cost-channel decomposition (DESIGN.md section 2)
------------------------------------------------
Per TY point ``i`` (``h`` = operating hours, h/yr)::

    R_i = MPSP_i * F_final_product_i * h                       [USD/yr]
        = sum_c E_c,i + phi * K_i + F + eps_i

Operating channels ``E_c,i`` [USD/yr] built from canonical columns:

=================  ============================================================
feedstock          F_feedstock [kg/h] * feedstock_price [USD/kg] * h
electricity        (electricity_consumption - electricity_production) [MW]
                   * 1000 [kW/MW] * electricity_price [USD/kWh] * h
natural_gas        natural_gas_flow [kg/h] * natural_gas_price [USD/kg] * h
neutralization     (base_cost + acid_cost) [USD/h] * h
materials_other    TEA material_cost [USD/yr] - feedstock - neutralization
                   (- natural_gas if natural gas is in the material pool)
utilities_other    TEA utility_cost [USD/yr] - electricity
                   (- natural_gas if natural gas is in the utility pool)
coproducts         - coproduct_sales [USD/yr] (main product sales excluded)
=================  ============================================================

Natural-gas accounting: in current BioSTEAM the boiler's natural gas is a
*stream utility* (``bst.stream_utility_prices['Fuel']``) whose cost is part of
``Unit.utility_cost`` (hence TEA ``utility_cost``) unless the stream itself
carries a price, in which case it is a priced feed and part of TEA
``material_cost``. When the harness records the NG cost per pool
(``ng_cost_utility``, ``ng_cost_material`` [USD/h], at the prices BioSTEAM
actually charges -- e.g. boiler 'Fuel' vs. dryer 'Natural gas'), each part
is subtracted from its own pool, the natural_gas channel is their sum, and
``natural_gas_price`` is the sample's cost-weighted price (pool 'split';
``int__ng`` is then the NG flow equivalent at that price). Otherwise
(legacy data: one recorded price) ``ng_pool='auto'`` decides per sample
(see :func:`detect_ng_pool`), :func:`dataset.build` makes the decision
consistent per configuration, and ``diag_ng_pool_admissible`` is False when
the chosen pool's remainder goes negative (mispriced or split NG; such
samples get ``channels_ok=False``); either way ``materials_other +
utilities_other`` is the same, only their split changes. Electricity: BioSTEAM's
``PowerUtility.cost = price * (consumption - production)`` is inside TEA
``utility_cost``, so ``utilities_other`` (the remainder) holds heat utilities,
stream utilities (ash disposal, boiler chemicals, ...) and any accounting
mismatch between the adapter's electricity columns and the power utilities.

Capital: the regressor ``K_i`` is TEA ``FCI`` [USD] when available
(``capital_basis='auto'``; this is what the DCF uses -- FCI/installed cost is
not constant across TY points because many TEAs apply site-development /
piping / warehouse factors to ISBL equipment only; on real smoke data FCI
gives lower closure and an intercept F that matches the labor part of FOC,
whereas installed cost gives a slightly negative F), else TEA
``installed_cost`` (else the sum over areas). ``phi`` is therefore the
annual capital charge per USD of FCI. The capital charge ``C_i`` (see
below) is allocated to the capital channels in proportion to the areas'
weights (:func:`capex_weights`: TDC contributions when the ISBL/OSBL split
of each area is recorded, else installed costs ``<area>__installed_cost``
[USD]): channel = ``C_i/W_i * sum(weights of its areas)``; installed cost
not covered by the areas (``I - sum_areas``) is put in ``capex_utilities``
(the harness maps unmapped unit groups to 'facilities' as well).
``fixed = F``; ``closure = R - sum(all other channels)`` so that the channels
sum exactly to the MPSP at every point.

Fixed costs and capital charge (``fixed_basis``, see :func:`decompose_sample`):
with TEA ``FOC`` recorded (default), ``F`` = intercept of FOC regressed on
``K`` (labor-type costs), ``C_i = phi_K*K_i + (FOC_i - F)`` with ``phi_K``
from the no-intercept least squares of ``R_i - sum_c E_c,i - FOC_i`` on
``K_i``; ``phi = phi_K + d FOC/dK``. Without FOC (legacy), ``phi`` and ``F``
are obtained by least squares of ``R_i - sum_c E_c,i`` on ``[K_i, 1]``
(``C_i = phi*K_i``). Both use the guards of :func:`estimate_phi_F`
(constant financial parameters within a sample); a negative ``F`` is
reported (``diag_F_negative``) and masks the channels in ``dataset``.

Each channel's MPSP contribution ``m_c,i = channel_c,i / (F_final_i*h)``
[USD/kg] is fitted with the same design matrix; by linearity of OLS the
channel coefficients sum to the total coefficients to machine precision. For
the price-explicit channels the intensities are fitted too:

* ``fs``: F_feedstock / F_final_product            [kg feedstock / kg final]
* ``el``: net electricity * 1000 / F_final_product  [kWh / kg final]
* ``ng``: natural_gas_flow / F_final_product        [kg NG / kg final]
  (pool 'split': NG cost / cost-weighted price / F_final_product)

so that ``ch__feedstock__k = feedstock_price * int__fs__k`` (exact when the
price is constant within the sample, which it is by construction).

All coefficients returned for training are in the primed (x g) basis.
"""
from __future__ import annotations

from dataclasses import dataclass, field

import numpy as np
import pandas as pd

__all__ = (
    'COEFF_KEYS', 'CHANNELS', 'OPEX_CHANNELS', 'CAPEX_CHANNEL_AREAS',
    'PRICE_CHANNELS', 'PRICE_COLUMNS', 'INTENSITY_KEYS', 'AREAS',
    'mpsp_equation', 'basis_matrix', 'FitResult', 'fit_ols', 'fit_ols_many',
    'r_squared', 'median_recovery', 'detect_ng_pool', 'ng_pool_admissibility',
    'opex_channels', 'estimate_phi_F', 'capex_weights', 'SampleDecomposition',
    'decompose_sample',
)

COEFF_KEYS = ('a', 'b', 'c', 'd')

AREAS = ('feedstock', 'conversion', 'fermentation', 'separation', 'upgrading',
         'wastewater', 'chp', 'facilities')

OPEX_CHANNELS = ('feedstock', 'electricity', 'natural_gas', 'neutralization',
                 'materials_other', 'utilities_other', 'coproducts')

#: capital channel -> canonical areas
CAPEX_CHANNEL_AREAS = {
    'capex_upstream': ('feedstock', 'conversion'),
    'capex_fermentation': ('fermentation',),
    'capex_downstream': ('separation', 'upgrading'),
    'capex_wastewater': ('wastewater',),
    'capex_utilities': ('chp', 'facilities'),
}

CHANNELS = ('feedstock', 'electricity', 'natural_gas', 'neutralization',
            'materials_other', 'utilities_other', 'capex_upstream',
            'capex_fermentation', 'capex_downstream', 'capex_wastewater',
            'capex_utilities', 'fixed', 'coproducts', 'closure')

#: price-explicit channel -> intensity key
PRICE_CHANNELS = {'feedstock': 'fs', 'electricity': 'el', 'natural_gas': 'ng'}
#: intensity key -> price column
PRICE_COLUMNS = {'fs': 'feedstock_price', 'el': 'electricity_price',
                 'ng': 'natural_gas_price'}
INTENSITY_KEYS = ('fs', 'el', 'ng')


# %% Equation and OLS fit

def mpsp_equation(y, t, a, b, c, d, g=1.0):
    """MPSP = (a + b/y + c/t + d/(y*t)) / g (works with numpy and torch)."""
    return (a + b / y + c / t + d / (y * t)) / g


def basis_matrix(y, t):
    """Design matrix ``[1, 1/y, 1/t, 1/(y*t)]`` of shape (n, 4)."""
    y = np.asarray(y, dtype=float).ravel()
    t = np.asarray(t, dtype=float).ravel()
    return np.column_stack([np.ones_like(y), 1. / y, 1. / t, 1. / (y * t)])


def r_squared(z, zhat):
    """Coefficient of determination as in the publication's ``get_Rsq``."""
    z = np.asarray(z, float)
    zhat = np.asarray(zhat, float)
    m = np.isfinite(z)
    rss = np.sum((z[m] - zhat[m]) ** 2)
    tss = np.sum((z[m] - z[m].mean()) ** 2)
    if tss <= 0:
        return 1.0 if rss <= 1e-24 * max(1., np.sum(z[m] ** 2)) else float('nan')
    return float(1. - rss / tss)


@dataclass
class FitResult:
    coef: np.ndarray        # (4,) A, B, C, D (unprimed, MPSP level)
    R2: float
    fit_mape: float         # mean |zhat - z| / |z| (fraction, not %)
    rmse: float
    n: int
    cond: float             # condition number of the column-scaled design


def _pinv(X):
    """Pseudo-inverse via column-scaled SVD (better conditioned); the result
    is a fixed linear operator, so fits of several right-hand sides with it
    are exactly additive up to floating-point rounding."""
    s = np.sqrt(np.sum(X ** 2, axis=0))
    s[s == 0] = 1.
    Xs = X / s
    P = np.linalg.pinv(Xs, rcond=1e-13)
    sv = np.linalg.svd(Xs, compute_uv=False)
    cond = float(sv[0] / sv[-1]) if sv[-1] > 0 else float('inf')
    return P / s[:, None], cond


def fit_ols_many(y, t, Z):
    """OLS coefficients (4, m) for each column of ``Z`` (n, m) with the same
    design matrix. Returns ``(coef, cond)``."""
    X = basis_matrix(y, t)
    P, cond = _pinv(X)
    Z = np.asarray(Z, float)
    return P @ Z, cond


def fit_ols(y, t, z):
    """Fit ``z = A + B/y + C/t + D/(y*t)`` by OLS."""
    y = np.asarray(y, float).ravel()
    t = np.asarray(t, float).ravel()
    z = np.asarray(z, float).ravel()
    coef, cond = fit_ols_many(y, t, z[:, None])
    coef = coef[:, 0]
    zhat = basis_matrix(y, t) @ coef
    with np.errstate(divide='ignore', invalid='ignore'):
        mape = float(np.mean(np.abs(zhat - z) / np.abs(z)))
    return FitResult(coef=coef, R2=r_squared(z, zhat), fit_mape=mape,
                     rmse=float(np.sqrt(np.mean((zhat - z) ** 2))), n=len(z),
                     cond=cond)


def median_recovery(recovery):
    r = np.asarray(recovery, float)
    r = r[np.isfinite(r)]
    return float(np.median(r)) if len(r) else float('nan')


# %% Channel construction

def _col(df, name, default=np.nan):
    """Column as float array; ``default`` where missing (column absent)."""
    if name in df.columns:
        return pd.to_numeric(df[name], errors='coerce').to_numpy(dtype=float)
    return np.full(len(df), default, dtype=float)


def ng_pool_admissibility(ng_cost_h, mat_rem_h, util_rem_h, tol=1e-6):
    """{'material': bool, 'utility': bool}: whether the remainder NG would be
    subtracted from (material_cost/h - feedstock - neutralization, or
    utility_cost/h - electricity) stays >= -tol*scale at all points."""
    ng = np.nan_to_num(np.asarray(ng_cost_h, float))
    if not len(ng):
        return {'material': True, 'utility': True}
    scale = np.nanmax(np.abs(np.concatenate([ng, mat_rem_h, util_rem_h]))) + 1e-12
    ok = lambda v: bool(np.all(np.nan_to_num(v, nan=0.) >= -tol * scale))
    return {'material': ok(mat_rem_h - ng), 'utility': ok(util_rem_h - ng)}


def detect_ng_pool(df, ng_cost_h, mat_rem_h, util_rem_h, tol=1e-6):
    """Decide whether purchased natural gas is booked in TEA material_cost
    ('material') or utility_cost ('utility'); 'none' if no NG is purchased.

    ``ng_cost_h`` NG cost [USD/h]; ``mat_rem_h`` material_cost/h - feedstock -
    neutralization [USD/h]; ``util_rem_h`` utility_cost/h - electricity.
    Rule: the CHP-area material cost decides when the per-area material
    costs are recorded and complete (a priced NG feed to the boiler appears
    in ``chp__material_cost``; a stream utility does not) -- this is direct
    evidence, whereas the sign test below is
    fooled by NG priced differently from the recorded ``natural_gas_price``
    (e.g. dryer NG at the BioSTEAM default price). Otherwise a hypothesis is
    admissible if the respective remainder minus NG stays >= -tol*scale at
    all points; if exactly one is admissible it wins, else 'utility' (current
    BioSTEAM stream-utility semantics). Per-sample decisions are made
    consistent per configuration in :func:`dataset.build`; whether the
    chosen pool is admissible is reported (``diag_ng_pool_admissible``)."""
    ng = np.nan_to_num(ng_cost_h)
    if not np.any(ng > 0):
        return 'none'
    chp_mat = _col(df, 'chp__material_cost')
    i = int(np.argmax(ng))
    # the area evidence is used only when the per-area material costs are
    # complete (sum to TEA material_cost within 1 %), i.e. feeds are attributed
    with np.errstate(divide='ignore', invalid='ignore'):
        mat_h = _col(df, 'material_cost')[i] / _col(df, 'operating_hours')[i]
    mat_areas = sum(_col(df, f'{a}__material_cost')[i] for a in AREAS)
    if np.isfinite(chp_mat[i]) and np.isfinite(mat_areas) and np.isfinite(mat_h) \
            and abs(mat_areas - mat_h) <= 0.01 * abs(mat_h):
        return 'material' if chp_mat[i] >= (1 - 1e-6) * ng[i] else 'utility'
    ok = ng_pool_admissibility(ng, mat_rem_h, util_rem_h, tol)
    if ok['material'] != ok['utility']:
        return 'material' if ok['material'] else 'utility'
    return 'utility'


def opex_channels(df, ng_pool='auto', capital_basis='auto'):
    """Annual operating-cost channels [USD/yr] per point.

    Returns ``(E, info)`` where ``E`` is a dict channel -> (n,) array for
    :data:`OPEX_CHANNELS` and ``info`` holds per-point helper arrays (R, I,
    h, F_final, intensities, area capex) and diagnostics."""
    n = len(df)
    h = _col(df, 'operating_hours')
    F_final = _col(df, 'F_final_product')
    mpsp = _col(df, 'MPSP')
    imputed = []

    def get(name, default=0.):
        v = _col(df, name)
        bad = ~np.isfinite(v)
        if bad.any():
            imputed.append(name)
            v = np.where(bad, default, v)
        return v

    F_fs = get('F_feedstock')
    p_fs = get('feedstock_price')
    e_cons = get('electricity_consumption')
    e_prod = get('electricity_production')
    p_el = get('electricity_price')
    ng_flow = get('natural_gas_flow')
    p_ng = get('natural_gas_price')
    base = get('base_cost')
    acid = get('acid_cost')
    mat = get('material_cost')
    util = get('utility_cost')
    cop = get('coproduct_sales')

    net_kW = (e_cons - e_prod) * 1000.  # MW -> kW
    fs_h = F_fs * p_fs                    # USD/h
    el_h = net_kW * p_el                  # USD/h
    neu_h = base + acid                   # USD/h
    with np.errstate(divide='ignore', invalid='ignore'):
        mat_rem_h = mat / h - fs_h - neu_h
        util_rem_h = util / h - el_h
    # NG cost recorded per pool (USD/h) at the prices BioSTEAM actually
    # charges: exact, used whenever available (ng_pool 'auto' or 'split')
    ng_u, ng_m = _col(df, 'ng_cost_utility'), _col(df, 'ng_cost_material')
    split = bool(n) and np.isfinite(ng_u).all() and np.isfinite(ng_m).all()
    q_ng_flow = ng_flow
    if split and ng_pool in ('auto', 'split'):
        pool = 'split'
        ng_h = ng_u + ng_m
        tot_flow = float(np.sum(ng_flow))
        # cost-weighted (effective) NG price of the sample; the intensity is
        # the NG flow equivalent at that price, so that channel = p * q exactly
        p_eff = float(np.sum(ng_h)) / tot_flow if tot_flow > 0 else float(np.median(p_ng))
        p_ng = np.full(n, p_eff)
        if p_eff > 0:
            q_ng_flow = ng_h / p_eff
        mat_rem_h = mat_rem_h - ng_m
        util_rem_h = util_rem_h - ng_u
        admissible = all(ng_pool_admissibility(np.zeros(n), mat_rem_h, util_rem_h).values())
    else:
        if ng_pool == 'split':
            raise ValueError('ng_pool="split" needs ng_cost_utility/ng_cost_material columns')
        ng_h = ng_flow * p_ng             # USD/h
        pool = detect_ng_pool(df, ng_h, mat_rem_h, util_rem_h) if ng_pool == 'auto' else ng_pool
        if pool not in ('material', 'utility', 'none'):
            raise ValueError(f'unknown ng_pool {ng_pool!r}')
        admissible = ng_pool_admissibility(ng_h, mat_rem_h, util_rem_h)[
            'material' if pool == 'material' else 'utility']
        if pool == 'material':
            mat_rem_h = mat_rem_h - ng_h
        else:
            util_rem_h = util_rem_h - ng_h
    # a missing TEA total must not turn its remainder into -(explicit channels)
    if 'material_cost' in imputed:
        mat_rem_h = np.zeros(n)
    if 'utility_cost' in imputed:
        util_rem_h = np.zeros(n)
    E = {
        'feedstock': fs_h * h,
        'electricity': el_h * h,
        'natural_gas': ng_h * h,
        'neutralization': neu_h * h,
        'materials_other': mat_rem_h * h,
        'utilities_other': util_rem_h * h,
        'coproducts': -cop,
    }
    # capital
    area_I = {a: _col(df, f'{a}__installed_cost') for a in AREAS}
    area_avail = np.array([np.isfinite(v).all() for v in area_I.values()])
    area_I = {a: np.nan_to_num(v) for a, v in area_I.items()}
    I_areas = sum(area_I.values())
    I_tea = _col(df, 'installed_cost')
    if np.isfinite(I_tea).all() and np.all(I_tea > 0):
        I = I_tea
    else:
        I = I_areas
        if not (np.isfinite(I_tea).all()):
            imputed.append('installed_cost')
    with np.errstate(divide='ignore', invalid='ignore'):
        q_fs = F_fs / F_final
        q_el = net_kW / F_final
        q_ng = q_ng_flow / F_final
        coverage = float(np.median(I_areas / I)) if np.all(I > 0) else float('nan')
    FCI = _col(df, 'FCI')
    fci_ok = bool(len(FCI) and np.isfinite(FCI).all() and np.all(FCI > 0))
    if capital_basis == 'auto':
        capital_basis = 'FCI' if fci_ok else 'installed_cost'
    if capital_basis == 'FCI':
        if not fci_ok:
            raise ValueError('capital_basis="FCI" but FCI is missing/non-positive')
        K = FCI
    elif capital_basis == 'installed_cost':
        K = I
    else:
        raise ValueError(f'unknown capital_basis {capital_basis!r}')
    scale = np.maximum(np.abs(mpsp * F_final * h), 1e-12)
    # consistency of the electricity columns with the per-area power data
    ea = sum(_col(df, f'{a}__electricity_consumption') - _col(df, f'{a}__electricity_production')
             for a in AREAS)
    with np.errstate(divide='ignore', invalid='ignore'):
        el_area_rel = float(np.nanmedian(np.abs(ea - (e_cons - e_prod))
                                         / np.maximum(np.abs(e_cons), 1e-12))) \
            if np.isfinite(ea).any() else np.nan
    info = dict(
        h=h, F_final=F_final, MPSP=mpsp, R=mpsp * F_final * h, I=I, K=K,
        I_areas=I_areas, area_I=area_I, q={'fs': q_fs, 'el': q_el, 'ng': q_ng},
        prices={'fs': p_fs, 'el': p_el, 'ng': p_ng},
        diagnostics=dict(
            ng_pool=pool,
            ng_pool_admissible=bool(admissible),
            capital_basis=capital_basis,
            area_coverage=coverage,
            area_data_complete=bool(area_avail.all()),
            el_area_mismatch_rel=el_area_rel,
            materials_other_min_rel=float(np.min(E['materials_other'] / scale)) if n else np.nan,
            utilities_other_min_rel=float(np.min(E['utilities_other'] / scale)) if n else np.nan,
            imputed=';'.join(sorted(set(imputed))),
        ),
    )
    return E, info


def estimate_phi_F(resid, I, phi_prior=None, min_points=3, min_rel_spread=1e-3,
                   phi_bounds=(0., np.inf), intercept=True):
    """Least squares ``resid_i = phi*I_i + F`` [USD/yr].

    Degenerate cases (fewer than ``min_points`` points, relative spread of
    ``I`` below ``min_rel_spread``, or ``phi`` outside ``phi_bounds``) fall
    back to ``phi = phi_prior`` (``F`` = mean residual) when a prior is
    given, else to ``F = 0`` (``phi`` = mean(resid)/mean(I)) for too-little
    spread, or to ``phi`` clipped into bounds with ``F`` re-estimated.
    With ``intercept=False`` the model is ``resid_i = phi*I_i`` (``F = 0``;
    identifiable from a single point; only the bounds apply).

    Returns ``(phi, F, diag)`` with diag keys ``phi_ls``, ``F_ls``,
    ``phi_flag`` ('ok' | 'prior' | 'F0' | 'clipped'), ``I_rel_spread``."""
    resid = np.asarray(resid, float)
    I = np.asarray(I, float)
    n = len(I)
    mI = float(np.mean(I)) if n else np.nan
    spread = float(np.std(I) / abs(mI)) if n and mI != 0 else 0.
    phi_ls = F_ls = np.nan
    if not intercept:
        ss = float(np.sum(I * I)) if n else 0.
        phi_ls = float(np.sum(I * resid)) / ss if ss > 0 else np.nan
        diag = dict(phi_ls=phi_ls, F_ls=0., I_rel_spread=spread)
        lo, hi = phi_bounds
        if not np.isfinite(phi_ls):
            diag['phi_flag'] = 'F0'
            return np.nan, 0., diag
        diag['phi_flag'] = 'ok' if lo <= phi_ls <= hi else 'clipped'
        return float(np.clip(phi_ls, lo, hi)), 0., diag
    if n >= max(2, min_points) and spread >= min_rel_spread:
        X = np.column_stack([I / mI, np.ones(n)])
        (p, F_ls), *_ = np.linalg.lstsq(X, resid, rcond=None)
        phi_ls = p / mI
    diag = dict(phi_ls=float(phi_ls), F_ls=float(F_ls), I_rel_spread=spread)
    lo, hi = phi_bounds
    if np.isfinite(phi_ls) and lo <= phi_ls <= hi:
        diag['phi_flag'] = 'ok'
        return float(phi_ls), float(F_ls), diag
    if phi_prior is not None and np.isfinite(phi_prior):
        phi = float(phi_prior)
        diag['phi_flag'] = 'prior'
        return phi, float(np.mean(resid - phi * I)), diag
    if not np.isfinite(phi_ls):
        diag['phi_flag'] = 'F0'
        phi = float(np.clip(np.mean(resid) / mI, lo, hi)) if n else np.nan
        return phi, float(np.mean(resid - phi * I)) if n else np.nan, diag
    phi = float(np.clip(phi_ls, lo, hi))
    diag['phi_flag'] = 'clipped'
    return phi, float(np.mean(resid - phi * I)), diag


# %% Per-sample decomposition

@dataclass
class SampleDecomposition:
    n_points: int
    g: float
    fit: FitResult                    # total fit (unprimed A..D)
    coef: np.ndarray                  # primed a', b', c', d'
    channel_coef: dict                # channel -> (4,) primed
    intensity_coef: dict              # 'fs'|'el'|'ng' -> (4,) primed
    phi: float
    F_fixed: float                    # USD/yr
    closure_rel: float
    sum_check: float                  # max |sum channels - total| / max|total|
    point_contrib: pd.DataFrame       # per point m__<channel> [USD/kg], q__*
    diagnostics: dict = field(default_factory=dict)
    prices: dict = field(default_factory=dict)

    def to_row(self) -> dict:
        row = {k: float(v) for k, v in zip(COEFF_KEYS, self.coef)}
        row.update(g=self.g, R2=self.fit.R2, n_points=self.n_points,
                   fit_mape=self.fit.fit_mape, fit_rmse=self.fit.rmse,
                   fit_cond=self.fit.cond, phi=self.phi, F_fixed=self.F_fixed,
                   closure_rel=self.closure_rel, sum_check=self.sum_check)
        for k, v in zip(COEFF_KEYS, self.fit.coef):
            row[f'{k.upper()}_mpsp'] = float(v)
        for ch in CHANNELS:
            for k, v in zip(COEFF_KEYS, self.channel_coef[ch]):
                row[f'ch__{ch}__{k}'] = float(v)
        for q in INTENSITY_KEYS:
            for k, v in zip(COEFF_KEYS, self.intensity_coef[q]):
                row[f'int__{q}__{k}'] = float(v)
        for q, col in PRICE_COLUMNS.items():
            row[col] = self.prices.get(q, np.nan)
        row.update({f'diag_{k}': v for k, v in self.diagnostics.items()})
        return row


def capex_weights(df, area_I, min_rel_spread=1e-3):
    """Per-area weights for allocating the capital charge (arrays per point)
    and diagnostics. Each area's weight is its TDC contribution
    ``I_isbl,a * (1 + f) + I_osbl,a`` when the harness recorded the OSBL part
    of each area's installed cost (``<area>__installed_cost_osbl`` [USD]),
    with ``f`` = warehouse + site development + additional piping factor
    (column ``isbl_factor``, else estimated per sample from
    ``FCI_i = c * (I_i + f * ISBL_i)`` by least squares); otherwise the
    installed costs themselves (``capex_weights='installed'``)."""
    osbl = {a: _col(df, f'{a}__installed_cost_osbl') for a in AREAS}
    out = dict(area_I)
    diag = dict(capex_weights='installed', isbl_factor=np.nan)
    if not len(df) or not all(np.isfinite(v).all() for v in osbl.values()):
        return out, diag
    isbl = {a: area_I[a] - osbl[a] for a in AREAS}
    f = float(np.nanmedian(_col(df, 'isbl_factor'))) if 'isbl_factor' in df else np.nan
    if not np.isfinite(f):
        FCI = _col(df, 'FCI')
        I_tot, S = sum(area_I.values()), sum(isbl.values())
        with np.errstate(divide='ignore', invalid='ignore'):
            frac = S / I_tot
        if len(df) >= 3 and np.isfinite(FCI).all() and np.isfinite(frac).all() \
                and np.std(frac) >= min_rel_spread * abs(np.mean(frac)):
            (c, cf), *_ = np.linalg.lstsq(np.column_stack([I_tot, S]) / np.mean(I_tot), FCI / np.mean(I_tot),
                                          rcond=None)
            f = float(cf / c) if c > 0 else np.nan
    if not np.isfinite(f) or f < 0:
        return out, diag
    out = {a: isbl[a] * (1. + f) + osbl[a] for a in AREAS}
    return out, dict(capex_weights='tdc', isbl_factor=f)


def decompose_sample(df, g=None, phi_prior=None, ng_pool='auto',
                     capital_basis='auto', fixed_basis='auto', **phi_kw):
    """Fit total coefficients and decompose them into cost channels.

    ``df`` holds the feasible points of ONE sample (canonical columns).
    ``g`` defaults to the median recovery over these points
    (:func:`dataset.build` passes the publication's definition: the median
    over all status-ok points of the sample, before inflection clipping).

    ``fixed_basis``: 'FOC' uses the TEA fixed operating cost as measured data:
    the labor-type fixed costs ``F`` are the intercept of FOC regressed on
    ``K`` (FOC = maintenance/insurance ~ K + labor), the capital-related part
    ``FOC - F`` is allocated with the capital channels, and the capital charge
    proper is the no-intercept least squares ``R - opex - FOC = phi_K * K``
    (``phi`` reported = phi_K + FOC slope; the remainder is closure).
    'intercept' is the legacy free fit ``R - opex = phi*K + F`` (it absorbs
    DCF weights != 1 and non-affine FOC into phi and F). 'auto' = 'FOC' when
    FOC is recorded at every point."""
    df = df.reset_index(drop=True)
    y = _col(df, 'y')
    t = _col(df, 't')
    mpsp = _col(df, 'MPSP')
    n = len(df)
    if g is None:
        g = median_recovery(_col(df, 'recovery'))
    fit = fit_ols(y, t, mpsp)
    coef = g * fit.coef

    E, info = opex_channels(df, ng_pool=ng_pool, capital_basis=capital_basis)
    h, F_final, R, I, K = info['h'], info['F_final'], info['R'], info['I'], info['K']
    opex_sum = sum(E.values())
    resid = R - opex_sum
    foc = _col(df, 'FOC')
    foc_ok = bool(n) and np.isfinite(foc).all()
    if fixed_basis == 'auto':
        fixed_basis = 'FOC' if foc_ok else 'intercept'
    if fixed_basis == 'FOC':
        if not foc_ok:
            raise ValueError('fixed_basis="FOC" but TEA FOC is missing')
        phi_K, _, cdiag = estimate_phi_F(resid - foc, K, intercept=False,
                                         phi_bounds=phi_kw.get('phi_bounds', (0., np.inf)))
        phi_K = phi_K if np.isfinite(phi_K) else 0.
        prior_foc = None if phi_prior is None else phi_prior - phi_K
        s_foc, F, pdiag = estimate_phi_F(foc, K, phi_prior=prior_foc, **phi_kw)
        phi = phi_K + s_foc
        charge_total = phi_K * K + (foc - F)     # USD/yr allocated to capex channels
        if pdiag['phi_flag'] == 'ok' and cdiag['phi_flag'] != 'ok':
            pdiag['phi_flag'] = cdiag['phi_flag']
        pdiag.update(phi_capital=float(phi_K), foc_slope=float(s_foc))
    elif fixed_basis == 'intercept':
        phi, F, pdiag = estimate_phi_F(resid, K, phi_prior=phi_prior, **phi_kw)
        charge_total = phi * K
    else:
        raise ValueError(f'unknown fixed_basis {fixed_basis!r}')
    pdiag['fixed_basis'] = fixed_basis
    annual = dict(E)
    # capital charge allocated to channels by the areas' weights (TDC
    # contributions when the ISBL/OSBL split is recorded, else installed
    # cost); installed cost not covered by the areas goes to capex_utilities
    w, wdiag = capex_weights(df, info['area_I'])
    pdiag.update(wdiag)
    uncovered = I - info['I_areas']
    W = sum(w.values()) + uncovered
    with np.errstate(divide='ignore', invalid='ignore'):
        charge = np.where(W > 0, charge_total / W, 0.)   # USD/yr per unit weight
    for ch, areas in CAPEX_CHANNEL_AREAS.items():
        annual[ch] = charge * sum(w[a] for a in areas)
    annual['capex_utilities'] = annual['capex_utilities'] + charge * uncovered
    annual['fixed'] = np.full(n, F)
    denom = F_final * h
    m = {}
    with np.errstate(divide='ignore', invalid='ignore'):
        for ch in CHANNELS:
            if ch != 'closure':
                m[ch] = annual[ch] / denom
    # closure computed at the MPSP level so that sum_c m_c == MPSP exactly
    m['closure'] = mpsp - sum(m[ch] for ch in CHANNELS if ch != 'closure')
    M = np.column_stack([m[ch] for ch in CHANNELS])
    Q = np.column_stack([info['q'][q] for q in INTENSITY_KEYS])
    C_all, _ = fit_ols_many(y, t, np.column_stack([M, Q]))
    C_all = g * C_all
    channel_coef = {ch: C_all[:, j] for j, ch in enumerate(CHANNELS)}
    intensity_coef = {q: C_all[:, len(CHANNELS) + j] for j, q in enumerate(INTENSITY_KEYS)}
    total_from_channels = sum(channel_coef.values())
    scale = max(np.max(np.abs(coef)), 1e-300)
    sum_check = float(np.max(np.abs(total_from_channels - coef)) / scale)
    rms = lambda v: float(np.sqrt(np.mean(np.square(v)))) if len(v) else np.nan
    closure_rel = rms(m['closure']) / rms(mpsp)

    prices = {}
    pspread = 0.
    for q in INTENSITY_KEYS:
        p = info['prices'][q]
        prices[q] = float(np.median(p)) if n else np.nan
        if n and prices[q] != 0:
            pspread = max(pspread, float(np.ptp(p) / abs(prices[q])))
    diag = dict(info['diagnostics'])
    diag.update(pdiag)
    diag['price_rel_spread'] = pspread
    diag['h'] = float(np.median(h)) if n else np.nan
    # R^2 of the (R - opex) ~ capital charge + F model
    diag['phiF_R2'] = r_squared(resid, charge_total + F)
    diag['F_negative'] = bool(F < 0)
    # cross-check: intercept of TEA FOC vs FCI (= labor-type fixed costs)
    foc, fci = _col(df, 'FOC'), _col(df, 'FCI')
    diag['F_foc_intercept'] = np.nan
    if n >= 3 and np.isfinite(foc).all() and np.isfinite(fci).all() and np.std(fci) > 0:
        X = np.column_stack([fci / np.mean(fci), np.ones(n)])
        diag['F_foc_intercept'] = float(np.linalg.lstsq(X, foc, rcond=None)[0][1])
    pc = pd.DataFrame({f'm__{ch}': m[ch] for ch in CHANNELS})
    for q in INTENSITY_KEYS:
        pc[f'q__{q}'] = info['q'][q]
    pc['installed_cost'] = I
    pc['capital_basis_value'] = K
    return SampleDecomposition(
        n_points=n, g=float(g), fit=fit, coef=coef, channel_coef=channel_coef,
        intensity_coef=intensity_coef, phi=float(phi), F_fixed=float(F),
        closure_rel=closure_rel, sum_check=sum_check, point_contrib=pc,
        diagnostics=diag, prices=prices)

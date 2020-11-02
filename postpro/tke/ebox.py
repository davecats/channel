import dns_channel.tke.read_uiuj as uiuj
import numpy as np
from scipy import integrate
from math import sqrt



def tennis(fdir, t_cnvrg, problem):

    # Calculates ebox terms with both extended Reynolds decomposition of mean field and large/small decomposition of fluctuation field.
    # Syntax:
    # phil, pl_s, pl_l, pd_s, pd_l, t_cross, eps_s, eps_l, phid = tennis(fdir, t_cnvrg, problem)
    # ------------------------------------------------------------------------------------------
    # - fdir is a string containing the path to the postprocessed files
    # - t_cnvrg is the instant of time (as of Runtimedata) at which statistical convergence is achieved
    # - problem is a string, either 'cou' or 'poi', which indicates the type of flow in analysis

    # rapid discussion of signs of returned quantities:
    # - production terms are a SOURCE FOR TKE, so they are positive if power is passed to TKE
    # - dissipation terms (laminar, deviation, turbulent large/small) here calculated are always positive
    # - sign of tcross is positive if energy flows from small to large


    # read data

    uvint,       uuint,       vvint,       wwint,       kkint,       mkint = uiuj.read_integrals(fdir)
    uvint_large, uuint_large, vvint_large, wwint_large, kkint_large, _     = uiuj.read_integrals(fdir,variant='large')
    uvint_small, uuint_small, vvint_small, wwint_small, kkint_small, _     = uiuj.read_integrals(fdir,variant='small')
    runtime, me, uv,  uu,  vv,  ww,  kk,  mk = uiuj.read(fdir, variant=None)

    _,       _,  uvS, uuS, vvS, wwS, kkS, _  = uiuj.read(fdir, variant='small')
    _,       _,  uvL, uuL, vvL, wwL, kkL, _  = uiuj.read(fdir, variant='large')

    dns_in = uiuj.read_dnsin('.//')


    # get accurate value for tauwall from runtimedata
    temp_tw_arr = (abs(runtime['Uy_bottom']) + abs(runtime['Uy_top']))/2
    temp_series = runtime['t']
    idx = temp_series[temp_series >= t_cnvrg].index[0]
    temp_tw_arr = temp_tw_arr.loc[idx:]
    tau_w = temp_tw_arr.mean() / dns_in['re'] # <--------------------


    # get ubulk from runtimedata or uw from BC
    if problem == 'cou':
        uk = abs(me['U'].iat[0]) # <--------------------
    elif problem == 'poi':
        temp_ub_arr = runtime['Ub']/2 # Ub provided by runtimedata is actually flow rate; hence divide by two
        temp_ub_arr = temp_ub_arr.loc[idx:]
        uk = temp_ub_arr.mean() # <--------------------


    # calculate power input and upi
    powt = tau_w * uk # <--------------------
    upi = powt**(1/3) # <--------------------


    # NB: produv from mke is more reliable than prod from tke, since the latter has contributions from vv and ww
    # which should be 0 but in this case only add numerical and statistical error


    # calculate alpha and beta
    uvpow = uv['var'].values/(upi**2)
    beta = integrate.simps(np.square(uvpow), uv['y'].values)/2
    if problem == 'cou':
        alpha = integrate.simps(-uvpow, uv['y'].values)/2
        rephi = dns_in['re']*upi*(alpha**2)
    elif problem == 'poi':
        integrand = np.multiply(-uvpow, (1-uv['y'].values))
        alpha = integrate.simps(integrand, uv['y'].values)/2
        rephi = 3*dns_in['re']*upi*(alpha**2)


    # hence calculate phid
    phil = 1 - rephi/2*(sqrt(1 + 4/rephi) - 1) # <--------------------
    phid = dns_in['re']*upi*beta - rephi # <--------------------


    # now read terms from integrals of large/small fluctuations
    eps_s = -kkint_small['psdiss'] / powt / 2 # <--------------------
    eps_l = -kkint_large['psdiss'] / powt / 2 # <--------------------
    t_cross = (kkint_large['tcross'] - kkint_small['tcross'])/2 / powt / 2 # <--------------------
                                                                           # sign: positive if power to large


    # now production terms are missing; for this, in addition to uv,
    # profiles of du_l/dy and du_delta/dy are needed
    dul = np.empty_like(uvpow)
    if problem == 'cou':
        dul[:] = uk/upi
    elif problem == 'poi':
        dul[:] = -3 * uk/upi * me['y'].values
    dud = me['Uy'].values/upi - dul

    # shortcuts for uv profiles
    uvsp = uvS['var']/(upi**2)
    uvlp = uvL['var']/(upi**2)

    # small production
    pl_s = -integrate.simps(np.multiply(dul,uvsp), uv['y'].values)/2 # <--------------------
    pd_s = -integrate.simps(np.multiply(dud,uvsp), uv['y'].values)/2 # <--------------------

    # large production
    pl_l = -integrate.simps(np.multiply(dul,uvlp), uv['y'].values)/2 # <--------------------
    pd_l = -integrate.simps(np.multiply(dud,uvlp), uv['y'].values)/2 # <--------------------


    return phil, pl_s, pl_l, pd_s, pd_l, t_cross, eps_s, eps_l, phid



def properties(fdir, t_cnvrg, problem):

    uvint,       uuint,       vvint,       wwint,       kkint,       mkint = uiuj.read_integrals(fdir)
    uvint_large, uuint_large, vvint_large, wwint_large, kkint_large, _     = uiuj.read_integrals(fdir,variant='large')
    uvint_small, uuint_small, vvint_small, wwint_small, kkint_small, _     = uiuj.read_integrals(fdir,variant='small')
    runtime, me, uv,  uu,  vv,  ww,  kk,  mk = uiuj.read(fdir, variant=None)

    _,       _,  uvS, uuS, vvS, wwS, kkS, _  = uiuj.read(fdir, variant='small')
    _,       _,  uvL, uuL, vvL, wwL, kkL, _  = uiuj.read(fdir, variant='large')

    dns_in = uiuj.read_dnsin('.//')

    # get accurate value for tauwall from runtimedata
    temp_tw_arr = (abs(runtime['Uy_bottom']) + abs(runtime['Uy_top']))/2
    temp_series = runtime['t']
    idx = temp_series[temp_series >= t_cnvrg].index[0]
    temp_tw_arr = temp_tw_arr.loc[idx:]
    tau_w = temp_tw_arr.mean() / dns_in['re'] # <--------------------


    # get ubulk from runtimedata or uw from BC
    if problem == 'cou':
        uk = abs(me['U'].iat[0]) # <--------------------
    elif problem == 'poi':
        temp_ub_arr = runtime['Ub']/2 # Ub provided by runtimedata is actually flow rate; hence divide by two
        temp_ub_arr = temp_ub_arr.loc[idx:]
        uk = temp_ub_arr.mean() # <--------------------


    # calculate power input and upi
    powt = tau_w * uk # <--------------------
    upi = powt**(1/3) # <--------------------

    # calculate alpha and beta
    uvpow = uv['var'].values/(upi**2)
    beta = integrate.simps(np.square(uvpow), uv['y'].values)/2
    if problem == 'cou':
        alpha = integrate.simps(-uvpow, uv['y'].values)/2
        rephi = dns_in['re']*upi*(alpha**2)
    elif problem == 'poi':
        integrand = np.multiply(-uvpow, (1-uv['y'].values))
        alpha = integrate.simps(integrand, uv['y'].values)/2
        rephi = 3*dns_in['re']*upi*(alpha**2)

    return rephi, alpha, beta
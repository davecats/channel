import channel as ch
import channel.postpro.tke as uiuj
import numpy as np
from scipy import integrate
from math import sqrt



has_been_set = False
local_power = 0


def get_properties(fdir, t_cnvrg, problem, uk, tau_w):

    _, me, uv,  uu,  vv,  ww,  kk,  mk = uiuj.read(fdir, variant=None)

    dns_in = ch.read_dnsin(fdir + 'dns.in')

    # calculate power input and upi
    powt = tau_w * uk # <--------------------
    upi = powt**(1/3) # <--------------------
    repi = dns_in['re']*upi
    local_power = powt
    has_been_set = True

    # calculate alpha for realpha
    if problem == 'cou':
        barealpha = integrate.simps(-uv['var'].values, uv['y'].values)/2
    elif problem == 'poi':
        integrand = np.multiply(-uv['var'].values, (1-uv['y'].values))
        barealpha = integrate.simps(integrand, uv['y'].values)/2
    realpha = dns_in['re'] * barealpha / uk

    # calculate alpha and beta
    uvpow = uv['var'].values/(upi**2)
    beta = integrate.simps(np.square(uvpow), uv['y'].values)/2
    if problem == 'cou':
        alpha = integrate.simps(-uvpow, uv['y'].values)/2
    elif problem == 'poi':
        integrand = np.multiply(-uvpow, (1-uv['y'].values))
        alpha = integrate.simps(integrand, uv['y'].values)/2

    return realpha, repi, alpha, beta



def get_ebox(fdir, t_cnvrg, problem, uk, tau_w):

    _, _, _, _, kkint, _ = uiuj.read_integrals(fdir)
    dns_in = ch.read_dnsin(fdir + 'dns.in')

    realpha, repi, alpha, beta = get_properties(fdir, t_cnvrg, problem, uk, tau_w)
    upi = repi/dns_in['re']
    powt = upi**3
    if has_been_set:
        powt = local_power

    #phil = 1/(realpha+1)
    #pl = realpha/(realpha+1)

    if problem == 'cou':
        pl = (repi*alpha**2)/2 * (sqrt(1+4/(repi*alpha**2)) - 1)
        phid = repi*(beta-alpha**2)
    elif problem == 'poi':
        pl = (3*repi*alpha**2)/2 * (sqrt(1+4/(3*repi*alpha**2)) - 1)
        phid = repi*(beta-3*alpha**2)
    
    phil = 1-pl
    pd = - phid

    eps = -kkint['psdiss'] / powt / 2

    return phil, phid, pl, pd, eps




def tennis(fdir, t_cnvrg, problem, uk, tau_w):

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

    _, _, _, _, kkint_large, _     = uiuj.read_integrals(fdir,variant='large')
    _, _, _, _, kkint_small, _     = uiuj.read_integrals(fdir,variant='small')
    _, me, uv,  _,  _,  _,  _,  _ = uiuj.read(fdir, variant=None)

    _,       _,  uvS, _, _, _, _, _  = uiuj.read(fdir, variant='small')
    _,       _,  uvL, _, _, _, _, _  = uiuj.read(fdir, variant='large')

    dns_in = ch.read_dnsin(fdir + 'dns.in')

    # get properties
    realpha, repi, alpha, beta = get_properties(fdir, t_cnvrg, problem, uk, tau_w)
    upi = repi/dns_in['re']
    powt = upi**3
    uk_upi = repi * alpha / realpha # = uk / upi

    # hence calculate phid
    phil = 1/(1+realpha) # <--------------------
    if problem == 'cou':
        phid = repi*(beta-alpha**2)
    elif problem == 'poi':
        phid = repi*(beta-3*alpha**2) # <--------------------


    # now read terms from integrals of large/small fluctuations
    eps_s = -kkint_small['psdiss'] / powt / 2 # <--------------------
    eps_l = -kkint_large['psdiss'] / powt / 2 # <--------------------
    t_cross = (kkint_large['tcross'] - kkint_small['tcross'])/2 / powt / 2 # <--------------------
                                                                           # sign: positive if power to large


    # now production terms are missing; for this, in addition to uv,
    # profiles of du_l/dy and du_delta/dy are needed
    dul = np.empty_like(uv['var'].values)
    if problem == 'cou':
        dul[:] = uk_upi
    elif problem == 'poi':
        dul[:] = -3 * uk_upi * me['y'].values
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
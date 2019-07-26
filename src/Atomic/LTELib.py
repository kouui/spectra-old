################################################################################
# this file defines functions for
#     calculations related to LTE process
#
################################################################################

import numpy as np
import numba as nb
from .. import Constants as Cst
from ..Structure import AtomCls

def Boltzmann_distribution(_gi, _gj, _Eji, _Te):
    r"""

    calculate the population ratio between upper level j and lower level i under LTE.

    with `nb.vectorize( [nb.float64(nb.uint8, nb.uint8, nb.float64, nb.float64)])`.

    Parameters
    ----------

    _gi  : np.uint8 or array-like
        statistical weight of lower level i, [-]
    _gj  : np.uint8 or array-like
        statistical weight of upper level j, [-]
    _Eji : np.double or array-like
        the gap of level energy between upper level j and lower level i,
        :math:`E_{ji}=E_j-E_i, \quad [erg]`
    _Te  : np.double or array-like
        electron temperature, [:math:`K`]

    Returns
    -------

    _rt : np.double or array-like
        population ratio of upper level j and lower level i.
        :math:`rt = n_j / n_i`, [-]

    Notes
    -----

    The population ratio according to Boltzmann distribution [1]_.

    .. math:: \frac{n_j}{n_i} = \frac{g_j}{g_i} e^{-(E_j-E_i)/{kT}}

    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 262, 2015.
    """

    _kT = Cst.k_ * _Te
    _rt = (_gj/_gi) * np.exp(-_Eji/_kT)

    return _rt

def Saha_distribution(_gi, _gk, _chi, _ne, _Te):
    r"""

    calculate the population ratio between the ground states of two subsequent ionization stage under LTE.

    with `nb.vectorize( [nb.float64(nb.uint8,nb.uint8,nb.float64,nb.float64,nb.float64)])`

    Parameters
    ----------

    _gk: np.uint8 or array-like
        statistical weight of the ground state in ionization stage I+1, [-]
    _gi: np.uint8 or array-like
        statistical weight of the ground state in ionization stage I, [-]
    _chi: np.double or array-like
        ionization energy from ground state in ionization stage I
        to ground state in ionization stage I, [:math:`erg`]
    _ne: np.double or array-like
        electron density, [:math:`cm^{-3}`]
    _Te: np.double or array-like
        electron temperature, [:math:`K`]

    Returns
    -------
    _rt: np.double or array-like
        :math:`n_k / n_i`, :math:`n_k` the population of the ground state of ionization stage I+1
        and :math:`n_i` the population of the ground state of ionization stage I. [-]

    Notes
    -----

    The population ratio according to Saha's equation [1]_.

    a factor contains physics constants only,

    .. math:: f = 2(\frac{2 \pi m_e k}{h^2})^{3/2}

    then the population ratio is,

    .. math:: \frac{n_{0,k}}{n_{0,i}} = f \cdot \frac{g_{0,k}}{g_{0,i}} e^{\frac{-\chi_{ki}}{kT}}  T^{3/2} / n_e

    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 94, 2015.
    """

    _kT = Cst.k_ * _Te
    _rt = Cst.saha_ * _Te**(1.5) * (_gk/_gi) * np.exp(-_chi/_kT) / _ne

    return _rt

def get_LTE_ratio(_erg, _g, _stage, _Te, _Ne):
    r"""
    Compute LTE population ratio relative to 1st level.

    Parameters
    ----------

    _erg : numpy.1darray of np.double
        level energy relative to 1st level, [:math:`erg`]

    _weight : numpy.1darray of np.uint8
        statistical weight, [-]

    _stage : numpy.1darray of np.uint8
        ionization stage, [-]

    _Te : np.double
        electron temperature, [:math:`K`]

    _Ne : np.double
        electron density, [:math:`cm^{-3}`]

    Returns
    --------
    _nRatio : numpy.1darray of np.double
        population ratio. [-]

    Notes
    -----
    We set the population of the 1st level to `1.0`, and then compute the
    population ratio level by level. `Saha_distribution` is applied whenever
    there is a jump of ionization stage.
    Finnaly, we normalize the _nRatio array with respec to the total population.

    """
    _nLevel = _erg.size
    _nRatio = np.empty(_nLevel, np.double)
    _nRatio[0] = 1.
    #kT = Cst.k_*Te
    for i in range(1, _nLevel):
        _gj, _gi = _g[i], _g[i-1]
        if _stage[i] - _stage[i-1] == 0:
            #nRatio[i] = nRatio[i-1] * gj/gi * np.exp(-(erg[i]-erg[i-1])/kT)
            _nRatio[i] = _nRatio[i-1] * Boltzmann_distribution(_gi, _gj, _erg[i]-_erg[i-1], _Te)
        elif _stage[i] - _stage[i-1] == 1:
            #nRatio[i] = nRatio[i-1] * gj/gi * np.exp(-(erg[i]-erg[i-1])/kT) * Cst.saha_ * Te**(1.5) / Ne
            _nRatio[i] = _nRatio[i-1] * Saha_distribution(_gi, _gj, _erg[i]-_erg[i-1], _Ne, _Te)

        _nRatio[:] /= _nRatio[:].sum()

    return _nRatio


def EinsteinA_to_EinsteinBs_hz(Aji, f0, gi, gj):
    r"""

    given Einstein A coefficient Aij,
    calculate Einstein B coefficients Bij and Bji.

    Parameters
    ----------

    Aji : np.double or array-like
        Einstein A coefficient Aji, [:math:`s^{-1}`]
    f0 : np.double or array-like
        central frequency of corresponding line transition, [:math:`Hz`]
    gi : np.uint8 or array-like
        statistical weight of lower level, [-]
    gj : np.uint8 or array-like
        statistical weight of upper level, [-]

    Returns
    -------

    Bji : np.double or array-like
        Einstein B coefficient Bji, [:math:`s^{-1}/(erg/cm^{2}/Sr/Hz/s)`]
    Bij : np.double or array-like
        Einstein B coefficient Bji, [:math:`s^{-1}/(erg/cm^{2}/Sr/Hz/s)`]

    Notes
    -----

    Einstein Relations for bound-bound transitions [1]_.

    .. math:: B_{ji} = A_{ji} / \frac{2h\nu^{3}}{c^{2}}
    .. math:: B_{ij} = B_{ji} \frac{g_{j}}{g_{i}}

    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 117, 2015.
    """
    factor_ = 2*Cst.h_*f0**3/Cst.c_**2
    Bji = Aji / factor_
    Bij = Bji * gj / gi

    return Bji, Bij


def EinsteinA_to_EinsteinBs_cm(Aji, w0, gi, gj):
    r"""

    given Einstein A coefficient Aij,
    calculate Einstein B coefficients Bij and Bji.

    Parameters
    ----------

    Aji : np.double or array-like
        Einstein A coefficient Aji, [:math:`s^{-1}`]
    f0 : np.double or array-like
        central frequency of corresponding line transition, [:math:`Hz`]
    gi : np.uint8 or array-like
        statistical weight of lower level, [-]
    gj : np.uint8 or array-like
        statistical weight of upper level, [-]

    Returns
    -------

    Bji : np.double or array-like
        Einstein B coefficient Bji, [:math:`s^{-1}/(erg/cm^{2}/Sr/Hz/s)`]
    Bij : np.double or array-like
        Einstein B coefficient Bji, [:math:`s^{-1}/(erg/cm^{2}/Sr/Hz/s)`]

    Notes
    -----

    Einstein Relations for bound-bound transitions.
    Since we are in wavelength unit, the relation between
    [:math:`A_{ji}`] and [:math:`B_{ji}`] should be evaluated using
    Planck function also in wavelength unit

    .. math:: B_{ji} = A_{ji} / (2 h c^{2} \lambda^{5} )
    .. math:: B_{ij} = B_{ji} \frac{g_{j}}{g_{i}}

    """
    factor_ = 2*Cst.h_*Cst.c_**2/w0**5
    Bji = Aji / factor_
    Bij = Bji * gj / gi

    return Bji, Bij


def Planck_hz(F,T):
    r"""
    given frequency and temperature,
    calculate the frequency based planck function.

    with `nb.vectorize( [nb.float64(nb.float64,nb.float64)])`.

    Parameters
    ----------

    F : np.double or array-like
        frequency, [:math:`Hz`]
    T : np.double or array-like
        temperature, [:math:`K`]

    Returns
    -------

    intensity : np.double or array-like
        frequency based intensity, [:math:`erg/cm^2/Sr/Hz/s`]

    Notes
    -----

    The Planck function [1]_.

    .. math:: B_{\nu}(T) = \frac{2h\nu^3}{c^2} \frac{1}{e^{h\nu/kT}-1} \quad [erg/cm^2/Sr/Hz/s]

    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 102, 2015.
    """
    F_T = F / T
    if F_T > 1.04183e+13:                                        # ignore exponential part, exponential > 500
        intensity = 0
    elif F_T > 2.08366e+12:                                      # Wien approximation, exponential > 100
        intensity = 2.0*Cst.h_*F*F*F/Cst.c_/Cst.c_ * np.exp( -Cst.h_*F_T/(Cst.k_) )
    elif F_T < 2.08366e+08:                                      # Rayleighâ€“Jeans approximation, exponential < 0.01
        intensity = 2.0*F*F/Cst.c_/Cst.c_ * Cst.k_ * T
    else:                                                        # normal formula
        intensity = 2.0*Cst.h_*F*F*F/Cst.c_/Cst.c_ / ( np.exp( Cst.h_*F_T/(Cst.k_) ) - 1.0 )

    return intensity


def Planck_cm(W,T):
    r"""
    given wavelength and temperature,
    calculate the wavelength based planck function.

    with `nb.vectorize( [nb.float64(nb.float64,nb.float64)])`.

    Parameters
    ----------

    W : np.double or array-like
        wavelength, [:math:`cm`]
    T : np.double or array-like
        temperature, [:math:`K`]

    Returns
    -------

    intensity : np.double or array-like
        wavelength based intensity, [:math:`erg/cm^2/Sr/cm/s`]

    Notes
    -----

    The Planck function [1]_.

    .. math:: B_{\lambda}(T) = \frac{2hc^2}{\lambda^5} \frac{1}{e^{hc/\lambda kT}-1} \quad [erg/cm^2/Sr/cm/s]

    References
    ----------

    .. [1] Robert J. Rutten, "Introduction to Astrophysical Radiative Transfer", pp. 43, 2015.
    """
    WT = W*T
    if WT < 2.87755e-03:                                         # ignore exponential part, exponential > 500
        intensity = 0
    elif WT < 1.43878e-02:                                       # Wien approximation, exponential > 100
        intensity = 2.0*Cst.h_*Cst.c_*Cst.c_ /(W)**5 * np.exp( -Cst.h_*Cst.c_/(Cst.k_*WT) )
    elif WT > 1.43878e+02:                                       # Rayleighâ€“Jeans approximation, exponential < 0.01
        intensity = 2.0*Cst.c_*Cst.k_*T / (W)**4
    else:                                                        # normal formula
        intensity = 2.0*Cst.h_*Cst.c_*Cst.c_ /(W)**5 /  ( np.exp( Cst.h_*Cst.c_/(Cst.k_*WT) ) - 1.0 )

    return intensity


################################################################################
# whether to compile them using numba's LLVM
################################################################################

if Cst.isJIT == True :
    Planck_cm = nb.vectorize( [nb.float64(nb.float64,nb.float64)])( Planck_cm )
    Planck_hz = nb.vectorize( [nb.float64(nb.float64,nb.float64)])( Planck_hz )

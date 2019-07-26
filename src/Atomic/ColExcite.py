import numpy as np
from .. import Constants as Cst

from scipy.interpolate import splrep, splev

def interpolate_CE_fac(_table, _Te, _Te_table, _f1, _f2):
    r"""
    given temperature, interpolate collisional excitation coefficient
    ( in most case is Effective collisional strength)

    Parameters
    ----------

    _table : np.double, array, (nLine, nTemperature)
        a table of CE coefficient as a function of temperature for interpolation

    _Te : scalar
        termperature, [:math:`K`]

    _Te_table : np.double, array, (nTemperature,)
        corresponding termperature points in _table, [:math:`K`]

    _f1 : int
        a factor needed to compute CE rate coefficient

    _f2 : int
        a factor needed to compute CE rate coefficient

    Returns
    -------

    _CE_fac : np.double, np.array, (nLine,)
        the CE coefficient we need to compute CE rate coefficient

    Notes
    -----

    From [1]_.

    .. math: n_{e} C_{ij} = n_{e} \frac{8.63e-6 \times (\Omega_{ij} f1 / f2) }{g_{i} * T_{e}^{1/2}}  \exp{\frac{-dE_{ji}}{kT_{e}} } \quad [s^{-1}]

    References
    ----------

    .. [1] John T. Mariska, "The Solar Transition Region",
        Cambridge University Press, pp. 22, 1992
    """

    _nTran = _table.shape[0]
    _CE_fac = np.empty(_nTran, dtype=np.double)
    for k in range(_nTran):
        #--- numpy linear interpolation
        #_CE_fac[k] = np.interp(_Te, _Te_table, _table[k,:]) * _f1[k] / _f2[k]
        #--- scipy B-spline interpolation
        _Bsp_obj = splrep(x=_Te_table[:], y=_table[k,:])
        _CE_fac[k] = splev(_Te, _Bsp_obj, ext=3) * _f1[k] / _f2[k]

    return _CE_fac

def Cij_to_Cji(_Cij,  _ni_LTE, _nj_LTE):
    r"""
    calculate Cji from Cij

    Parameters
    -----------

    _Cij : np.double, np.array, (nLine); scalar
        collisional upward rate coefficient, [:math:`cm^{-3} s^{-1}`]

    _ni_LTE : np.double, np.array, (nLine); scalar
        population in lower level, [:math:`cm^{-3}`]

    _nj_LTE : np.double, np.array, (nLine); scalar
        population in upper level, [:math:`cm^{-3}`]

    Returns
    -------
    _Cji : np.double, np.array, (nLine); scalar
        collisional downward rate coefficient, [:math:`cm^{-3} s^{-1}`]

    Notes
    -----

    .. math: n_j^{LTE} C_{ji} = n_i^{LTE} C_{ij}

    """

    _Cji = _Cij * _ni_LTE / _nj_LTE

    return _Cji

def get_CE_rate_coe(_CE_fac, _Te, _gi, _dEij, _type):
    r"""
    compute the CE rate coefficient for a specific type.

    Parameters
    ----------

    _CE_fac : np.double, np.array, (nLine,); scalar
        the coefficient we interpolate from data

    _Te : scalar
        termperature, [:math:`K`]

    _gi : np.uint8, np.array, (nLine,); scalar
        statistical weight of lower level of the transition, [-]

    _dEij : np.double, np.array, (nLine); scalar
        excitation energy, [:math:`erg`]

    _type : str
        type of the data for interpolating _CE_fac, which decides the formula we will use.

    Returns
    -------

    _CEij : np.double, np.array, (nLine); scalar
        collisional excitation rate coefficient, [:math:`cm^{-3} s^{-1}`]

    Notes
    -----

    For "ECS" (Effective Collisional Strength) [1]_.

    .. math: n_{e} C_{ij} = n_{e} \frac{8.63e-6 \times (\Omega_{ij} f1 / f2) }{g_{i} * T_{e}^{1/2}}  \exp{\frac{-dE_{ji}}{kT_{e}} } \quad [s^{-1}]

    References
    ----------

    .. [1] John T. Mariska, "The Solar Transition Region",
        Cambridge University Press, pp. 22, 1992
    """

    _kT = Cst.k_ * _Te

    if _type == "ECS":

        _CEij = (8.63E-06 * _CE_fac) / (_gi * _Te**0.5) * np.exp( - _dEij / _kT )

    else:
        return None

    return _CEij

import numpy as np
from . import Constants as Cst

def interpolate_CE_fac(_table, _Te, _Te_table, _f1, _f2):

    _nTran = _table.shape[0]
    _CE_fac = np.empty(_nTran, dtype=np.double)
    for k in range(_nTran):
        _CE_fac[k] = np.interp(_Te, _Te_table, _table[k,:]) * _f1[k] / _f2[k]

    return _CE_fac

def Cij_to_Cji(_Cij,  _ni_LTE, _nj_LTE):

    _Cji = _Cij * _ni_LTE / _nj_LTE

    return _Cji

def get_CE_rate_coe(_CE_fac, _Te, _gi, _dEij, _type):

    _kT = Cst.k_ * _Te

    if _type == "ECS":

        _CEij = (8.63E-06 * _CE_fac) / (_gi * _Te**0.5) * np.exp( - _dEij / _kT )

    else:
        return None

    return _CEij

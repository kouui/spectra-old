import numpy as np
import Constants as Cst

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


if __name__ == "__main__":

    import AtomCls
    import LTELib

    file = "/Users/liu/kouui/workspace/statistical_equilibrium/atom/C_III_Be_like.txt"
    atom = AtomCls.Atom(file)
    Te = 2E+04
    ne = 1E+10

    #--- compute LTE population ratio for each CE transition
    n_LTE = LTELib.get_LTE_ratio(_erg=atom.Level.erg[:], _g=atom.Level.g[:],
                    _stage=atom.Level.stage[:], _Te=Te, _Ne=ne)
    nTran = atom.CE_table.shape[0]
    ni_LTE = np.empty(nTran, np.double)
    nj_LTE = np.empty(nTran, np.double)

    for k in range(nTran):
        ni_LTE[k] = n_LTE[atom.CE_coe.idxI[k]]
        nj_LTE[k] = n_LTE[atom.CE_coe.idxJ[k]]
    #---
    CE_fac = interpolate_CE_fac(_table=atom.CE_table[:,:], _Te=Te, _Te_table=atom.CE_Te_table[:],
                            _f1=atom.CE_coe.f1[:], _f2=atom.CE_coe.f2[:])
    CEij = get_CE_rate_coe(_CE_fac=CE_fac, _Te=Te, _gi=atom.CE_coe.gi[:],
                            _dEij=atom.CE_coe.dEij[:], _type=atom.CE_type[0])
    CEji = Cij_to_Cji(_Cij=CEij,  _ni_LTE=ni_LTE, _nj_LTE=nj_LTE)

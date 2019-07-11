import numpy as np


def setMatrixC(_Cmat, _Cji, _Cij, _idxI, _idxJ, _ne):
    r"""
    Compute the collisional rate matrix.

    Parameters
    ----------

    _Cmat : numpy.2darray of np.double
        collisional rate matrix, a 2D Array to store computed results.

    _Cji : numpy.1darray of np.double
        downward collisional transition rate, [:math:`s^{-1} \cdot cm^{-3}`]

    _Cij : numpy.1darray of np.double
        upward collisional transition rate, [:math:`s^{-1} \cdot cm^{-3}`]

    _idxI : numpy.1darray of np.uint16
        level index of lower level i, [-]

    _idxJ : numpy.1darray of np.uint16
        level index of upper level j, [-]

    _ne: np.double
        electron density, [:math:`cm^{-3}`]

    Notes
    -----
    Refer to [1]_ Equation(9.80).

    .. math:: \sum_{j \neq i} n_j (R_{ji}+C_{ji}) - n_i \sum_{j \neq i}(R_{ij}+C_{ij}) = 0


    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 282, 2015.

    """

    _n_row, _n_col = _Cmat.shape
    assert _n_row == _n_col, '_Cmat should be a squared matrix.'

    _nTran = _Cji.size
    for k in range(_nTran):
        i, j = _idxI[k], _idxJ[k]
        _Cmat[i,j] += _ne * _Cji[k]
        _Cmat[j,i] += _ne * _Cij[k]


def setMatrixR(_Rmat, _Rji_spon, _Rji_stim, _Rij, _idxI, _idxJ):
    r"""
    Compute the radiative rate matrix.

    Parameters
    ----------

    _Rmat : numpy.2darray of np.double
        radiative rate matrix, a 2D Array to store computed results.

    _Rji_spon : numpy.1darray of np.double
        spontaneous radiative transition rate, [:math:`s^{-1} \cdot cm^{-3}`]

    _Rji_stim : numpy.1darray of np.double
        spontaneous radiative transition rate, [:math:`s^{-1} \cdot cm^{-3}`]

    _Rij : numpy.1darray of np.double
        upward radiative transition rate, [:math:`s^{-1} \cdot cm^{-3}`]

    _idxI : numpy.1darray of np.uint16
        level index of lower level i, [-]

    _idxJ : numpy.1darray of np.uint16
        level index of upper level j, [-]

    Notes
    -----
    Refer to [1]_ Equation(9.80).

    .. math:: \sum_{j \neq i} n_j (R_{ji}+C_{ji}) - n_i \sum_{j \neq i}(R_{ij}+C_{ij}) = 0


    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 282, 2015.

    """
    _n_row, _n_col = _Rmat.shape
    assert _n_row == _n_col, '_Rmat should be a squared matrix.'

    _nTran = _Rji_spon.size
    for k in range(_nTran):
        i, j = _idxI[k], _idxJ[k]
        _Rmat[i,j] += _Rji_spon[k] + _Rji_stim[k]
        _Rmat[j,i] += _Rij[k]


def solveSE(_Rmat, _Cmat):

    _nLevel = _Rmat.shape[0]
    _A = _Cmat[:,:] + _Rmat[:,:]
    _b = np.zeros(_nLevel, dtype=np.double)

    #-------------------------------------------------------------
    # diagnal components
    #-------------------------------------------------------------
    for k in range(_nLevel):
        _A[k,k] = -_A[:,k].sum()

    #-------------------------------------------------------------
    # abundance definition equation
    #-------------------------------------------------------------
    _A[-1,:] = 1.
    _b[-1] = 1.

    _nArr = np.linalg.solve(_A,_b)

    return _nArr


if __name__ == "__main__":

    import AtomCls
    import LTELib
    import ColExcite
    import OpticallyThin

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

    #--- compute collision excitation/de-excitation rate coefficient
    CE_fac = ColExcite.interpolate_CE_fac(_table=atom.CE_table[:,:], _Te=Te, _Te_table=atom.CE_Te_table[:],
                            _f1=atom.CE_coe.f1[:], _f2=atom.CE_coe.f2[:])
    CEij = ColExcite.get_CE_rate_coe(_CE_fac=CE_fac, _Te=Te, _gi=atom.CE_coe.gi[:],
                            _dEij=atom.CE_coe.dEij[:], _type=atom.CE_type[0])
    CEji = ColExcite.Cij_to_Cji(_Cij=CEij,  _ni_LTE=ni_LTE, _nj_LTE=nj_LTE)

    #--- solve SE equations
    nLevel = atom.nLevel
    Cmat = np.zeros((nLevel, nLevel), np.double)
    Rmat = np.zeros((nLevel, nLevel), np.double)
    setMatrixC(_Cmat=Cmat[:,:], _Cji=CEji[:], _Cij=CEij[:],
                _idxI=atom.CE_coe.idxI[:], _idxJ=atom.CE_coe.idxJ[:], _ne=ne)

    Rji_stim = np.zeros(atom.Line.AJI[:].shape, np.double)
    Rij = np.zeros(atom.Line.AJI[:].shape, np.double)
    setMatrixR(_Rmat=Rmat[:,:], _Rji_spon=atom.Line.AJI[:],
        _Rji_stim=Rji_stim[:], _Rij=Rij, _idxI=atom.Line.idxI[:], _idxJ=atom.Line.idxJ[:])

    n_SE = solveSE(_Rmat=Rmat[:,:], _Cmat=Cmat[:,:])

    #-- compute optically thin relative flux
    nj_SE = np.empty(nTran, np.double)
    for k in range(nTran):
        nj_SE[k] = n_SE[atom.CE_coe.idxJ[k]]

    rel_flux = OpticallyThin.get_relative_flux(_AJI=atom.Line.AJI[:], _f0=atom.Line.f0[:], _nj=nj_SE[:])

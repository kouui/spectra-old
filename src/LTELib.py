import numpy as np
from . import Constants as Cst
from . import AtomCls

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

if __name__ == "__main__":

    import AtomCls

    file = "/Users/liu/kouui/workspace/spectra/atom/C_III_Be_like.txt"
    atom = AtomCls.Atom(file)
    nRatio = get_LTE_ratio(atom.Level.erg[:], atom.Level.g[:], atom.Level.stage[:], 2E+04, 1E+10)
    #print(nRatio)

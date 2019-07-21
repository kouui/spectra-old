import numpy as np
from .. import Constants as Cst

def get_relative_flux(_AJI, _f0, _nj):
    r"""
    calculate optically thin relative flux (some constants are removed)

    Parameters
    ----------

    _AJI : np.double, np.array, (nLine,) ; scalar
        Einstein A coefficient, [:math:`s^{-1}`]

    _f0 : np.double, np.array, (nLine,); scalar
        Transition line frequency, [:math:`hz`]

    _nj : np.double, np.array, (nLine); scalar
        Population of the upper level of line transition, [:math:`cm^{-3}`]

    Returns
    -------

    _rel_flux : np.double, np.array, (_atom.nLine,); scalar
        relative flux under the assumption of optically thin. [:math:`erg \; cm^{-3} \; s^{-1}`]

    Notes
    -----

    The absolute flux under the assumption of optically thin [1]_.

    .. math:: F_{ji} = \frac{1}{4 \pi R^{2}} \int_{\Delta V} \epsilon_{ji} dV \quad [erg \; cm^{-2} \; s^{-1}]

    where the emissivity :math:`\epsilon_{ji}` is

    .. math:: \epsilon_{ji} = h \nu n_{j} A_{ji} \quad [erg \; cm^{-3} \; s^{-1}]

    and

    .. math:: \epsilon_{ji} = \int_{\nu} \epsilon_{\nu} d \nu = \int_{\nu} h \nu n_{j} A_{ji} \psi d \nu

    So our relative flux is given by

    .. math: h \nu n_{j} A_{ji} \quad [erg \; cm^{-3} \; s^{-1}]

    References
    ----------

    .. [1] John T. Mariska, "The Solar Transition Region",
        Cambridge University Press, pp. 19, 1992

    """

    _rf = Cst.h_ * _f0[:] * _nj[:] * _AJI[:]

    return _rf

################################################################################
# this file defines functions for
#     calculations of opacities of various atomic process
#
################################################################################

import numpy as np

from .. import Constants as Cst


################################################################################
# Thomson scattering (wavelength independent)
################################################################################

def Thomson_scattering(n_e):
    r"""
    Thomson scattering of free electrons (non-relativistic).
    Returns absorption coefficient instead of absorption cross section.

    Parameters
    ----------

    n_e : np.double or array-like
        electron density, [:math:`cm^{-3}`]

    Returns
    -------

    kappa : np.double or array-like
        absorption coefficient, [:math:`cm^{-1}`]

    Notes
    -----

    .. math:: \kappa = \sigma_{T} n_{e}

    where :math:`\sigma_{T}` is the absorption cross section for Thomson scattering
    .. math:: \sigma_{T} = \frac{8 \pi e^{4}}{3 m_{e}^{2} c^{4}} = 6.6524 \times 10^{-25} \quad cm^{2}

    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 149, 2015.
    """

    kappa = 6.6524E-25 * n_e
    return kappa

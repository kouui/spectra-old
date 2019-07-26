################################################################################
# this file defines functions for
#     calculations related to spectral profile
#
################################################################################

import numpy as np
import numba as nb

from .. import Constants as Cst


################################################################################
# absorption profile
#    - Voigt
#    - Gaussian
################################################################################

def Voigt(a,x):
    r"""
    Calculate Doppler width normalized voigt function using polynomial fitting formula.

    `voigt(a,x)` function itself is not a vectorized function, only available to scalar operation.
    so we applied

    `nb.vectorize([nb.float64(nb.float64,nb.float64)],nopython=True)`.

    Parameters
    ----------

    a : np.double or array-like
        damping constant normalized by Doppler width, [-]
    x : np.double or array-like
        Doppler width normalized mesh, [-]

    Returns
    -------

    res : np.double or array-like
        voigt function, normalized to 1, [-]

    Notes
    -----

    The Voigt function is not normalized but has area :math:`\sqrt{\pi}` in `x` unit

    This is a combination of Humlicek(1982) [1]_ and Hui et al.(1978) [2]_ methods.

    When :math:`a > 0.1`, one can not ignore the wing component of Voigt function.
    That is, to guarantee its normalization, one has to take care of
    whether the mesh points are wide enough.

    References
    ----------

    .. [1] J.Humlicek, 'Optimized computation of the voigt and complex probability functions',
        Journal of Quantitative Spectroscopy and Radiative Transfer (JQSRT),
        Volume 27, Issue 4, April 1982, Pages 437-444.

    .. [2] A.K.Hui, B.H.Armstrong, A.A.Wray, 'Rapid computation of the Voigt and complex error functions',
        Journal of Quantitative Spectroscopy and Radiative Transfer (JQSRT),
        Volume 19, Issue 5, May 1978, Pages 509-516.
    """

    if a < 0.01:
        Z = a - 1j*x
        U = Z * Z
        A0 = 36183.31; B0 = 32066.6
        A1 = 3321.9905;B1 = 24322.84
        A2 = 1540.787; B2 = 9022.228
        A3 = 219.0313; B3 = 2186.181
        A4 = 35.76683; B4 = 364.2191
        A5 = 1.320522; B5 = 61.57037
        A6 = .56419  ; B6 = 1.841439
        F = (np.exp(U) - Z * ( A0 - U * ( A1 - U * ( A2 - U * ( A3 - U * ( A4 - U * ( A5 - U * A6 )))))) /
         ( B0 - U * ( B1 - U * ( B2 - U * ( B3 - U * ( B4 - U * ( B5 - U * ( B6 - U ))))))))
    else:
        A0 = 122.607931777; B0 = 122.607931774
        A1 = 214.382388695; B1 = 352.730625111
        A2 = 181.928533092; B2 = 457.334478784
        A3 = 93.1555804581; B3 = 348.703917719
        A4 = 30.1801421962; B4 = 170.354001821
        A5 = 5.91262620977; B5 = 53.9929069129
        A6 = .564189583563; B6 = 10.4798571143

        Z = a - 1j*abs(x)# * C.k[0]
        F = ( ( A0 + Z * ( A1 + Z * ( A2 + Z * ( A3 + Z * ( A4 + Z * ( A5 + Z * A6 ) ) ) ) ) ) /
             ( B0 + Z * ( B1 + Z * ( B2 + Z * ( B3 + Z * ( B4 + Z * ( B5 + Z * ( B6 + Z ) ) ) ) ) ) ) )

    res = F.real / Cst.sqrtPi_
    return res

def Gaussian(x):
    r"""
    Calculate Doppler width normalized gaussian profile.

    Parameters
    ----------

    x : np.double or array-like
        Doppler width normalized mesh, [-]

    Returns
    -------

    res : np.double or array-like
        gaussian profile, normalized to 1, [-]

    Notes
    -----

    This formula refers to [1]_.

    References
    ----------

    .. [1] Ivan Hubeny, Dimitri Mihalas, "Theory of Stellar Atmosphere:
        An Introduction to Astrophysical Non-equilibrium
        Quantitative Spectroscopic Analysis",
        Princeton University Press, pp. 203, Eq(8.24), 2015.
    """
    res = np.exp(-x*x) / Cst.sqrtPi_

    return res

################################################################################
# whether to compile them using numba's LLVM
################################################################################

if Cst.isJIT == True:
    Voigt = nb.vectorize( [nb.float64(nb.float64,nb.float64)],nopython=True)( Voigt )

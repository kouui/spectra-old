################################################################################
# this file defines functions for
#     calculations related to naive/basic physics process
################################################################################

from .. import Constants as Cst


################################################################################
# wavelength [cm] <--> frequency [hz]
################################################################################

def wave_to_freq(wave):
    r"""
    convert Wavelength to Frequency.

    Parameters
    ----------

    wave : np.double or array-like
        Wavelength, [:math:`cm`]

    Returns
    -------

    freq : np.double or array-like
        Frequency, [:math:`Hz`]

    Notes
    -----

    .. math:: \nu = \frac{c}{\lambda}
    """
    freq = Cst.c_ / wave
    return freq

def freq_to_wave(freq):
    r"""
    convert Frequency to Wavelength.

    Parameters
    ----------

    freq : np.double or array-like
        Frequency, [:math:`Hz`]

    Returns
    -------

    wave : np.double or array-like
        Wavelength, [:math:`cm`]

    Notes
    -----

    .. math:: \lambda = \frac{c}{\nu}
    """
    wave = Cst.c_ / freq
    return wave

################################################################################
# Doppler shift and Doppler width
################################################################################

def Dvlocity_to_Dshift(p0, v):
    r"""
    given Doppler velocity and the line central wavelength/frequency,
    compute Doppler shift in wavelength/frequency.

    Parameters
    ----------

    p0 : np.double or array-like
        Frequency in any frequency unit or Wavelength in any length unit
    v : np.double or array-like
        line of sight velocity, [:math:`cm/s`]

    Returns
    -------

    dp : np.double or array-like
        Frequency or Wavelength, same unit with input `p0`

    Notes
    -----

    in wavelength,

    .. math:: d\lambda = \lambda_0 \frac{v}{c}

    in frequency,

    .. math:: d\nu = \nu_0 \frac{v}{c}
    """
    dp = p0 * v / Cst.c_
    return dp

def get_Doppler_width(p0, Te, Vt, am):
    r"""
    Given central wavelength/frequency, relative atomic mass of a line,
    and the temperature, turbulent velocity, compute the corresponding
    Doppler Width.

    Parameters
    ----------

    p0 : np.double or array-like
        Frequency in any frequency unit or Wavelength in any length unit
    Te : np.double or array-like
        Temperature, [:math:`K`]
    Vt : np.double or array-like
        Turbulent velocity, [:math:`cm/s`]
    am : np.double or array-like
        atomic mass relative to hydrogen atom. [-]

    Returns
    -------

    dp : np.double or array-like
        Frequency/Wavelength Doppler width, same unit with input argument `p0`

    Notes
    -----

    According to [1]_.

    .. math:: \eta_{0} = (2kT/m + Vt^2)^{1/2}

    in wavelength,

    .. math:: d\lambda = \lambda_0 \frac{\eta_{0}}{c}

    in frequency,

    .. math:: d\nu = \nu_0 \frac{\eta_{0}}{c}

    References
    -----------

    .. [1] Robert J. Rutten, "Radiative Transfer in Stellar Atmosphere", pp. 82, 2003.
    """
    eta0_ = (2.*Cst.k_*Te/(Cst.mH_*am) + Vt*Vt)**(0.5)
    dp = p0 * eta0_ / Cst.c_
    return dp

import numpy as np
from .. import Constants as Cst

def get_relative_flux(_AJI, _f0, _nj):
    r"""

    """

    _rf = Cst.h_ * _f0[:] * _nj[:] * _AJI[:]

    return _rf

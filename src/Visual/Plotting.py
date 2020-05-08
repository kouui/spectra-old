"""
Library to support plotting using matplotlib
"""

import numpy as np

def set_imshow_ticks(axe, arr, axis, points=None, fmt='%1.3f', rot=0, fontsize=None):
    r"""
    customize ticks and ticklabels for a specific axe.

    Parameters
    ----------
    axe : matplotlib.pyplot.Axes
        the axe whose ticks and ticklabels we are going to modify.

    arr : array-like, numpy.ndarray,
        the array to create ticklabels.

    axis : 'x' or 'y'
        modify ticks and ticklabels of 'x' or 'y' axis.

    points : int or list of int or None, optional
        int : nbins of ticks
        list of int : list of ticks
        None : use the default ticks
        default : None

    fmt : string Formattor
        format to for string formatting ticklabels.

    rot : angle, [:math:`^\circ`]
        rotation angle of ticklabels.


    Returns
    -------
    None

    Note
    -----
    None

    """

    assert axis in ('x', 'y')

    #-- if integer format requied, convert array to np.int64
    if fmt[-1] == 'd':
        arr_ = arr.astype(np.int64)
    else:
        arr_ = arr.copy()

    #-- if points is integer, create equally spaced points
    if isinstance(points, int):
        temp = np.linspace(0, arr_.shape[0]-1, points)
        #points = ((temp[1:] + temp[:-1]) * 0.5).astype(np.int64)
        points = temp.astype(np.int64)
    #-- if points is None, use the default ticks
    elif points is None:
        if axis in ('x',):
            points = axe.get_xticks()[:-1].astype(np.int64)
        elif axis in ('y',):
            points = axe.get_yticks()[:-1].astype(np.int64)
        points = points[points>=0]

    #-- format ticklabels
    ticklabels = [("{:" + f"{fmt[1:]}" + "}").format(arr_[i]) for i in points]

    #-- set ticks and ticklabels
    if axis in ('x',):
        axe.set_xticks(points)
        axe.set_xticklabels(ticklabels, rotation=rot, fontsize=fontsize)
    elif axis in ('y',):
        axe.set_yticks(points)
        axe.set_yticklabels(ticklabels, rotation=rot, fontsize=fontsize)

    return None

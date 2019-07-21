import matplotlib.pyplot as plt
from collections import Counter
from .. import Constants as Cst
import math

def prepare_dict(_atom, _conf_duplicate):
    r"""
    separate singlet and multiplet

        - singlet   : { (conf_clean,term) : { J : (E[eV], L, conf_origin) } }
        - multiplet : { (conf_clean,term) : { J : (E[eV], L, conf_origin) } }

    Parameters
    ----------
    _atom : AtomCls.Atom
        object of the atomic model

    _conf_duplicate : str
        common configuration string of the inner shell

    Returns
    -------

    _singlet : { (conf_clean,term) : { J : (E[eV], L, conf_origin) } }

    _multiplet : { (conf_clean,term) : { J : (E[eV], L, conf_origin) } }

    _Lset : { "singlet" : {int,}, "multiplet" : {int,} }
        set of Quantum number L in integer for singlet and multiplet, respectively.
    """
    _L_duplicate = len(_conf_duplicate)

    #-------------------------------------------------------------------------
    # create and count list of (conf, term)
    #-------------------------------------------------------------------------
    _conf_term = list( zip( _atom.Level_info["configuration"], _atom.Level_info["term"] ) )
    _count = Counter( _conf_term )
    #-------------------------------------------------------------------------

    #-------------------------------------------------------------------------
    # separate singlet and multiplet
    #-------------------------------------------------------------------------
    _L_duplicate = len(_conf_duplicate)

    _singlet = {}
    _multiplet = {}
    _Lset = {"singlet" : set(), "multiplet" : set()}
    for k in range(_atom.nLevel):

        #---------------------------------------------------------------------
        # - remove inner shell electron configuration
        # - remove all '.' in conf
        #---------------------------------------------------------------------
        _conf = _atom.Level_info["configuration"][k]
        _conf_clean = _conf[_L_duplicate:].replace('.','')
        #---------------------------------------------------------------------
        _term = _atom.Level_info["term"][k]
        _J = _atom.Level_info["J"][k]
        _L = Cst.L_s2i[ _term[-1] ]

        if _count[ (_conf, _term) ] == 1:
            _d = _singlet
            _Lset["singlet"].add(_L)
        else:
            _d = _multiplet
            _Lset["multiplet"].add(_L)

        _key = (_conf_clean ,_term)
        if _key not in _d.keys():
            _d[_key] = {}
        _d[_key][_J] = ( _atom.Level.erg[k] / Cst.eV2erg_, _L, _conf )
    #-------------------------------------------------------------------------
    return _singlet, _multiplet, _Lset

def line_with_text(_ax, _lp, _rp, _text, _tsize, _r):
    r"""
    Parameters
    ----------

    _ax : matlotlib.pyplot.Axe
        the axe to plot a line

    _lp : tuple of float/int, (x,y)
        left point of the line

    _rp : tuple of float/int, (x,y)
        right point of the line

    _text : str
        string of the text

    _tsize : int
        texture size

    _r : float, in the range of [0 : 1]
        relative position along the line starting from left point

    Returns
    -------

    _line_obj : matplotlib.lines.Line2D
        object of the line we plot

    _text_obj : matplotlib.text.Text
        object of the text we plot
    """
    # {'color': 'black', 'fontsize': 6, 'ha': 'center', 'va': 'center',
    #    'bbox': dict(boxstyle="round", fc="white", ec="white", pad=0.2)}
    _line_obj, = _ax.plot([_lp[0], _rp[0]], [_lp[1], _rp[1]], "-", linewidth=0.7)
    _angle = math.atan( (_rp[1]-_lp[1]) / (_rp[0]-_lp[0]) )
    _angle = math.degrees( _angle )
    _tx = _lp[0] + (_rp[0]-_lp[0]) * _r
    _ty = _lp[1] + (_rp[1]-_lp[1]) * _r
    _text_obj = _ax.text( _tx, _ty, _text, {'color': 'black', 'fontsize': _tsize, 'ha': 'center', 'va': 'center',
                        'bbox': dict(boxstyle="round", fc="white", ec="white", pad=0.2)} )

    return _line_obj, _text_obj

class Grotrian:

    def __init__(self, _atom, _conf_duplicate="1s2."):

        self.atom = _atom

        #---------------------------------------------------------------------
        # prepare structures for plotting
        #---------------------------------------------------------------------
        singlet, multiplet, Lset = prepare_dict(_atom=_atom, _conf_duplicate=_conf_duplicate)
        self.singlet = singlet
        self.multiplet = multiplet
        self.Lset = Lset
        #---------------------------------------------------------------------

        #---------------------------------------------------------------------
        # a dictionary to store level position of levels in the plot.
        #---------------------------------------------------------------------
        self.pos_level = {}
        #---------------------------------------------------------------------

    def make_fig(self, _figsize=(6,8), _dpi=120, _f=200):
        # _f : enlarge factor insdie fine-structure

        Lset = self.Lset
        singlet = self.singlet
        multiplet = self.multiplet
        pos_level = self.pos_level

        #---------------------------------------------------------------------
        # config
        #---------------------------------------------------------------------
        _hw = 0.5                       # half width of a level in the plot
        _b = len(Lset["singlet"])       # bias in the x axis of multiplet
        _fontsize = 14                  # fontsize of labels
        _textsize = 10                  # fontsize of text of term
        _Jsize = 7                      # fontsize of text of J
        _st = 0.1                       # horizontal spacing between term and text
        #---------------------------------------------------------------------

        fig = plt.figure(figsize=_figsize, dpi=_dpi)
        self.fig = fig

        #---------------------------------------------------------------------
        # singlet
        #---------------------------------------------------------------------
        for k0, v0 in singlet.items():
            for k1, v1 in v0.items():
                #plot level
                xs_level = v1[1]-_hw, v1[1]+_hw
                ys_level = v1[0], v1[0]
                plt.plot(xs_level, ys_level, "-k", linewidth=1)
                # store level posiiton, (conf_origin, term, J)
                pos_level[ (v1[2], k0[1], k1) ] = {}
                pos_level[ (v1[2], k0[1], k1) ]["xs"] = xs_level
                pos_level[ (v1[2], k0[1], k1) ]["ys"] = ys_level
                # plot text
                x_text = xs_level[1] + _st
                y_text = ys_level[0]
                plt.text(x_text, y_text, "{} {}".format(k0[0],k0[1]), fontsize=_textsize, color="k")

        # a vertical line separates singlet panel and multiplet panel
        plt.axvline(x=_b+1, linestyle="--", linewidth=0.5, color="k")
        #---------------------------------------------------------------------

        #---------------------------------------------------------------------
        # plotting multiplet
        #---------------------------------------------------------------------
        _hw = 0.6                       # multiplet need wider
        _sf = 0.1                    # space for connecting fine structure
        for k0, v0 in multiplet.items():
            # compute mean term energy
            y_count = 0.
            y_mean = 0
            for k1, v1 in v0.items():
                y_mean += v1[0]
                y_count += 1
            y_mean /= y_count

            for k1, v1 in v0.items():
                # plot level
                y_pos = y_mean + _f * (v1[0]-y_mean)            # enlarge space
                ys_level = y_pos, y_pos
                xs_level = v1[1]+1-_hw+_b, v1[1]+1+_hw+_b-_sf
                plt.plot(xs_level, ys_level, "-k", linewidth=1)
                # store level posiiton, (conf_origin, term, J)
                pos_level[ (v1[2], k0[1], k1) ] = {}
                pos_level[ (v1[2], k0[1], k1) ]["xs"] = xs_level
                pos_level[ (v1[2], k0[1], k1) ]["ys"] = ys_level
                # connect fine structure
                plt.plot([xs_level[1],xs_level[1]+_sf], [y_pos,y_mean], "--k", linewidth=0.5)
                # plot text of term
                x_pos_text = xs_level[1]+_sf + _st
                plt.text(x_pos_text, y_mean, "{} {}".format(k0[0],k0[1]), fontsize=_textsize, color="k")
                # plot text of J
                plt.text(xs_level[0]-2*_st, y_pos, k1, fontsize=_Jsize, color="k")
        #---------------------------------------------------------------------

        #---------------------------------------------------------------------
        # config of the plot
        #---------------------------------------------------------------------

        #--- tick of singlet
        xtick1 = list(Lset["singlet"])
        xticklabel1 = [Cst.L_i2s[xt] for xt in xtick1]
        #--- tick of multiplet
        xtick2 = [ x+_b+1 for x in list(Lset["multiplet"]) ]
        xticklabel2 = [Cst.L_i2s[xt-_b-1] for xt in xtick2]
        #--- customize xticks and xticklabels
        plt.xticks( xtick1 + xtick2, xticklabel1 + xticklabel2, fontsize=_fontsize )

        #--- x,y label; title
        plt.xlabel("L", fontsize=_fontsize)
        plt.ylabel("E [eV]", fontsize=_fontsize, rotation=90)
        plt.title(self.atom.title, fontsize=_fontsize, y=0.93)

        #--- change x limit
        plt.xlim(-1, xtick2[-1]+2)
        #---------------------------------------------------------------------

    def save_fig(self, _filename, _dpi=120):

        self.fig.savefig(_filename, dpi=_dpi)

    def show_fig(self):

        plt.show()

    def connect_line(self, _cfj1, _cfj2, _r1, _r2, _c):
        r"""

        Parameters
        ----------

        _cfj1 : tuple of str, (conf,term,J)
            to specify level 1

        _cfj2 : tuple of str, (conf,term,J)
            to specify level 2

        _r1 : float, in the range of [0 : 1]

            relative position along the line starting from left point

        """
        pass

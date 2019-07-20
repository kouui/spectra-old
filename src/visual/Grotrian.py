import matplotlib.pyplot as plt
from collections import Counter
from .. import Constants as Cst

def prepare_dict(_atom, _conf_duplicate):
    r"""
    separate singlet and multiplet

        - singlet   : { (conf,term) : { J : (E[eV], L) } }
        - multiplet : { (conf,term) : { J : (E[eV], L) } }
    """
    _L_duplicate = len(_conf_duplicate)

    #---
    # create and count list of (conf, term)
    #---
    _conf_term = list( zip( _atom.Level_info["configuration"], _atom.Level_info["term"] ) )
    _count = Counter( _conf_term )

    #---
    # separate singlet and multiplet
    #---
    _L_duplicate = len(_conf_duplicate)

    _singlet = {}
    _multiplet = {}
    _Lset = {"singlet" : set(), "multiplet" : set()}
    for k in range(_atom.nLevel):

        #---
        # - remove inner shell electron configuration
        # - remove all '.' in conf
        #---
        _conf = _atom.Level_info["configuration"][k]
        _conf_clean = _conf[_L_duplicate:].replace('.','')
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
        _d[_key][_J] = ( _atom.Level.erg[k] / Cst.eV2erg_, _L )


    return _singlet, _multiplet, _Lset

def line_with_text(_text):
    # {'color': 'black', 'fontsize': 6, 'ha': 'center', 'va': 'center',
    #    'bbox': dict(boxstyle="round", fc="white", ec="white", pad=0.2)}
    pass

class Grotrian:

    def __init__(self, _atom, _conf_duplicate="1s2."):

        self.atom = _atom

        #---
        # prepare structures for plotting
        #---
        singlet, multiplet, Lset = prepare_dict(_atom=_atom, _conf_duplicate=_conf_duplicate)
        self.singlet = singlet
        self.multiplet = multiplet
        self.Lset = Lset

    def make_fig(self, _figsize=(6,8), _dpi=120, _f=200):
        # _f : enlarge factor insdie fine-structure

        Lset = self.Lset
        singlet = self.singlet
        multiplet = self.multiplet

        #---
        # config
        #---
        _hw = 0.5                       # half width of a level in the plot
        _b = len(Lset["singlet"])       # bias in the x axis of multiplet
        _fontsize = 14                  # fontsize of labels
        _textsize = 10                  # fontsize of text of term
        _Jsize = 7                      # fontsize of text of J


        fig = plt.figure(figsize=_figsize, dpi=_dpi)

        #---
        # singlet
        #---
        for k0, v0 in singlet.items():
            for k1, v1 in v0.items():
                plt.plot([v1[1]-_hw,v1[1]+_hw], [v1[0],v1[0]], "-k", linewidth=1)
                x_pos_text = v1[1]+_hw + 0.1
                plt.text(x_pos_text, v1[0], "{} {}".format(k0[0],k0[1]), fontsize=_textsize, color="k")


        plt.axvline(x=_b+1, linestyle="--", linewidth=0.5, color="k")

        #---
        # multiplet
        #---
        _hw = 0.6                       # multiplet need wider
        for k0, v0 in multiplet.items():
            #---
            # compute mean term energy
            #---
            y_count = 0.
            y_mean = 0
            for k1, v1 in v0.items():
                y_mean += v1[0]
                y_count += 1
            y_mean /= y_count

            for k1, v1 in v0.items():
                y_pos = y_mean + _f * (v1[0]-y_mean)            # enlarge space
                plt.plot([v1[1]+1-_hw+_b,v1[1]+1+_hw+_b-0.1], [y_pos,y_pos], "-k", linewidth=1)
                plt.plot([v1[1]+1+_hw+_b-0.1,v1[1]+1+_hw+_b], [y_pos,y_mean], "--k", linewidth=0.5)
                x_pos_text = v1[1]+1+_hw+_b + 0.1
                plt.text(x_pos_text, y_mean, "{} {}".format(k0[0],k0[1]), fontsize=_textsize, color="k")
                plt.text(v1[1]+1-_hw+_b-0.2,y_pos, k1, fontsize=_Jsize, color="k")

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

        self.fig = fig

    def save_fig(self, _filename, _dpi=120):

        self.fig.savefig(_filename, dpi=_dpi)

    def show_fig(self):

        plt.show()

    def connect_line(self):

        pass

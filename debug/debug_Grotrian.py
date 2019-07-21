

if __name__ == "__main__":

    import sys
    sys.path.append("..")

    from src.Structure import AtomCls
    from src.Visual import Grotrian
    #-------------------------------------------------------------------------
    # C III Be like
    #-------------------------------------------------------------------------
    ##filename = "/Users/liu/kouui/workspace/spectra/atom/C_III_Be_like.txt"
    ##atom = AtomCls.Atom(filename)
    ##gro = Grotrian.Grotrian(_atom=atom, _conf_duplicate="1s2.")
    ##gro.make_fig(_figsize=(10,6))
    ##
    ##line_plot = (
    ##    (0, 4, "977", 0.3, 0.5),
    ##    (0, 2, "1909", 0.7, 0.1),
    ##    (4, 9, "1247", 0.3, 0.5),
    ##    (4, 8, "2297", 0.8, 0.5),
    ##    (2, 5, "1175.99", 0.1, 0.1),
    ##    (1, 6, "1175.26", 0.25, 0.25),
    ##    (2, 6, "1175.59", 0.40, 0.40),
    ##    (3, 6, "1176.38", 0.55, 0.55),
    ##    (2, 7, "1174.93", 0.70, 0.70),
    ##    (3, 7, "1175.71", 0.85, 0.85),
    ##)
    #-------------------------------------------------------------------------

    #-------------------------------------------------------------------------
    # O V Be like
    #-------------------------------------------------------------------------
    ##filename = "/Users/liu/kouui/workspace/spectra/atom/O_V_Be_like.txt"
    ##atom = AtomCls.Atom(filename)
    ##gro = Grotrian.Grotrian(_atom=atom, _conf_duplicate="1s2.")
    ##gro.make_fig(_figsize=(10,6), _f=50)
    ##
    ##line_plot = (
    ##    (0, 4, "629.7", 0.3, 0.5),
    ##    (0, 2, "1218.3", 0.7, 0.1),
    ##    (4, 9, "774.5", 0.3, 0.5),
    ##    (4, 8, "1371.3", 0.8, 0.5),
    ##    (2, 5, "761.1", 0.1, 0.1),
    ##    (1, 6, "759.4", 0.25, 0.25),
    ##    (2, 6, "760.2", 0.40, 0.40),
    ##    (3, 6, "762.0", 0.55, 0.55),
    ##    (2, 7, "758.7", 0.70, 0.70),
    ##    (3, 7, "760.4", 0.85, 0.85),
    ##)
    #-------------------------------------------------------------------------

    #-------------------------------------------------------------------------
    # Si III Mg like, problem with 3s.4s 3S and 3s.3d 3D term
    #-------------------------------------------------------------------------
    ##filename = "/Users/liu/kouui/workspace/spectra/atom/Si_III_Mg_like.txt"
    ##atom = AtomCls.Atom(filename)
    ##gro = Grotrian.Grotrian(_atom=atom, _conf_duplicate="")
    ##gro.make_fig(_figsize=(10,6), _f=50)
    ##
    ##line_plot = (
    ##    (0, 4, "1206.5", 0.5, 0.5),
    ##)
    #-------------------------------------------------------------------------
    for i, j, wl, _r1, _r2 in line_plot:
        _cfj1 = atom.Level_info_table[i]
        _cfj2 = atom.Level_info_table[j]
        gro.connect_line(_cfj1=_cfj1, _cfj2=_cfj2, _r1=_r1, _r2=_r2, _c="black", _text=wl, _tsize=7, _r=0.4)
    gro.show_fig()
    #gro.save_fig("../grotrian_diagram/" + filename.split('/')[-1].replace(".txt",".png"))

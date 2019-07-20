
if __name__ == "__main__":

    import sys
    sys.path.append("..")

    from src import AtomCls, LTELib

    file = "/Users/liu/kouui/workspace/spectra/atom/C_III_Be_like.txt"
    atom = AtomCls.Atom(file)
    nRatio = LTELib.get_LTE_ratio(atom.Level.erg[:], atom.Level.g[:], atom.Level.stage[:], 2E+04, 1E+10)

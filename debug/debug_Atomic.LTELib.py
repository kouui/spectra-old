
if __name__ == "__main__":

    import sys
    sys.path.append("..")

    from src.Structure import AtomCls
    from src.Atomic import LTELib

    file = "/Users/liu/kouui/workspace/spectra/atom/C_III_Be_like.txt"
    atom = AtomCls.Atom(file)
    nRatio = LTELib.get_LTE_ratio(atom.Level.erg[:], atom.Level.g[:], atom.Level.stage[:], 2E+04, 1E+10)

    wave = 6.563E-5 # [cm]
    Te = 1.E+4 # [K]
    print(LTELib.EinsteinA_to_EinsteinBs_cm(0, wave, 2,8))
    print(LTELib.Planck_cm(wave, Te))

    

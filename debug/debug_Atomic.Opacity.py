
if __name__ == "__main__":

    import sys
    sys.path.append("..")

    from src.Atomic import Opacity

    Ne = 1.E+10 # cm^-3
    print(Opacity.Thomson_scattering(Ne))

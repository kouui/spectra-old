if __name__ == "__main__":

    import sys
    sys.path.append("..")

    from src.Math import Integrate

    import numpy as np

    x = np.arange(5)
    y = np.arange(5)

    print(Integrate.Trapze(y,x))

if __name__ == "__main__":

    import sys
    sys.path.append("..")

    from src.RadiativeTransfer import Profile
    #from src import Constants as Cst
    #from src.Math import Integrate

    import numpy as np
    import matplotlib.pyplot as plt


    x = np.linspace(-3, 3, 61, endpoint=True, dtype=np.double)
    a = 0.2
    y = {
        "Voigt"    : (Profile.Voigt(a, x), "red"),
        "Gaussian" : (Profile.Gaussian(x), "blue")
    }

    fig = plt.figure(figsize=(8,5), dpi=80)
    for key, value in y.items():
        plt.plot(x, value[0], color=value[1], label=key)
    plt.legend(loc="best")
    plt.show()

    #print(Integrate.Trapze(y["Voigt"][0],x), Cst.sqrtPi_)

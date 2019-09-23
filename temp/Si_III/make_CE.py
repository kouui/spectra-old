
import pandas as pd
from read_tables import readAll

import numpy as np

if __name__ == "__main__":

    df = readAll()
    df = df["CE"]

    TeLog_list = ["4.10","4.30","4.50","4.70","4.90","5.10","5.30","5.50","5.70","5.90"]
    Te_arr = np.array([float(i) for i in TeLog_list])
    Te_arr = 10**Te_arr

    _s_out  = "#   "
    _s_out += "{:<20s}".format("configuration_i") + " "*4
    _s_out += "{:<8s}".format("term_i") + " "*4
    _s_out += "{:<4s}".format("J_i") + " "*4
    _s_out += "{:<20s}".format("configuration_j") + " "*4
    _s_out += "{:<8s}".format("term_j") + " "*4
    _s_out += "{:<4s}".format("J_j") + " "*4
    for k in range(10):
        _s_out += "{:<1.2E}".format(Te_arr[k]) + " "*4
    _s_out += "{:<2s}".format("f1") + " "*4
    _s_out += "{:<2s}".format("f2") + " "*4
    _s_out += "\n"

    for i in range(len(df)):
        _s_out += "    "
        _s_out += "{:<20s}".format(df["i_Configuration"].iloc[i]) + " "*4
        _s_out += "{:<8s}".format(df["i_Term"].iloc[i]) + " "*4
        _s_out += "{:<4s}".format(df["i_J"].iloc[i]) + " "*4
        _s_out += "{:<20s}".format(df["j_Configuration"].iloc[i]) + " "*4
        _s_out += "{:<8s}".format(df["j_Term"].iloc[i]) + " "*4
        _s_out += "{:<4s}".format(df["j_J"].iloc[i]) + " "*4
        for k in range(10):
            _s_out += "{:<1.2E}".format(df[TeLog_list[k]].iloc[i]) + " "*4
        _s_out += "{:<2s}".format(str(1)) + " "*4
        _s_out += "{:<2s}".format(str(1)) + " "*4

        _s_out += '\n'

    _s_out += '\n'


    with open("./Kanti_2017.Electron", "w") as _f:
        _f.write(_s_out)





    #print(df.head())

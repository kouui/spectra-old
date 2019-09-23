
import pandas as pd
from read_tables import readAll


if __name__ == "__main__":

    df = readAll()
    df = df["Aji"]


    # format output string
    _s_out  = "#   "
    _s_out += "{:<40s}".format("configuration_i") + " "*4
    _s_out += "{:<8s}".format("term_i") + " "*4
    _s_out += "{:<4s}".format("J_i") + " "*4
    _s_out += "{:<40s}".format("configuration_j") + " "*4
    _s_out += "{:<8s}".format("term_j") + " "*4
    _s_out += "{:<4s}".format("J_j") + " "*4
    _s_out += "{:<9s}".format("Aji[s^-1]") + " "*4
    _s_out += "{:<14s}".format("Wavelength[AA]") + " "*4
    _s_out += "\n"

    for i in range(len(df)):
        _s_out += "    "
        _s_out += "{:<40s}".format(df["i_Configuration"].iloc[i]) + " "*4
        _s_out += "{:<8s}".format(df["i_Term"].iloc[i]) + " "*4
        _s_out += "{:<4s}".format(df["i_J"].iloc[i]) + " "*4
        _s_out += "{:<40s}".format(df["j_Configuration"].iloc[i]) + " "*4
        _s_out += "{:<8s}".format(df["j_Term"].iloc[i]) + " "*4
        _s_out += "{:<4s}".format(df["j_J"].iloc[i]) + " "*4
        _s_out += "{:<1.3e}".format(df["Aji_E1"].iloc[i]+df["Aji_E2"].iloc[i]+df["Aji_M1"].iloc[i]+df["Aji_M2"].iloc[i]) + " "*4
        _s_out += "{:<1.7e}".format(df["lamda_ij(AA)"].iloc[i]) + " "*4

        _s_out += '\n'
    _s_out += '\n'

    with open("./Kanti_2017.Aji", "w") as _f:
        _f.write(_s_out)

    #print(df.head())

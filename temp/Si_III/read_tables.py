
import pandas as pd

def readAll():

    df = {
        "Level" : None,
        "Aji" : None,
        "CE" : None,
    }

    file = "./Kanti_2017_Level.txt"
    dtype = {
        "Index" : int,
        "Configuration" : str,
        "Level" : str,
        "NIST" : float,
        "GRASP1a" : float,
        "GRASP1b" : float,
        "GRASP2" : float,
        "FAC1" : float,
        "FAC2" : float,
        "AS" : float,
    }
    df["Level"] = pd.read_csv(file, sep=' ', header=0, dtype=dtype, na_values=["-------",], index_col="Index")

    add_prefix = lambda x : "1s2.2s2.2p6." + x
    df["Level"]["Configuration"] = df["Level"]["Configuration"].apply(add_prefix)

    get_term = lambda x : x[:-1]
    get_J = lambda x : x[-1]
    df["Level"]["Term"] = df["Level"]["Level"].apply(get_term)
    df["Level"]["J"] = df["Level"]["Level"].apply(get_J)

    file = "./Kanti_2017_Aji.txt"
    dtype = {
        "i" : int,
        "j" : int,
        "lamda_ij(AA)" : float,
        "Aji_E1" : float,
        "fij_E1" : float,
        "S_E1"  : float,
        "Aji_E2" : float,
        "Aji_M1" : float,
        "Aji_M2" : float,
    }

    df["Aji"] = pd.read_csv(file, sep=' ', header=0, dtype=dtype)

    file = "./Kanti_2017_CE.txt"
    dtype = {
        "i" : int,
        "j" : int,
        "4.10" : float,
        "4.30" : float,
        "4.50" : float,
        "4.70" : float,
        "4.90" : float,
        "5.10" : float,
        "5.30" : float,
        "5.50" : float,
        "5.70" : float,
        "5.90" : float,
    }

    df["CE"] = pd.read_csv(file, sep=' ', header=0, dtype=dtype)

    for k1 in ("Aji","CE"):
        for k3 in ("i","j"):
            for k2 in ("Configuration", "Term", "J"):
                df[k1][f"{k3}_{k2}"] = df["Level"][k2].loc[df[k1][k3]].tolist()

    return df

if __name__ == "__main__":


    df = readAll()


    print(df["CE"].head())
    print(df["Aji"].head())
    print(df["Level"].head())

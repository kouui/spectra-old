
if __name__ == "__main__":

    import sys
    sys.path.append("..")

    from src.Atomic import Gaunt

    print(Gaunt.g0(1), Gaunt.g1(1), Gaunt.g2(1))

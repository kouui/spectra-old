

if __name__ == "__main__":

    import sys
    sys.path.append("..")

    from src import AtomCls
    from src.visual import Grotrian

    atom = AtomCls.Atom("/Users/liu/kouui/workspace/spectra/atom/C_III_Be_like.txt")
    gro = Grotrian.Grotrian(_atom=atom)
    gro.make_fig(_figsize=(10,6))
    gro.show_fig()
    #gro.save_fig("./test.png")

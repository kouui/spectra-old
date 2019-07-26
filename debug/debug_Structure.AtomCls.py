
if __name__ == "__main__":

    import sys
    sys.path.append("..")

    from src.Structure import AtomCls

    file = "/Users/liu/kouui/workspace/spectra/atom/C_III_Be_like.txt"
    atom = AtomCls.Atom(file)
    #--- assert that the index - configuration search method works well
    idx_line = 3
    conf_line = atom.line_idx_to_conf(idx_line)
    assert idx_line == atom.conf_to_line_idx(conf_lower=conf_line["lower"], conf_upper=conf_line["upper"])

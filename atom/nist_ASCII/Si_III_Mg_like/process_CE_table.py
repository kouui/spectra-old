
def _set_float(_s, _p, _m):

    assert _s[-3] in (_p,_m), "bad format"

    if _s[-3] == _m:
        _s1 = _s[:-3] + 'e-' + _s[-2:]
    elif _s[-3] == _p:
        _s1 = _s[:-3] + 'e+' + _s[-2:]
    else:
        assert False

    return float(_s1)

if __name__ == "__main__":

    with open("./CE_table_origin.txt", "r") as handle:
        lines = handle.readlines()

    plus  = lines[1].split()[2][-3]
    minus = lines[0].split()[2][-3]

    count = 0
    s_out = ""
    for li in lines:
        words = li.split()

        i = int(words[0]) - 1
        j = int(words[1]) - 1
        s_out += f"{i:<2d}\t{j:<2d}\t"

        for k in range(2,12):
            v =  _set_float(_s=words[k], _p=plus, _m=minus)
            s_out += f"{v:<1.2e}"
            if k == 11:
                s_out += "\n"
            else:
                s_out += "\t"

        count += 1

    print(f"processed : {count} lines")

    with open("./CE_table_processed.txt", "w") as handle:
        handle.write(s_out)

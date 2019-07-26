
ctj_dict = {
    "0" : ("1s2.2s2", "1S", "0"),
    "1" : ("1s2.2s.2p", "3P", "0"),
    "2" : ("1s2.2s.2p", "3P", "1"),
    "3" : ("1s2.2s.2p", "3P", "2"),
    "4" : ("1s2.2s.2p", "1P", "1"),
    "5" : ("1s2.2p2", "3P", "0"),
    "6" : ("1s2.2p2", "3P", "1"),
    "7" : ("1s2.2p2", "3P", "2"),
    "8" : ("1s2.2p2", "1D", "2"),
    "9" : ("1s2.2p2", "1S", "0"),
}


def skip_line(_ln):
    r"""
    skip
        1. comment line. start with '#'
        2. empty line. _ln.strip() is ''
    """
    if _ln[0] == "#" or _ln.strip() == '':
        return True

    return False

def check_end(_ln):
    r"""
    end data reading
        1. _ln starts with "END" end data reading
    """
    if _ln[:3].upper() == "END":
        return True

    return False


def read_data_file(_path):
    r"""

    Parameters
    ----------

    _path : str
        the path to data file

    Returns
    -------

    _result : dict
        a hash dictionary stores necessary information of CE coe

    """

    # output result
    _result = {
        "ctj_i" : [],
        "ctj_j" : [],
        "coe" : [],
        "f1" : [],
        "f2" : [],
        "nLine" : 0
    }

    # read file to list of ... (every line in ASCII file)
    with open(_path, 'r') as _f:
        _tlines = _f.readlines()

    # read line by line
    for kk, _ln in enumerate(_tlines):

        # whether skip line
        if skip_line(_ln):
            continue

        # whether end data reading
        if check_end(_ln):
            break

        _words = _ln.split()
        _i, _j = _words[0].strip(), _words[1].strip()
        _coes = [v.strip() for v in _words[3:-2]]
        _f1 = _words[-2].strip()
        _f2 = _words[-1].strip()

        _result["ctj_i"].append( ctj_dict[_i] )
        _result["ctj_j"].append( ctj_dict[_j] )
        _result["coe"].append( _coes )
        _result["f1"].append( _f1 )
        _result["f2"].append( _f2 )

        _result["nLine"] += 1

    return _result

def format_CE_information(_result, _path, _prefix=""):
    r"""
    given the resotored CE information in dictionary `_result`,
    format and output text file to _path.

    Parameters
    -----------
    _result : dict
        a hash dictionary stores necessary information of CE

    _path : str
        the path to the output file file *.Electron which stores the formatted CE information.

    _prefix : str
        always the configuration string does not contains the configuration of inner shells.
        In these case we have to add the prefix to them, default: ''

    Returns
    -------

    None
    """

    # check lengths of lists
    for _key, _value in _result.items():
        if isinstance(_value, list):
            assert len(_value) == _result["nLine"], "bad length in `_result`."

    # format output string
    _s_out  = "#   "
    _s_out += "{:<15s}".format("configuration_i") + " "*4
    _s_out += "{:<8s}".format("term_i") + " "*4
    _s_out += "{:<4s}".format("J_i") + " "*4
    _s_out += "{:<15s}".format("configuration_j") + " "*4
    _s_out += "{:<8s}".format("term_j") + " "*4
    _s_out += "{:<4s}".format("J_j") + " "*4
    for k in range(9):
        _s_out += "{:<8s}".format("--------") + " "*4
    _s_out += "{:<2s}".format("f1") + " "*4
    _s_out += "{:<2s}".format("f2") + " "*4
    _s_out += "\n"

    for i in range(_result["nLine"]):
        _s_out += "    "
        _s_out += "{:<15s}".format(_prefix+_result["ctj_i"][i][0]) + " "*4
        _s_out += "{:<8s}".format(_result["ctj_i"][i][1]) + " "*4
        _s_out += "{:<4s}".format(_result["ctj_i"][i][2]) + " "*4
        _s_out += "{:<15s}".format(_prefix+_result["ctj_j"][i][0]) + " "*4
        _s_out += "{:<8s}".format(_result["ctj_j"][i][1]) + " "*4
        _s_out += "{:<4s}".format(_result["ctj_j"][i][2]) + " "*4
        for k in range(9):
            _s_out += "{:<8s}".format(_result["coe"][i][k]) + " "*4
        _s_out += "{:<2s}".format(_result["f1"][i]) + " "*4
        _s_out += "{:<2s}".format(_result["f2"][i]) + " "*4

        _s_out += '\n'

    _s_out += '\n'

    with open(_path, "w") as _f:
        _f.write(_s_out)

    return None


if __name__ == "__main__":

    path = "./Berrington_et_al_1985.txt"
    result = read_data_file(_path=path)
    print("restored ", result["nLine"], " Lines.")
    path_out = "../../atom/C_III/Collisional_Excitation/Berrington_et_al_1985.Electron"
    format_CE_information(_result=result, _path=path_out, _prefix="")

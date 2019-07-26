
def is_bad_line(_line):
    r"""
    if
        1. empty line
        2. -------- line
        3. Configuration  | ...  title line : lines start with Alphabet
    return True
    """
    # check -------- line
    if _line[0] == '-':
        return True

    # check empty line
    # if empty line elements in list _ele are all ""
    _ele = [i.strip() for i in _line.split('|')]
    # then check whether whether they are the same --> ""
    if len(set(_ele)) <= 1:
        return True

    # check title line
    if _ele[0][0].isalpha():
        return True

    # if not bad line return False
    return False

def str2num(_s):
    r"""
    convert
        1. string of integer --> int
        2. string of float --> float
    """
    try:
        return int(_s)
    except ValueError:
        return float(_s)
    else :
        return None

def read_Aji_information(_path):
    r"""

    Parameters
    ----------

    _path : str
        the path to NIST ASCII Level file *.NistLine

    Returns
    -------

    _result : dict
        a hash dictionary stores necessary information of Ajis

    """

    # output result
    _result = {
        "ctj_i" : [],
        "ctj_j" : [],
        "Aji[s^-1]" : [],
        "Wavelength[AA]" : [],
        "nLine" : 0
    }

    # read *.NistLine file to list of ... (every line in ASCII file)
    with open(_path, 'r') as _f:
        _tlines = _f.readlines()

    _count = 0
    # read line by line
    for kk, _tline in enumerate(_tlines):

        # check bad line
        if is_bad_line(_tline):
            continue

        _words = _tline.split('|')

        #  get cfj for upper/lower levelss
        _conf_i, _term_i, _J_i = _words[5].strip(), _words[6].strip()[:2], _words[7].strip()
        _conf_j, _term_j, _J_j = _words[8].strip(), _words[9].strip()[:2], _words[10].strip()

        # Aji and wavelength

        _Aji = float(_words[2].strip())
        _wl = float(_words[0].strip())

        # append into _result
        _result["ctj_i"].append( (_conf_i, _term_i, _J_i) )
        _result["ctj_j"].append( (_conf_j, _term_j, _J_j) )
        _result["Aji[s^-1]"].append( _Aji )
        _result["Wavelength[AA]"].append( _wl )

        _count += 1

    _result["nLine"] = _count

    return _result

def format_Aji_information(_result, _path, _prefix=""):
    r"""
    given the resotored Aji information in dictionary `_result`,
    format and output text file to _path.

    Parameters
    -----------
    _result : dict
        a hash dictionary stores necessary information of Ajis

    _path : str
        the path to the output text file *.Aji which stores the formatted Aji information.

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
    _s_out += "{:<40s}".format("configuration_i") + " "*4
    _s_out += "{:<8s}".format("term_i") + " "*4
    _s_out += "{:<4s}".format("J_i") + " "*4
    _s_out += "{:<40s}".format("configuration_j") + " "*4
    _s_out += "{:<8s}".format("term_j") + " "*4
    _s_out += "{:<4s}".format("J_j") + " "*4
    _s_out += "{:<9s}".format("Aji[s^-1]") + " "*4
    _s_out += "{:<14s}".format("Wavelength[AA]") + " "*4
    _s_out += "\n"

    for i in range(_result["nLine"]):
        _s_out += "    "
        _s_out += "{:<40s}".format(_prefix+_result["ctj_i"][i][0]) + " "*4
        _s_out += "{:<8s}".format(_result["ctj_i"][i][1]) + " "*4
        _s_out += "{:<4s}".format(_result["ctj_i"][i][2]) + " "*4
        _s_out += "{:<40s}".format(_prefix+_result["ctj_j"][i][0]) + " "*4
        _s_out += "{:<8s}".format(_result["ctj_j"][i][1]) + " "*4
        _s_out += "{:<4s}".format(_result["ctj_j"][i][2]) + " "*4
        _s_out += "{:<1.3e}".format(_result["Aji[s^-1]"][i]) + " "*4
        _s_out += "{:<1.7e}".format(_result["Wavelength[AA]"][i]) + " "*4

        _s_out += '\n'
    _s_out += '\n'

    with open(_path, "w") as _f:
        _f.write(_s_out)

    return None

if __name__ == "__main__":

    #-------------------------------------------------------------------------
    # C III
    #-------------------------------------------------------------------------
    #path = "../atom/NIST_ASCII/C_III/C_III.NistLine"
    #path_out = "../atom/C_III/Einstein_A/Nist.Aji"
    #-------------------------------------------------------------------------

    #-------------------------------------------------------------------------
    # O V
    #-------------------------------------------------------------------------
    path = "../atom/NIST_ASCII/O_V/O_V.NistLine"
    path_out = "../atom/O_V/Einstein_A/Nist.Aji"
    #-------------------------------------------------------------------------
    result = read_Aji_information(_path=path)
    print("restored ", result["nLine"], " Lines.")
    format_Aji_information(_result=result, _path=path_out, _prefix="")


# a hash dictionary mapping symbolic quantum number L to its integer value

L_s2i = { "S" : 0, "P" : 1, "D" : 2, "F" : 3, "G" : 4, "H" : 5, "I" : 6 }
L_i2s = { 0 : "S", 1 : "P", 2 : "D", 3 : "F", 4 : "G", 5 : "H", 6 : "I" }

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
    if _ele[0].isalpha():
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

def store_level_info(_conf, _term, _J, _g, _E, _2s_plus_1, _L, _stage, _result):

    _result["configuration"].append( _conf )
    _result["term"].append( _term )
    _result["n"].append( _conf.split('.')[-1][0] )
    _result["J"].append( _J )
    _result["L"].append( _L )
    _result["2S+1"].append( _2s_plus_1 )
    _result["g=2J+1"].append( _g )
    _result["stage"].append( _stage )
    _result["E[eV]"].append( _E )

    return None

def read_level_information(_term_ulim, _path):
    r"""
    Nist Level txt is listed in a order of increasing level energy.
    So we select Levels up to a specific term `_term_ulim`.

    Parameters
    ----------

    _term_ulim : tuple of str, (conf, term)
        the upper limit term of our selection

    _path : str
        the path to NIST ASCII Level file *.NistLevel

    Returns
    -------

    _result : dict
        a hash dictionary stores necessary information of Levels
    """

    # output result
    _result = {
        "configuration" : [],
        "term" : [],
        "n" : [],
        "J" : [],
        "L" : [],
        "2S+1" : [],
        "g=2J+1" : [],
        "stage" : [],
        "E[eV]" : [],
        "nLevel" : 0
    }

    # read *.NistLevel file to list of str (every line in ASCII file)
    with open(_path, 'r') as _f:
        _tlines = _f.readlines()

    # read information line by line
    _count = 0                      # count levels we include
    _i = 0                          # count text line we read
    _duplicated = False             # the boolean variable for not double checking bad lines
    while _i < len(_tlines):

        # already checked bad line
        if _duplicated:
            _duplicated = False
            _i += 1
            continue

        _tline = _tlines[_i]

        # skip if bad line
        if is_bad_line(_tline):
            _i += 1
            continue

        # start to extract level information from good lines
        _words = _tline.split('|')
        _conf, _term = _words[0].strip(), _words[1].strip()[:2]

        # if reached the upper limit, break while loop
        if (_conf, _term) == _term_ulim:
            break

        # otherwise, this is the term we need, then store information
        _J, _g, _E = _words[2].strip(), _words[3].strip(), float(_words[4].strip())
        _2s_plus_1 = _term[0]
        _L = str(L_s2i[_term[1]])
        _stage = "-"

        store_level_info(_conf=_conf, _term=_term, _J=_J, _g=_g, _E=_E,
                        _2s_plus_1=_2s_plus_1, _L=_L, _stage=_stage, _result=_result)
        # count of text line +1; count of level +1
        # now _i points to the next line
        _i += 1; _count += 1

        # check whether we have reached to the end before asserting next line
        if _i >= len(_tlines) : break

        # if next text line is not bad line, then it must belong to the same fine structure
        while not is_bad_line(_tlines[_i]):
            _words = _tlines[_i].split('|')
            _J, _g, _E = _words[2].strip(), _words[3].strip(), float(_words[4].strip())
            store_level_info(_conf=_conf, _term=_term, _J=_J, _g=_g, _E=_E,
                                _2s_plus_1=_2s_plus_1, _L=_L, _stage=_stage, _result=_result)
            _i += 1; _count += 1
            if _i >= len(_tlines) : break

        # next line is bad line, this bad line is already checked, so set _duplicated to True
        _duplicated = True

    _result["nLevel"] = _count

    return _result

def format_level_information(_result, _path, _prefix=""):
    r"""
    given the resotored Level information in dictionary `_result`,
    format and output text file to _path.

    Parameters
    -----------
    _result : dict
        a hash dictionary stores necessary information of Levels

    _path : str
        the path to the output text file *.Level which stores the formatted necessary level information.

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
            assert len(_value) == _result["nLevel"], "bad length in `_result`."

    # format output string
    _s_out  = "#   "
    _s_out += "{:<40s}".format("configuration") + " "*4
    _s_out += "{:<8s}".format("term") + " "*4
    _s_out += "{:<4s}".format("J") + " "*4
    _s_out += "{:<4s}".format("n") + " "*4
    _s_out += "{:<4s}".format("L") + " "*4
    _s_out += "{:<4s}".format("2S+1") + " "*4
    _s_out += "{:<6s}".format("g=2J+1") + " "*4
    _s_out += "{:<5s}".format("stage") + " "*4
    _s_out += "{:<14s}".format("E[eV]") + " "*4
    _s_out += '\n'

    for i in range(_result["nLevel"]):
        _s_out += "    "
        _s_out += "{:<40s}".format(_prefix+_result["configuration"][i]) + " "*4
        _s_out += "{:<8s}".format(_result["term"][i]) + " "*4
        _s_out += "{:<4s}".format(_result["J"][i]) + " "*4
        _s_out += "{:<4s}".format(_result["n"][i]) + " "*4
        _s_out += "{:<4s}".format(_result["L"][i]) + " "*4
        _s_out += "{:<4s}".format(_result["2S+1"][i]) + " "*4
        _s_out += "{:<6s}".format(_result["g=2J+1"][i]) + " "*4
        _s_out += "{:<5s}".format(_result["stage"][i]) + " "*4
        _s_out += "{:<1.7e}".format(_result["E[eV]"][i]) + " "*4

        _s_out += '\n'

    _s_out += '\n'

    with open(_path, "w") as _f:
        _f.write(_s_out)

    return None


if __name__ == "__main__":

    # Nist Level txt is listed in a order of increasing level energy

    #-------------------------------------------------------------------------
    # C III
    #-------------------------------------------------------------------------
    #path = "../atom/NIST_ASCII/C_III/C_III.NistLevel"
    #term_ulim = ("1s2.2s.3s","3S")
    #-------------------------------------------------------------------------

    #-------------------------------------------------------------------------
    # O V
    #-------------------------------------------------------------------------
    path = "../atom/NIST_ASCII/O_V/O_V.NistLevel"
    term_ulim = ("1s2.2s.3s","3S")
    #-------------------------------------------------------------------------
    result = read_level_information(_term_ulim=term_ulim, _path=path)
    print("restored ", result["nLevel"], " Levels.")
    path_out = path.split("NIST_ASCII")[0] + "config/" + path.split('/')[-1].replace("Nist", "")
    format_level_information(_result=result, _path=path_out, _prefix="")

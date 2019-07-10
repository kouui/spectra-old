

L_dict = {
    "S" : "0",
    "P" : "1",
    "D" : "2",
    "F" : "3",
}

def is_bad_line(_line):
    r"""
    if
        1. empty line
        2. -------- line
        3. Configuration  | ...  title line
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
    if _line[0] == 'C':
        return True

    # if not bad line return False
    return False

def str2num(_s):
    r"""
    convert string of integer to int and string of float to float
    """
    try:
        return int(_s)
    except ValueError:
        return float(_s)

def store_level_info(_Index, _conf, _term, _J, _g, _E, _2s_plus_1, _L, _stage, result_Level):

    result_Level["Index"].append( _Index )
    result_Level["configuration"].append( _conf )
    result_Level["term"].append( _term )
    result_Level["n"].append( _conf.split('.')[-1][0] )
    result_Level["J"].append( _J )
    result_Level["L"].append( _L )
    result_Level["2S+1"].append( _2s_plus_1 )
    result_Level["g=2J+1"].append( _g )
    result_Level["stage"].append( _stage )
    result_Level["E[eV]"].append( _E )

    return None

def format_level_info(result_level):

    # check lengths of lists
    _length = len(result_level["Index"])
    for key, value in result_level.items():
        assert len(value) == _length

    #output string
    _s_out = "\tIndex\tconfiguration\t\tterm\tn\tJ\tL\t2S+1\tg=2J+1\tstage\tE[eV]\n"
    for i in range(_length):
        _s_out += "\t{0:<2d}\t\t{1:<16s}\t{2:<4s}\t{3:1s}\t{4:3s}\t{5:2s}\t{6:2s}\t\t{7:2s}\t\t{8:1s}\t\t+{9:1.7e}\n".format(
        result_level["Index"][i],result_level["configuration"][i],result_level["term"][i],
        result_level["n"][i], result_level["J"][i],result_level["L"][i],result_level["2S+1"][i],
        result_level["g=2J+1"][i],result_level["stage"][i],result_level["E[eV]"][i]
        )
    _s_out += "\n"

    return _s_out

def read_level_infomation(_levels, _path):

    _fnames = {}
    # read file name
    _spl = _path.split('/')
    if _spl[-1] == '':
        _fnames["Folder"] = _spl[-2]
    else:
        _fnames["Folder"] = '/' + _spl[-1]
    _fnames["Level"] = _path + _fnames["Folder"] + ".Level"

    # levels --> configurations
    #_configurations = [i.split()[0] for i in _levels]
    #_terms = [i.split()[1] for i in _levels]

    # read .Level file to list of lines
    with open(_fnames["Level"],"r") as _f:
        _lines = _f.readlines()

    # output result
    _result = {
        "Level" : {
            "Index" : [],
            "configuration" : [],
            "term" : [],
            "n" : [],
            "J" : [],
            "L" : [],
            "2S+1" : [],
            "g=2J+1" : [],
            "stage" : [],
            "E[eV]" : []
        },
        "Line" : {
            "i" : [],
            "j" : [],
            "Aji[s^-1]" : [],
            "Wavelength[AA]" : []
        }
    }

    # read information line by line
    _count = 0
    i = 0
    _duplicated = False             # the boolean variable for not double checking bad lines
    while i < len(_lines):
        # already checked
        if _duplicated:
            _duplicated = False
            i += 1
            continue

        _line = _lines[i]
        #  skip bad lines
        if is_bad_line(_line):
            i += 1
            continue

        # standard case
        _words = _line.split('|')
        _conf, _term = _words[0].strip(), _words[1].strip()[:2]
        # if not the level we need
        if (_conf + ' ' + _term) not in _levels:
            i += 1
            continue
        # if it is the level we need, store infomation
        _J, _g, _E = _words[2].strip(), _words[3].strip(), float(_words[4].strip())
        _2s_plus_1 = _term[0]
        _L = L_dict[_term[1]]
        _stage = "-"
        store_level_info(_Index=_count, _conf=_conf, _term=_term, _J=_J, _g=_g, _E=_E,
                            _2s_plus_1=_2s_plus_1, _L=_L, _stage=_stage, result_Level=_result["Level"])
        i += 1; _count += 1
        if i >= len(_lines) : break
        # next line is not bad line
        while not is_bad_line(_lines[i]):
            _words = _lines[i].split('|')
            _J, _g, _E = _words[2].strip(), _words[3].strip(), float(_words[4].strip())
            store_level_info(_Index=_count, _conf=_conf, _term=_term, _J=_J, _g=_g, _E=_E,
                                _2s_plus_1=_2s_plus_1, _L=_L, _stage=_stage, result_Level=_result["Level"])
            i += 1; _count += 1
            if i >= len(_lines) : break
        # next line is bad line, this bad line is already checked, so set _duplicated to True
        _duplicated = True

    return _fnames, _result

def is_bad_line2(_line):
    r"""
    if
        1. empty line
        2. -------- line
        3.    Ritz  | ...  title line
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
    if _line.split('|')[0].strip()[0].isalpha():
        return True

    # if not bad line return False
    return False

def format_line_info(result_line):

    # check lengths of lists
    _length = len(result_line["i"])
    for key, value in result_line.items():
        assert len(value) == _length

    # output string
    _s_out = "\t i\t\t j\t\tWavelength[AA]\t\tAji[s^-1]\n"
    for _i in range(_length):
        _s_out += "\t{0:2d}\t\t{1:2d}\t\t{2:1.8e}\t\t{3:1.4e}\n".format(
        result_line["i"][_i], result_line["j"][_i], result_line["Wavelength[AA]"][_i], result_line["Aji[s^-1]"][_i]
        )
    _s_out += "\n"

    return _s_out

def read_line_information(_fnames, _result):
    r"""
    """
    # fname of .Line file
    _fnames["Line"] = _fnames["Level"][:-5] + "Line"

    # read .Level file to list of lines
    with open(_fnames["Line"],"r") as _f:
        _lines = _f.readlines()

    # make check list for comparison
    _idxs = _result["Level"]["Index"]
    _lvls = []
    for k in _idxs:
        _lvls.append( (_result["Level"]["configuration"][k], _result["Level"]["term"][k], _result["Level"]["J"][k]) )

    # check line by line
    _As = []
    _idx_pair = []
    for kk, _line in enumerate(_lines):

        # check bad line
        if is_bad_line2(_line):
            continue

        # check whether conf and Term are that we need
        _words = _line.split('|')
        _conf_i, _term_i, _J_i = _words[5].strip(), _words[6].strip()[:2], _words[7].strip()
        _conf_j, _term_j, _J_j = _words[8].strip(), _words[9].strip()[:2], _words[10].strip()
        if ( (_conf_i, _term_i, _J_i) in _lvls ) and ( (_conf_j, _term_j, _J_j) in _lvls ):
            # there are several transition not inside .Line file (forbidden line)
            _i = _lvls.index( (_conf_i, _term_i, _J_i) )
            _j = _lvls.index( (_conf_j, _term_j, _J_j) )
            _idx_pair.append( (_i, _j) )
            #_w = 1.23984176 * 1E+4 / ( _result["Level"]["E[eV]"][_j] - _result["Level"]["E[eV]"][_i] ) # [AA]
            _As.append( float(_words[2].strip()) )

    # sort and supplement transitions whose Aji equals to 0
    for _i in _idxs:
        for _j in range(_i+1, _idxs[-1]+1):
            _result["Line"]["i"].append( _i )
            _result["Line"]["j"].append( _j )

            _w = 1.23984176 * 1E+4 / ( _result["Level"]["E[eV]"][_j] - _result["Level"]["E[eV]"][_i] ) # [AA]
            _result["Line"]["Wavelength[AA]"].append( _w )

            if (_i, _j) in _idx_pair:
                ii = _idx_pair.index( (_i, _j) )
                _result["Line"]["Aji[s^-1]"].append( _As[ii] )
            else:
                _result["Line"]["Aji[s^-1]"].append( 0.0 )

    return _result

def make_result_file(t_levels, t_path):

    fnames, result = read_level_infomation(_levels=t_levels, _path=t_path)
    result = read_line_information(_fnames=fnames, _result=result)
    s_out = format_level_info(result_level=result["Level"])
    s_out += format_line_info(result_line=result["Line"])
    #print(s_out)

    # output .result file
    fnames["Result"] = fnames["Level"][:-5] + "Result"
    with open(fnames["Result"], "w") as f:
        f.write(s_out)

    return None


if __name__ == "__main__":

    # column order :
    #       Level : Configuration   Term    J   g   Level(eV)
    #       Line :  i   j   Aji[s^-1]   Wavelength[AA]

    levels = ["1s2.2s2 1S", "1s2.2s.2p 3P", "1s2.2s.2p 1P", "1s2.2p2 3P", "1s2.2p2 1D", "1s2.2p2 1S"]
    path = "/Users/liu/kouui/workspace/statistical_equilibrium/atom/nist_ASCII/C_III_Be_like/"
    make_result_file(t_levels=levels, t_path=path)

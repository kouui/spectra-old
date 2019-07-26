
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


def read_general_info(_rs, _lns):
    r"""
    read general information in *.Level config file
        1. Title
        2. Z
        3. Element
    """

    for _i, _ln in enumerate(_lns[_rs:]):

        if skip_line(_ln):
            continue
        if check_end(_ln):
            break

        if _ln.split()[0].strip().lower() == "title:":
            _Title = ' '.join( _ln.split()[1:] )

        if _ln.split()[0].strip().lower() == "z":
            _Z = _ln.split()[-1].strip()

        if _ln.split()[0].strip().lower() == "element":
            _Element = _ln.split()[-1].strip()

        if _ln.split()[0].strip().lower() == "nlevel":
            _nLevel = int( _ln.split()[-1].strip() )

    _re = _rs + _i + 1

    return _re, _Title, _Z, _Element, _nLevel

def read_level_info(_rs, _lns, _Level_info, _erg, _g, _stage):
    r"""
    read level information to
        1. _Level_info
            - ["configuration"]
            - ["term"]
            - ["J"]
            - ["2S+1"]
        2. _erg
        3. _g
        4. _stage
    """

    _idx = 0
    for _i, _ln in enumerate(_lns[_rs:]):

        if skip_line(_ln):
            continue
        elif check_end(_ln):
            break

        _words = _ln.split()
        _words = [_v.strip() for _v in _words]

        _Level_info["configuration"].append( _words[0] )
        _Level_info["term"].append( _words[1] )
        _Level_info["J"].append( _words[2] )
        _Level_info["2S+1"].append( _words[5] )

        _erg[_idx] = float( _words[8] )
        _g[_idx] = int( _words[6] )
        _stage[_idx] = int( _words[7] )

        _idx += 1

    _re = _rs + _i + 1

    return _re

def read_line_info(_lns, _Aji, _line_ctj_table):
    r"""
    read line information
    """
    _count = 0
    for _i, _ln in enumerate(_lns[:]):

        if skip_line(_ln):
            continue
        elif check_end(_ln):
            break

        _words = _ln.split()
        _words = [_v.strip() for _v in _words]

        # get ctj pair
        ctj_ij = ( (_words[0],_words[1],_words[2]), (_words[3],_words[4],_words[5]) )
        # get line_index
        """
        try :
            line_index = _line_ctj_table.index( ctj_ij )
        except ValueError:
            continue
        else:
            assert False, "Error in function `AtomIO.read_line_info()``"

        _Aji[line_index] = float( _words[6] )
        _w0_AA[line_index] = float( _words[7] )

        _count += 1

        if _count == _Aji.size:
            break
        """
        if ctj_ij in _line_ctj_table:
            line_index = _line_ctj_table.index( ctj_ij )
            _Aji[line_index] += float( _words[6] )
            _count += 1

    return None

def read_CE_Temperature(_lns):
    r"""
    read Temperature grid for interpolation
    """
    for _i, _ln in enumerate(_lns[:]):

        if skip_line(_ln):
            continue
        elif check_end(_ln):
            break

        _words = _ln.split()
        _words = [_v.strip() for _v in _words]

        if _words[0].lower() == "type":
            _type = _words[1]

        if _words[0].lower() == "temperature":
            _Te = [float(v) for v in _words[1:]]

    _nTe = len( _Te )
    _re = _i + 1

    return _re, _nTe, _Te, _type

def read_CE_table(_rs, _lns, _CE_table, _f1, _f2, _line_ctj_table):
    r"""
    read CE table for interpolation
    """
    _count = 0
    for _i, _ln in enumerate(_lns[_rs:]):

        if skip_line(_ln):
            continue
        elif check_end(_ln):
            break

        _words = _ln.split()
        _words = [_v.strip() for _v in _words]

        # get ctj pair
        _ctj_ij = ( (_words[0],_words[1],_words[2]), (_words[3],_words[4],_words[5]) )

        if _ctj_ij in _line_ctj_table:
            line_index = _line_ctj_table.index( _ctj_ij )
            _CE_table[line_index,:] += [float(v) for v in _words[6:-2]]
            _f1[line_index] = float(_words[-2])
            _f2[line_index] = float(_words[-1])

            _count += 1

    return None


def skip_line(_ln):

    if _ln[0] == "#":
        return True

    return False

def check_end(_ln):

    if _ln[:3] == "END":
        return True

    return False


def read_title(_rs, _lns):

    _title = ""
    for _i, _ln in enumerate(_lns[_rs:]):

        if skip_line(_ln):
            continue

        _title += _ln + '\n'
        break

    _re = _rs + _i + 1

    return _re, _title

def read_general_info(_rs, _lns):

    for _i, _ln in enumerate(_lns[_rs:]):

        if skip_line(_ln):
            continue
        if check_end(_ln):
            break

        if _ln.split(':')[0].strip() == "Z":
            _Z = _ln.split(':')[-1].strip()

        if _ln.split(':')[0].strip() == "Element":
            _Element = _ln.split(':')[-1].strip()

        if _ln.split(':')[0].strip() == "nLevel":
            _nLevel = int( _ln.split(':')[-1].strip() )

        if _ln.split(':')[0].strip() == "nLine":
            _nLine = int( _ln.split(':')[-1].strip() )

    _re = _rs + _i + 1

    return _re, _Z, _Element, _nLevel, _nLine


def read_level_info(_rs, _lns, _Level_info, _erg, _g, _stage):

    _idx = 0
    for _i, _ln in enumerate(_lns[_rs:]):

        if skip_line(_ln):
            continue
        elif check_end(_ln):
            break

        _words = _ln.split()
        _words = [_v.strip() for _v in _words]
        assert _idx == int(_words[0]), "bad indexing."

        _Level_info["configuration"].append( _words[1] )
        _Level_info["term"].append( _words[2] )
        _Level_info["J"].append( _words[4] )
        _Level_info["2S+1"].append( _words[6] )

        _erg[_idx] = float( _words[9] )
        _g[_idx] = int( _words[7] )
        _stage[_idx] = int( _words[8] )

        _idx += 1

    _re = _rs + _i + 1

    return _re

def read_line_info(_rs, _lns, _idxI, _idxJ, _Aji, _w0_AA):

    _idx = 0
    for _i, _ln in enumerate(_lns[_rs:]):

        if skip_line(_ln):
            continue
        elif check_end(_ln):
            break

        _words = _ln.split()
        _words = [_v.strip() for _v in _words]

        _idxI[_idx] = int( _words[0] )
        _idxJ[_idx] = int( _words[1] )
        _Aji[_idx] = float( _words[3] )
        _w0_AA[_idx] = float( _words[2] )

        _idx += 1

    _re = _rs + _i + 1

    return _re

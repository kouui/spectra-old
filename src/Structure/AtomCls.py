
import numpy as np
from .. import Constants as Cst
from . import AtomIO as AIO

class Atom:

    def __init__(self, filepath):
        r"""
        initial method of class Atom
        """

        self.filepath = filepath
        self.__read_data_file(filepath)


    def __read_data_file(self, filepath):
        r"""
        restore

            self.Level
            self.Level_info
            self.Level_info_table
            self.Line
            self.Line_idx_table
            self.CE_Te_table
            self.CE_table
            self.CE_type
            self.CE_coe

        from the filepath.
        """

        with open(filepath, 'r') as file:
            fLines = file.readlines()

        #--- read general info
        rs, self.title = AIO.read_title(_rs=0, _lns=fLines)
        rs, self.Z, self.Element, self.nLevel, self.nLine, self.paramCE = AIO.read_general_info(_rs=rs, _lns=fLines)

        #--- read Level info
        dtype  = np.dtype([
                          ('erg',np.double),            #: level energy, erg
                          ('g',np.uint8),               #: g=2J+1, statistical weight
                          ('stage',np.uint8),           #: ionization stage
                          ])
        self.Level = np.recarray(self.nLevel, dtype=dtype)
        self.Level_info = {"configuration" : [], "term" : [], "J": [], "2S+1": []}
        rs = AIO.read_level_info(rs, _lns=fLines, _Level_info=self.Level_info,
                            _erg=self.Level.erg[:], _g=self.Level.g[:], _stage=self.Level.stage[:])
        self.Level.erg[:] *= Cst.eV2erg_

        #--- make tuple of tuple (configuration, term, J)
        self.Level_info_table = []
        for k in range(self.nLevel):
            self.Level_info_table.append((self.Level_info["configuration"][k],
                                          self.Level_info["term"][k],
                                          self.Level_info["J"][k]))
        self.Level_info_table = tuple(self.Level_info_table)

        #--- read line info
        dtype = np.dtype([('idxI',np.uint16),           #: level index, the Level index of lower level
                           ('idxJ',np.uint16),          #: level index, the Level index of lower level
                           ('AJI',np.double),           #: Einstein Aji coefficient
                           ('f0',np.double),            #: central frequency
                           ('w0',np.double),            #: central wavelength in cm
                           ('w0_AA',np.double),         #: central wavelength in Angstrom
                           ])
        self.Line = np.recarray(self.nLine, dtype=dtype)
        rs = AIO.read_line_info(_rs=rs, _lns=fLines, _idxI=self.Line.idxI[:], _idxJ=self.Line.idxJ[:],
                            _Aji=self.Line.AJI[:], _w0_AA=self.Line.w0_AA[:])
        self.Line.w0[:] = self.Line.w0_AA[:] * 1E-8
        self.Line.f0[:] = Cst.c_ / self.Line.w0[:]

        #--- make tuple of tuple (self.Line.idxI, self.Line.idxJ)
        self.Line_idx_table = []
        for k in range(self.nLine):
            self.Line_idx_table.append( (self.Line.idxI[k], self.Line.idxJ[k]) )
        self.Line_idx_table = tuple(self.Line_idx_table)

        #--- read Collision Excitation
        nL, nTemp = self.paramCE
        self.CE_Te_table = np.empty(nTemp, dtype=np.double)   #: Temperature mesh [:math:`K`]
        self.CE_table = np.empty((nL,nTemp), dtype=np.double)
        dtype  = np.dtype([
                          ('idxI',np.uint8),      #: level index, the Level index of lower level
                          ('idxJ',np.uint8),      #: level index, the Level index of lower level
                          ('f1',np.uint8),
                          ('f2',np.uint8),
                          ('gi',np.uint8),        #: statistical weight of lower level
                          ('gj',np.uint8),        #: statistical weight of upper level
                          ('dEij',np.double)      #: excitation energy, [:math:`erg`]
                          ])
        self.CE_coe = np.recarray(nL, dtype=dtype)
        self.CE_type = []
        rs = AIO.read_CE_info(_rs=rs, _lns=fLines, _idxI=self.CE_coe.idxI[:], _idxJ=self.CE_coe.idxJ[:],
                        _Te=self.CE_Te_table[:], _table=self.CE_table[:,:], _CE_type=self.CE_type,
                        _f1=self.CE_coe.f1[:], _f2=self.CE_coe.f2[:])

        for k in range(nL):
            self.CE_coe.gi[k] = self.Level.g[self.CE_coe.idxI[k]]
            self.CE_coe.gj[k] = self.Level.g[self.CE_coe.idxJ[k]]
            self.CE_coe.dEij[k] = self.Level.erg[self.CE_coe.idxJ[k]] - self.Level.erg[self.CE_coe.idxI[k]]

    def line_idx_to_conf(self, line_idx):
        r"""
        given the index of a line transition,
        returnthe  configuration and term information for lower and upper levels, respectively.
        """
        idxI = self.Line.idxI[line_idx]
        idxJ = self.Line.idxJ[line_idx]

        conf = {
            "lower" : "",
            "upper" : ""
        }
        for key, idx in zip( ("lower","upper"), (idxI, idxJ) ):
            conf[key] = (self.Level_info["configuration"][idx],
                         self.Level_info["term"][idx],
                         self.Level_info["J"][idx])

        return conf

    def conf_to_level_idx(self, conf, term, J):
        r"""
        given the configuration, term and J of a specific level,
        return the index of that level.
        """

        return self.Level_info_table.index( (conf, term, J) )

    def idxIJ_to_line_idx(self, idxI, idxJ):
        r"""
        given the index of lower and upper levels, respectively,
        return the index of that line.
        """

        return self.Line_idx_table.index( (idxI, idxJ) )

    def conf_to_line_idx(self, conf_lower, conf_upper):
        r"""
        given the configuration tuple of lower and upper level, respectively,
        return the index of that transition
        """
        idxI = self.Level_info_table.index( conf_lower )
        idxJ = self.Level_info_table.index( conf_upper )

        return self.idxIJ_to_line_idx(idxI, idxJ)

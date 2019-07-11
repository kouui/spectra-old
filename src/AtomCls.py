
import numpy as np
import Constants as Cst
import AtomIO as AIO

class Atom:

    def __init__(self, filepath):

        self.filepath = filepath
        self.__read_data_file(filepath)


    def __read_data_file(self, filepath):

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

        #--- read Collision Excitation
        nL, nTemp = self.paramCE
        self.CE_Te = np.empty(nTemp, dtype=np.double)   #: Temperature mesh [:math:`K`]
        self.CE_table = np.empty((nL,nTemp), dtype=np.double)
        dtype  = np.dtype([
                          ('idxI',np.uint8),      #: level index, the Level index of lower level
                          ('idxJ',np.uint8),      #: level index, the Level index of lower level
                          ('f1',np.uint8),
                          ('f2',np.uint8),
                          ])
        self.CE_coe = np.recarray(nL, dtype=dtype)
        self.CE_type = []
        rs = AIO.read_CE_info(_rs=rs, _lns=fLines, _idxI=self.CE_coe.idxI[:], _idxJ=self.CE_coe.idxJ[:],
                        _Te=self.CE_Te[:], _table=self.CE_table[:,:], _CE_type=self.CE_type,
                        _f1=self.CE_coe.f1[:], _f2=self.CE_coe.f2[:])







if __name__ == "__main__":

    file = "/Users/liu/kouui/workspace/statistical_equilibrium/atom/C_III_Be_like.txt"
    atom = Atom(file)

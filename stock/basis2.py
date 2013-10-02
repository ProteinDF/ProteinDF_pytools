#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy

from basisset import *

class Basis2(object):
    """
    """
    def __init__(self):
        self._basisset = {}
        self._basisset_J = {}
        self._basisset_XC = {}

    def load(self, path):
        """
        basis2ファイルを読み込む
        """
        name = ''
        cgto_types = []
        num_CGTOs = 0
        cgto_id = 0
        num_PGTOs = 0
        pgto_id = 0
        CGTO = None
        AUX_MODE = False
        AUX_XC = False
        
        f = open(path, 'r')
        for line in f:
            line = line.strip()
            if len(line) == 0:
                continue
            if (line[0] == '#') or (line[0:2] == '//'):
                continue

            print(line)
            
            if len(name) == 0:
                name = line
                if name[0] == 'A':
                    AUX_MODE = True
                    AUX_XC = False
                else:
                    AUX_MODE = False
                continue

            if num_CGTOs == 0:
                CGTOs = line.split()
                num_s = int(CGTOs[0])
                num_p = int(CGTOs[1])
                num_d = int(CGTOs[2])
                cgto_types = ['s'] * num_s + ['p'] * num_p + ['d'] * num_d
                num_CGTOs = num_s + num_p + num_d
                cgto_id = 0
                continue
            
            if num_PGTOs == 0:
                num_PGTOs = int(line)
                CGTO = ContractedGTO(cgto_types[cgto_id],
                                            num_PGTOs)
                continue
                
            values = line.split()
            exp = values[0]
            coef = 1.0
            if len(values) == 2:
                coef = values[1]
            CGTO[pgto_id] = PrimitiveGTO(exp, coef)
            pgto_id = pgto_id +1
            
            if pgto_id == num_PGTOs:
                num_PGTOs = 0
                pgto_id = 0
                cgto_id = cgto_id + 1

            if cgto_id == num_CGTOs:
                num_CGTOs = 0
                cgto_id = 0

                if AUX_MODE == False:
                    self._basisset[name] = copy.deepcopy(CGTO)
                    name = ''
                else:
                    if AUX_XC == False:
                        self._basisset_J[name] = copy.deepcopy(CGTO)
                        AUX_XC = True
                    else:
                        self._basisset_XC[name] = copy.deepcopy(CGTO)
                        name = ''
                
        f.close()

    def get_basisset(self, name):
        """
        指定された名前のbasissetを返す
        """
        answer = self._basisset.get(name, BasisSet())
        return answer

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()

    #bs = Basis2()
    #bs.load("/home/hirano/local/intel/ProteinDF/data/basis2")

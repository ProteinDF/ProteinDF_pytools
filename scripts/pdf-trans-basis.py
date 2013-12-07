#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2002-2014 The ProteinDF project
# see also AUTHORS and README.
# 
# This file is part of ProteinDF.
# 
# ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import os
import re
import argparse
from orbinfo import *

class BasisSetParser(object):
    """
    gaussian94 format text data (eg. ESML html) parser
    """
    _re_atom_line = re.compile('\s*(\S+)\s+0')
    _re_CGTO_line = re.compile('\s*(\S+)\s+(\S+)\s+(\S+)')
    _re_PGTO_line = re.compile('\s*(\S+)\s+(\S+)')
    
    def __init__(self):
        self._basissets = {}
        
    def load(self, file_path, mode = 'Gaussian94'):
        self._load_gaussian94(file_path)

    def _load_gaussian94(self, file_path):
        if (os.path.isfile(file_path) != True):
            return

        fin = open(file_path, "r")
        basissets = []
        atom_basis = []
        in_section = False
        for line in fin:
            line = line.rstrip('\n')
            
            if (line[0:4] == '****'):
                if in_section:
                    # parse & store data
                    basissets.append(atom_basis)
                    atom_basis = []
                else:
                    in_section = True
                    atom_basis = []
                continue

            if in_section:
                atom_basis.append(line)
        fin.close()

        for atom_basis in basissets:
            (atom, bs) = self._parse_block_gaussian94(atom_basis)
            self._basissets[atom] = bs

    def _parse_block_aussian94(self, lines):
        has_error = False
        
        # [atom] 0
        matchObj = self._re_atom_line.match(lines[0])
        atom = matchObj.group(1)
        #print(">>>> %s" % (atom))
        
        line_index = 1
        cgtos = []
        while line_index < len(lines):
            # CGTO
            matchObj = self._re_CGTO_line.match(lines[line_index])
            type = matchObj.group(1)
            num_of_PGTOs = int(matchObj.group(2))
            CGTO_coef = float(matchObj.group(3))
            line_index += 1
            cgto = ContractedGTO(type.lower(), num_of_PGTOs)
            
            for pgto_index in range(num_of_PGTOs):
                pgto_line = lines[line_index + pgto_index]
                matchObj = self._re_PGTO_line.match(pgto_line)
                exponent = float(matchObj.group(1))
                coef = float(matchObj.group(2))
                cgto[pgto_index] = PrimitiveGTO(exponent, coef)
            line_index += num_of_PGTOs
            cgtos.append(cgto)
            
        basisset = BasisSet('', len(cgtos))
        for cgto_index in range(len(cgtos)):
            basisset[cgto_index] = cgtos[cgto_index]
            
        return (atom, basisset)

    def __str__(self):
        answer = ''
        for atom, basisset in self._basissets.items():
            answer += ">>>> %s\n" % (atom)
            answer += str(basisset)
            answer += '\n'

        return answer
    
def main():
    parser = argparse.ArgumentParser(description='parse EMSL BasisSet Data')
    parser.add_argument('FILE',
                        nargs=1,
                        help='EMSL basisset html file')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default = False)
    args = parser.parse_args()

    # setting
    file_path = args.FILE[0]
    verbose = args.verbose
    
    obj = BasisSetParser()
    if verbose:
        print('reading: %s' % (file_path))
    obj.load(file_path)
    print(obj)
    

if __name__ == '__main__':
    main()


#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2014 The ProteinDF development team.
# see also AUTHORS and README if provided.
# 
# This file is a part of the ProteinDF software package.
# 
# The ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# The ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import os
import re
import argparse

import pdfpytools as pdf

class BasisSetParser(object):
    """
    gaussian94 format text data (eg. ESML html) parser
    """
    _re_atom_line = re.compile('\s*(\S+)\s+0')
    _re_CGTO_line = re.compile('\s*(\S+)\s+(\S+)\s+(\S+)')
    _re_PGTO_line1 = re.compile('\s*(\S+)\s+(\S+)')
    _re_PGTO_line2 = re.compile('\s*(\S+)\s+(\S+)\s+(\S+)')
    
    def __init__(self):
        pass
        
    def load(self, file_path, mode = 'Gaussian94'):
        with open(file_path, 'r') as f:
            lines = f.readlines()
            self._parse(lines)
        
        self._load_gaussian94(file_path)

    def parse(self, lines):
        lines = lines.splitlines()
        basissets, basissets_density, basissets_xc = self._parse_gaussian94(lines)
        
        return (basissets, basissets_density, basissets_xc)
        
    def _parse_gaussian94(self, lines):
        basissets = {}
        basissets_density = {}
        basissets_xc = {}
        
        basisset_blocks = []
        atom_block = []
        for line in lines:
            line = line.rstrip('\n')
            
            if (line[0:4] == '****'):
                # parse & store data
                if len(atom_block) > 0:
                    basisset_blocks.append(atom_block)
                atom_block = []
            else:
                atom_block.append(line)

        for atom_block in basisset_blocks:
            (atom, bs) = self._parse_block_gaussian94(atom_block)
            if atom != None:
                if atom not in basissets:
                    basissets[atom] = bs
                elif atom not in basissets_density:
                    basissets_density[atom] = bs
                elif atom not in basissets_xc:
                    basissets_xc[atom] = bs
                else:
                    sys.stderr.write('too many setup basis: {}\n{}'.format(atom, repr(bs)))

        return (basissets, basissets_density, basissets_xc)

    def _parse_block_gaussian94(self, lines):
        # [atom] 0
        matchObj = self._re_atom_line.match(lines[0])
        if matchObj == None:
            return (None, None)
        atom = matchObj.group(1)
        atom = atom.capitalize()
        
        line_index = 1
        cgtos = []
        while line_index < len(lines):
            # CGTO
            matchObj = self._re_CGTO_line.match(lines[line_index])
            cgto_type = matchObj.group(1).lower()
            num_of_PGTOs = int(matchObj.group(2))
            CGTO_coef = float(matchObj.group(3))
            line_index += 1

            if len(cgto_type) == 1:
                cgto = pdf.ContractedGTO(cgto_type, num_of_PGTOs)
            
                for pgto_index in range(num_of_PGTOs):
                    pgto_line = lines[line_index + pgto_index]
                    matchObj = self._re_PGTO_line1.match(pgto_line)
                    exponent = float(matchObj.group(1).replace('D', 'E'))
                    coef = float(matchObj.group(2).replace('D', 'E'))
                    cgto[pgto_index] = pdf.PrimitiveGTO(exponent, coef)
                line_index += num_of_PGTOs
                cgtos.append(cgto)
            elif len(cgto_type) == 2:
                cgto1 = pdf.ContractedGTO(cgto_type[0], num_of_PGTOs)
                cgto2 = pdf.ContractedGTO(cgto_type[1], num_of_PGTOs)
                for pgto_index in range(num_of_PGTOs):
                    pgto_line = lines[line_index + pgto_index]
                    matchObj = self._re_PGTO_line2.match(pgto_line)
                    exponent = float(matchObj.group(1).replace('D', 'E'))
                    coef1 = float(matchObj.group(2).replace('D', 'E'))
                    coef2 = float(matchObj.group(3).replace('D', 'E'))
                    cgto1[pgto_index] = pdf.PrimitiveGTO(exponent, coef1)
                    cgto2[pgto_index] = pdf.PrimitiveGTO(exponent, coef2)
                line_index += num_of_PGTOs
                cgtos.append(cgto1)
                cgtos.append(cgto2)
            else:
                raise    
                
        basisset = pdf.BasisSet('', len(cgtos))
        for cgto_index in range(len(cgtos)):
            basisset[cgto_index] = cgtos[cgto_index]

        return (atom, basisset)

    #def __str__(self):
    #    answer = ''
    #    for atom, basisset in self._basissets.items():
    #        answer += ">>>> %s\n" % (atom)
    #        answer += str(basisset)
    #        answer += '\n'
    #
    #    return answer
    
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


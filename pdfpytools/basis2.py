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
import logging

import pdfbridge as bridge
import pdfpytools as pdf

class Basis2(object):
    _shell_chars = list("spdfghijklmnopqr")
    _data = None

    def __new__(cls, *args, **kwargs):
        '''for singleton'''
        if '_inst' not in vars(cls):
            cls._inst = super(Basis2, cls).__new__(cls, *args, **kwargs)
        return cls._inst

    def __init__(self):
        nullHandler = bridge.NullHandler()
        self._logger = logging.getLogger(__name__)
        self._logger.addHandler(nullHandler)

        self._initialize()
        
    def _load(self):
        db_path = '%s/data/basis2' % (pdf.pdf_home())

        self._data = {}
        self._data['basis'] = {}
        self._data['basis_j'] = {}
        self._data['basis_xc'] = {}
        
        self._line_count = 0
        basisset = None
        name = None
        #SPD = None
        #SPD_order = 0
        cgto = None
        numOfCGTPs = 0
        numOfPGTOs = None
        PGTO_index = 0
        with open(db_path) as f:
            for line in f:
                self._line_count += 1
                line = line.strip()
                if (len(line) == 0):
                    continue
                if (line[0] == '#'):
                    continue

                if (name == None):
                    self._logger.debug('%6d: %s' % (self._line_count, line))

                    name = line
                    if (name[0] == 'O'):
                        #basisset = copy.deepcopy(self._load_basis(f))
                        basisset = pdf.BasisSet(self._load_basis(f))

                        self._logger.debug(str(basisset))
                        #print('>>>> {0}'.format(name))
                        #print(basisset)
                        #print('')
                        
                        basisset = basisset.expand()
                        basisset.name = name
                        if name in self._data['basis']:
                            self._logger.warning('duplicate basisset name: {0}'.format(name))
                        self._data['basis'][name] = basisset
                        name = None
                    elif (name[0] == 'A'):
                        #basisset = copy.deepcopy(self._load_basis_aux(f))
                        basisset = pdf.BasisSet(self._load_basis_aux(f))
                        basisset = basisset.expand()
                        basisset.name = name
                        self._data['basis_j'][name] = basisset
                        
                        #basisset = copy.deepcopy(self._load_basis_aux(f))
                        basisset = pdf.BasisSet(self._load_basis_aux(f))
                        basisset = basisset.expand()
                        basisset.name = name
                        self._data['basis_xc'][name] = basisset
                        name = None
                    else:
                        self._logger.critical('basisset name is not found: {}'.format(name))
                        raise
                    
    def _load_basis(self, fh):
        basisset = None
        SPDFG = None
        SPDFG_order = 0
        cgto = None
        numOfPGTOs = None
        PGTO_index = 0
        
        for line in fh:
            self._line_count += 1
            line = line.strip()
            if (len(line) == 0):
                continue
            if (line[0] == '#'):
                continue

            #self._logger.debug('%6d: %s' % (self._line_count, line))
            if (SPDFG == None):
                SPDFG = []
                input_SPDFG = line.split()
                if len(input_SPDFG) > len(self._shell_chars):
                    self._logger.warning(
                        "not support shell: {} > {}".format(len(input_SPDFG),
                                                            len(self._shell_chars)))
                numOfCGTOs = 0
                for i in range(len(input_SPDFG)):
                    value = int(input_SPDFG[i])
                    SPDFG.append(value)
                    numOfCGTOs += value
                basisset = pdf.BasisSet("", numOfCGTOs)
            elif (numOfPGTOs == None):
                numOfPGTOs = int(line)
                shell_type = ''

                order = 0
                for i in range(len(SPDFG)):
                    if (order <= SPDFG_order) and (SPDFG_order < order + SPDFG[i]):
                        shell_type = self._shell_chars[i]
                        break
                    order += SPDFG[i]
                
                # self._logger.debug("create CGTO SPD=%d shell_type=%s size=%d" % (SPD_order, shell_type, numOfPGTOs))
                cgto = pdf.ContractedGTO(shell_type, numOfPGTOs)
                assert(cgto.shell_type == shell_type)
                assert(len(cgto) == numOfPGTOs)
                PGTO_index = 0
            else:
                values = line.split()
                assert(len(values) == 2)
                exponent = float(values[0])
                coefficient = float(values[1])
                cgto[PGTO_index] = pdf.PrimitiveGTO(exponent, coefficient)
                PGTO_index += 1
                if (PGTO_index == numOfPGTOs):
                    # end of PGTO list
                    assert(len(cgto) == numOfPGTOs)
                    basisset[SPDFG_order] = cgto
                    SPDFG_order += 1
                    numOfPGTOs = None
                if (SPDFG_order == numOfCGTOs):
                    # end of basis set
                    assert(basisset.get_num_of_CGTOs('s') == SPDFG[0])
                    assert(basisset.get_num_of_CGTOs('p') == SPDFG[1])
                    assert(basisset.get_num_of_CGTOs('d') == SPDFG[2])
                    assert(basisset.get_num_of_CGTOs('f') == SPDFG[3])
                    assert(basisset.get_num_of_CGTOs('g') == SPDFG[4])
                    assert(len(basisset) == numOfCGTOs)
                    break

        return basisset
                    
    def _load_basis_aux(self, fh):
        basisset = None
        SPDFG = None
        SPDFG_order = 0
        cgto = None
        numOfPGTOs = None
        PGTO_index = 0
        
        for line in fh:
            self._line_count += 1
            line = line.strip()
            if (len(line) == 0):
                continue
            if (line[0] == '#'):
                continue

            #self._logger.debug('%6d: %s' % (self._line_count, line))
            if (SPDFG == None):
                SPDFG = []
                input_SPDFG = line.split()
                if len(input_SPDFG) > len(self._shell_chars):
                    self._logger.warning(
                        "not support shell: {} > {}".format(len(input_SPDFG),
                                                            len(self._shell_chars)))
                numOfCGTOs = 0
                for i in range(len(input_SPDFG)):
                    value = int(input_SPDFG[i])
                    SPDFG.append(value)
                    numOfCGTOs += value
                basisset = pdf.BasisSet("", numOfCGTOs)
            elif (numOfPGTOs == None):
                numOfPGTOs = int(line)
                shell_type = ''

                order = 0
                for i in range(len(SPDFG)):
                    if (order <= SPDFG_order) and (SPDFG_order < order + SPDFG[i]):
                        shell_type = self._shell_chars[i]
                        break
                    order += SPDFG[i]
                
                #print("create CGTO SPD=%d shell_type=%s size=%d" % (SPD_order, shell_type, numOfPGTOs))
                cgto = pdf.ContractedGTO(shell_type, numOfPGTOs)
                assert(cgto.shell_type == shell_type)
                assert(len(cgto) == numOfPGTOs)
                PGTO_index = 0
            else:
                values = line.split()
                assert(len(values) > 0)
                exponent = float(values[0])
                coefficient = 1.0
                #print("PGTO_index=%d" % (PGTO_index))
                cgto[PGTO_index] = pdf.PrimitiveGTO(exponent, coefficient)
                PGTO_index += 1
                if (PGTO_index >= numOfPGTOs):
                    # print('>>>> end of pGTO: spd=%d' % (SPD_order))
                    # end of PGTO list
                    # basisset[SPDFG_order] = copy.deepcopy(cgto)
                    basisset[SPDFG_order] = pdf.ContractedGTO(cgto)
                    SPDFG_order += 1
                    numOfPGTOs = None
                if (SPDFG_order >= numOfCGTOs):
                    # end of basis set
                    break

        return basisset
                        
    def get_basisset(self, name):
        return self._data['basis'].get(name, pdf.BasisSet())

    def get_basisset_j(self, name):
        return self._data['basis_j'].get(name, pdf.BasisSet())
        
    def get_basisset_xc(self, name):
        return self._data['basis_xc'].get(name, pdf.BasisSet())

    #@property
    #def basis(self):
    #    return self._basis

    #@property
    #def basis_j(self):
    #    return self._basis_j
        
    #@property
    #def basis_xc(self):
    #    return self._basis_xc

    def get_basis2(self):
        output = ''
        # basis
        for name, basis in self._data['basis'].items():
            output += '# {}\n'.format(name)
            output += str(basis)
            output += '\n'

        # basis for J, XC
        for name, basis_j in self._data['basis_j'].items():
            if name in self._data['basis_xc']:
                basis_xc = self._basis_xc[name]
                output += name + '\n'
                output += str(basis_j)
                output += '#\n'
                output += str(basis_xc)
                output += '\n'
            
        return output
    
    def _initialize(self):
        if self._data == None:
            self._load()
        
        
if __name__ == "__main__":
    #logging.basicConfig(level=logging.DEBUG)
    bs2 = pdf.Basis2()
    print(bs2.get_basis2())

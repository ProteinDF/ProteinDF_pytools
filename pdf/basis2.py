#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import re
import logging

import bridge
from basisset import *
from pdfcommon import *

class Basis2(object):
    def __init__(self):
        nullHandler = bridge.NullHandler()
        self._logger = logging.getLogger(__name__)
        self._logger.addHandler(nullHandler)

        self._initialize()
        
    def _load(self):
        db_path = '%s/data/basis2' % (pdf_home())
        #self._basisset_db = {}

        self._line_count = 0
        basisset = None
        name = None
        SPD = None
        SPD_order = 0
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
                        basisset = copy.deepcopy(self._load_basis(f))

                        self._logger.debug(str(basisset))
                        #print('>>>> {0}'.format(name))
                        #print(basisset)
                        #print('')
                        
                        basisset = basisset.expand()
                        basisset.name = name
                        if name in self._basis:
                            self._logger.warning('duplicate basisset name: {0}'.format(name))
                        self._basis[name] = basisset
                        name = None
                    elif (name[0] == 'A'):
                        basisset = copy.deepcopy(self._load_basis_aux(f))
                        basisset = basisset.expand()
                        basisset.name = name
                        self._basis_j[name] = basisset
                        
                        basisset = copy.deepcopy(self._load_basis_aux(f))
                        basisset = basisset.expand()
                        basisset.name = name
                        self._basis_xc[name] = basisset
                        name = None
                    else:
                        print(name)
                        abort()
                    
    def _load_basis(self, fh):
        basisset = None
        SPD = None
        SPD_order = 0
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

            if (SPD == None):
                SPD = [0, 0, 0]
                input_SPD = line.split()
                numOfCGTOs = 0
                for i in range(len(input_SPD)):
                    value = int(input_SPD[i])
                    SPD[i] = value
                    numOfCGTOs += value
                basisset = BasisSet("", numOfCGTOs)
            elif (numOfPGTOs == None):
                numOfPGTOs = int(line)
                shell_type = ''
                if (SPD_order < SPD[0]):
                    shell_type = 's'
                elif (SPD_order < (SPD[0] + SPD[1])):
                    shell_type = 'p'
                elif (SPD_order < (SPD[0] + SPD[1] + SPD[2])):
                    shell_type = 'd'
                else:
                    abort()
                #self._logger.debug("create CGTO SPD=%d shell_type=%s size=%d" % (SPD_order, shell_type, numOfPGTOs))
                cgto = ContractedGTO(shell_type, numOfPGTOs)
                PGTO_index = 0
            else:
                values = line.split()
                assert(len(values) == 2)
                exponent = float(values[0])
                coefficient = float(values[1])
                cgto[PGTO_index] = PrimitiveGTO(exponent, coefficient)
                PGTO_index += 1
                if (PGTO_index == numOfPGTOs):
                    # end of PGTO list
                    assert(len(cgto) == numOfPGTOs)
                    basisset[SPD_order] = cgto
                    SPD_order += 1
                    numOfPGTOs = None
                if (SPD_order == numOfCGTOs):
                    # end of basis set
                    assert(basisset.get_num_of_CGTOs('s') == SPD[0])
                    assert(basisset.get_num_of_CGTOs('p') == SPD[1])
                    assert(basisset.get_num_of_CGTOs('d') == SPD[2])
                    assert(len(basisset) == numOfCGTOs)
                    break

        return basisset
                    
    def _load_basis_aux(self, fh):
        basisset = None
        SPD = None
        SPD_order = 0
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

            self._logger.debug('%6d: %s' % (self._line_count, line))
                
            if (SPD == None):
                SPD = [0, 0, 0, 0]
                input_SPD = line.split()
                numOfCGTOs = 0
                for i in range(len(input_SPD)):
                    value = int(input_SPD[i])
                    SPD[i] = value
                    numOfCGTOs += value
                basisset = BasisSet("", numOfCGTOs)
            elif (numOfPGTOs == None):
                numOfPGTOs = int(line)
                shell_type = ''
                if (SPD_order < SPD[0]):
                    shell_type = 's'
                elif (SPD_order < (SPD[0] + SPD[1])):
                    shell_type = 'p'
                elif (SPD_order < (SPD[0] + SPD[1] + SPD[2])):
                    shell_type = 'd'
                elif (SPD_order < (SPD[0] + SPD[1] + SPD[2] + SPD[3])):
                    shell_type = 'spd'
                else:
                    abort()
                #print("create CGTO SPD=%d shell_type=%s size=%d" % (SPD_order, shell_type, numOfPGTOs))
                cgto = ContractedGTO(shell_type, numOfPGTOs)
                PGTO_index = 0
            else:
                values = line.split()
                assert(len(values) == 1)
                exponent = float(values[0])
                coefficient = 1.0
                #print("PGTO_index=%d" % (PGTO_index))
                cgto[PGTO_index] = PrimitiveGTO(exponent, coefficient)
                PGTO_index += 1
                if (PGTO_index >= numOfPGTOs):
                    # print('>>>> end of pGTO: spd=%d' % (SPD_order))
                    # end of PGTO list
                    basisset[SPD_order] = copy.deepcopy(cgto)
                    SPD_order += 1
                    numOfPGTOs = None
                if (SPD_order >= numOfCGTOs):
                    # end of basis set
                    break

        return basisset
                        
    def get_basisset(self, name):
        return self._basis.get(name, BasisSet())

    def get_basisset_j(self, name):
        return self._basis_j.get(name, BasisSet())
        
    def get_basisset_xc(self, name):
        return self._basis_xc.get(name, BasisSet())

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
        for name, basis in self.basis.items():
            output += name + '\n'
            output += str(basis)
            output += '\n'

        # basis for J, K
        for name, basis_j in self.basis_j.items():
            if name in self.basis_k:
                basis_k = self.basis_k[name]
                output += name + '\n'
                output += str(basis_j)
                output += '#\n'
                output += str(basis_k)
                output += '\n'
            
        return output
    
    def _initialize(self):
        self._basis = {}
        self._basis_j = {}
        self._basis_xc = {}

        self._load()
        
if __name__ == "__main__":
    bs2 = Basis2()
    #bs2.debug = True
    bs2.load()

    print(bs2.get_basis2())

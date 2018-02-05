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
logger = logging.getLogger(__name__)

import proteindf_bridge as bridge
from .basisset import BasisSet, ContractedGTO, PrimitiveGTO


class Basis2(object):
    """
    read ProteinDF basis2 file

    >>> bs2 = Basis2("data/basis2.sample")
    >>> print(bs2.get_basis2())
    O-CARBON (621/41) by FS
     3 2 0 0 0 0 0 0 0 0 0
     6
         2.808064e+03  2.017830e-03
         4.211383e+02  1.543320e-02
         9.558662e+01  7.558155e-02
         2.673900e+01  2.478282e-01
         8.432827e+00  4.793725e-01
         2.760582e+00  3.338344e-01
     2
         5.447004e+00 -7.784077e-02
         4.792422e-01  5.689560e-01
     1
         1.461565e-01  1.000000e+00
     4
         1.813085e+01  1.585473e-02
         4.099883e+00  9.568277e-02
         1.185837e+00  3.049119e-01
         3.685974e-01  4.935017e-01
     1
         1.097200e-01  1.000000e+00
    <BLANKLINE>
    O-6-31G.H
     2 0 0 0 0 0 0 0 0 0 0
     3
         1.873114e+01  3.349460e-02
         2.825394e+00  2.347269e-01
         6.401217e-01  8.137573e-01
     1
         1.612778e-01  1.000000e+00
    <BLANKLINE>
    O-6-31G.He
     2 0 0 0 0 0 0 0 0 0 0
     3
         3.842163e+01  2.376600e-02
         5.778030e+00  1.546790e-01
         1.241774e+00  4.696300e-01
     1
         2.979640e-01  1.000000e+00
    <BLANKLINE>
    O-6-31G.C
     3 2 0 0 0 0 0 0 0 0 0
     6
         3.047525e+03  1.834700e-03
         4.573695e+02  1.403730e-02
         1.039487e+02  6.884260e-02
         2.921016e+01  2.321844e-01
         9.286663e+00  4.679413e-01
         3.163927e+00  3.623120e-01
     3
         7.868272e+00 -1.193324e-01
         1.881288e+00 -1.608542e-01
         5.442493e-01  1.143456e+00
     1
         1.687144e-01  1.000000e+00
     3
         7.868272e+00  6.899910e-02
         1.881288e+00  3.164240e-01
         5.442493e-01  7.443083e-01
     1
         1.687144e-01  1.000000e+00
    <BLANKLINE>
    O-6-31G.N
     3 2 0 0 0 0 0 0 0 0 0
     6
         4.173511e+03  1.834800e-03
         6.274579e+02  1.399500e-02
         1.429021e+02  6.858700e-02
         4.023433e+01  2.322410e-01
         1.282021e+01  4.690700e-01
         4.390437e+00  3.604550e-01
     3
         1.162636e+01 -1.149610e-01
         2.716280e+00 -1.691180e-01
         7.722180e-01  1.145852e+00
     1
         2.120313e-01  1.000000e+00
     3
         1.162636e+01  6.758000e-02
         2.716280e+00  3.239070e-01
         7.722180e-01  7.408950e-01
     1
         2.120313e-01  1.000000e+00
    <BLANKLINE>
    O-6-31G.O
     3 2 0 0 0 0 0 0 0 0 0
     6
         5.484672e+03  1.831100e-03
         8.252350e+02  1.395010e-02
         1.880470e+02  6.844510e-02
         5.296450e+01  2.327143e-01
         1.689757e+01  4.701930e-01
         5.799635e+00  3.585209e-01
     3
         1.553962e+01 -1.107775e-01
         3.599934e+00 -1.480263e-01
         1.013762e+00  1.130767e+00
     1
         2.700058e-01  1.000000e+00
     3
         1.553962e+01  7.087430e-02
         3.599934e+00  3.397528e-01
         1.013762e+00  7.271586e-01
     1
         2.700058e-01  1.000000e+00
    <BLANKLINE>
    O-DZVP2.H
     2 1 0 0 0 0 0 0 0 0 0
     4
         5.099918e+01  9.660500e-03
         7.483218e+00  7.372890e-02
         1.777468e+00  2.958581e-01
         5.193295e-01  7.159053e-01
     1
         1.541100e-01  1.000000e+00
     1
         7.500000e-01  1.000000e+00
    <BLANKLINE>
    O-DZVP2.C
     3 2 1 0 0 0 0 0 0 0 0
     7
         5.784157e+03  8.190000e-04
         8.693035e+02  6.293500e-03
         1.985116e+02  3.178120e-02
         5.642990e+01  1.172734e-01
         1.828546e+01  3.034763e-01
         6.448715e+00  4.535214e-01
         2.341860e+00  2.430591e-01
     2
         5.459533e+00 -7.780440e-02
         4.781968e-01  5.714947e-01
     1
         1.457301e-01  1.000000e+00
     5
         3.425856e+01  5.804300e-03
         7.863895e+00  4.064030e-02
         2.344519e+00  1.550219e-01
         7.961715e-01  3.531444e-01
         2.726804e-01  4.550062e-01
     1
         8.926050e-02  1.000000e+00
     1
         6.000000e-01  1.000000e+00
    <BLANKLINE>
    O-DZVP2.N
     3 2 1 0 0 0 0 0 0 0 0
     7
         8.104176e+03  7.969000e-04
         1.217314e+03  6.128900e-03
         2.777399e+02  3.104710e-02
         7.884760e+01  1.153682e-01
         2.553716e+01  3.025738e-01
         9.004571e+00  4.557913e-01
         3.283528e+00  2.430208e-01
     2
         7.849357e+00 -7.763640e-02
         6.862239e-01  5.679815e-01
     1
         2.035026e-01  1.000000e+00
     5
         4.901461e+01 -5.900700e-03
         1.131667e+01 -4.164440e-02
         3.403405e+00 -1.610249e-01
         1.161111e+00 -3.583538e-01
         3.953358e-01 -4.488415e-01
     1
         1.268981e-01  1.000000e+00
     1
         7.000000e-01  1.000000e+00
    <BLANKLINE>
    O-DZVP2.O
     3 2 1 0 0 0 0 0 0 0 0
     7
         1.081440e+04  7.809000e-04
         1.623753e+03  6.010200e-03
         3.701827e+02  3.052220e-02
         1.049748e+02  1.140089e-01
         3.398442e+01  3.019574e-01
         1.198431e+01  4.571107e-01
         4.385970e+00  2.432478e-01
     2
         1.063003e+01 -7.876540e-02
         9.398526e-01  5.706303e-01
     1
         2.766213e-01  1.000000e+00
     5
         6.154422e+01  6.623800e-03
         1.427619e+01  4.646420e-02
         4.331768e+00  1.744229e-01
         1.476604e+00  3.666115e-01
         4.959857e-01  4.369361e-01
     1
         1.544836e-01  1.000000e+00
     1
         8.000000e-01  1.000000e+00
    <BLANKLINE>
    A-CARBON (7/2;7/2) by FS
     7 2 0 0 0 0 0 0 0 0 0
     1
         3.467611e+02  1.000000e+00
     1
         7.966888e+01  1.000000e+00
     1
         2.201073e+01  1.000000e+00
     1
         6.473901e+00  1.000000e+00
     1
         2.508147e+00  1.000000e+00
     1
         7.781041e-01  1.000000e+00
     1
         2.288027e-01  1.000000e+00
     1
         2.508147e+00  1.000000e+00
     1
         2.288027e-01  1.000000e+00
    A-CARBON (7/2;7/2) by FS
     7 2 0 0 0 0 0 0 0 0 0
     1
         1.160000e+02  1.000000e+00
     1
         2.655629e+01  1.000000e+00
     1
         7.336909e+00  1.000000e+00
     1
         2.157967e+00  1.000000e+00
     1
         8.360489e-01  1.000000e+00
     1
         2.593680e-01  1.000000e+00
     1
         7.626760e-02  1.000000e+00
     1
         8.360489e-01  1.000000e+00
     1
         7.626760e-02  1.000000e+00
    <BLANKLINE>
    A-DZVP2.H
     4 1 1 0 0 0 0 0 0 0 0
     1
         4.500000e+01  1.000000e+00
     1
         7.500000e+00  1.000000e+00
     1
         3.000000e-01  1.000000e+00
     1
         1.500000e+00  1.000000e+00
     1
         1.500000e+00  1.000000e+00
     1
         1.500000e+00  1.000000e+00
    A-DZVP2.H
     4 1 1 0 0 0 0 0 0 0 0
     1
         1.500000e+01  1.000000e+00
     1
         2.500000e+00  1.000000e+00
     1
         1.000000e-01  1.000000e+00
     1
         5.000000e-01  1.000000e+00
     1
         5.000000e-01  1.000000e+00
     1
         5.000000e-01  1.000000e+00
    <BLANKLINE>
    A-DZVP2.C
     8 4 4 0 0 0 0 0 0 0 0
     1
         1.500000e+03  1.000000e+00
     1
         3.300000e+02  1.000000e+00
     1
         9.432000e+01  1.000000e+00
     1
         2.700000e+01  1.000000e+00
     1
         9.920000e+00  1.000000e+00
     1
         2.200000e+00  1.000000e+00
     1
         6.300000e-01  1.000000e+00
     1
         1.800000e-01  1.000000e+00
     1
         9.920000e+00  1.000000e+00
     1
         2.200000e+00  1.000000e+00
     1
         6.300000e-01  1.000000e+00
     1
         1.800000e-01  1.000000e+00
     1
         9.920000e+00  1.000000e+00
     1
         2.200000e+00  1.000000e+00
     1
         6.300000e-01  1.000000e+00
     1
         1.800000e-01  1.000000e+00
    A-DZVP2.C
     8 4 4 0 0 0 0 0 0 0 0
     1
         5.000000e+02  1.000000e+00
     1
         1.100000e+02  1.000000e+00
     1
         3.140000e+01  1.000000e+00
     1
         9.000000e+00  1.000000e+00
     1
         3.300000e+00  1.000000e+00
     1
         7.300000e-01  1.000000e+00
     1
         2.100000e-01  1.000000e+00
     1
         6.000000e-02  1.000000e+00
     1
         3.300000e+00  1.000000e+00
     1
         7.300000e-01  1.000000e+00
     1
         2.100000e-01  1.000000e+00
     1
         6.000000e-02  1.000000e+00
     1
         3.300000e+00  1.000000e+00
     1
         7.300000e-01  1.000000e+00
     1
         2.100000e-01  1.000000e+00
     1
         6.000000e-02  1.000000e+00
    <BLANKLINE>
    A-DZVP2.N
     8 4 4 0 0 0 0 0 0 0 0
     1
         2.066000e+03  1.000000e+00
     1
         4.590000e+02  1.000000e+00
     1
         1.310000e+02  1.000000e+00
     1
         3.750000e+01  1.000000e+00
     1
         1.380000e+01  1.000000e+00
     1
         3.060000e+00  1.000000e+00
     1
         8.800000e-01  1.000000e+00
     1
         2.500000e-01  1.000000e+00
     1
         1.380000e+01  1.000000e+00
     1
         3.060000e+00  1.000000e+00
     1
         8.800000e-01  1.000000e+00
     1
         2.500000e-01  1.000000e+00
     1
         1.380000e+01  1.000000e+00
     1
         3.060000e+00  1.000000e+00
     1
         8.800000e-01  1.000000e+00
     1
         2.500000e-01  1.000000e+00
    A-DZVP2.N
     8 4 4 0 0 0 0 0 0 0 0
     1
         6.880000e+02  1.000000e+00
     1
         1.530000e+02  1.000000e+00
     1
         4.370000e+01  1.000000e+00
     1
         1.250000e+01  1.000000e+00
     1
         4.600000e+00  1.000000e+00
     1
         1.020000e+00  1.000000e+00
     1
         2.900000e-01  1.000000e+00
     1
         8.300000e-02  1.000000e+00
     1
         4.600000e+00  1.000000e+00
     1
         1.020000e+00  1.000000e+00
     1
         2.900000e-01  1.000000e+00
     1
         8.300000e-02  1.000000e+00
     1
         4.600000e+00  1.000000e+00
     1
         1.020000e+00  1.000000e+00
     1
         2.900000e-01  1.000000e+00
     1
         8.300000e-02  1.000000e+00
    <BLANKLINE>
    A-DZVP2.O
     8 4 4 0 0 0 0 0 0 0 0
     1
         2.566000e+03  1.000000e+00
     1
         5.700000e+02  1.000000e+00
     1
         1.630000e+02  1.000000e+00
     1
         4.650000e+01  1.000000e+00
     1
         1.700000e+01  1.000000e+00
     1
         3.800000e+00  1.000000e+00
     1
         1.080000e+00  1.000000e+00
     1
         3.100000e-01  1.000000e+00
     1
         1.700000e+01  1.000000e+00
     1
         3.800000e+00  1.000000e+00
     1
         1.080000e+00  1.000000e+00
     1
         3.100000e-01  1.000000e+00
     1
         1.700000e+01  1.000000e+00
     1
         3.800000e+00  1.000000e+00
     1
         1.080000e+00  1.000000e+00
     1
         3.100000e-01  1.000000e+00
    A-DZVP2.O
     8 4 4 0 0 0 0 0 0 0 0
     1
         8.550000e+02  1.000000e+00
     1
         1.900000e+02  1.000000e+00
     1
         5.400000e+01  1.000000e+00
     1
         1.550000e+01  1.000000e+00
     1
         5.660000e+00  1.000000e+00
     1
         1.270000e+00  1.000000e+00
     1
         3.600000e-01  1.000000e+00
     1
         1.000000e-01  1.000000e+00
     1
         5.660000e+00  1.000000e+00
     1
         1.270000e+00  1.000000e+00
     1
         3.600000e-01  1.000000e+00
     1
         1.000000e-01  1.000000e+00
     1
         5.660000e+00  1.000000e+00
     1
         1.270000e+00  1.000000e+00
     1
         3.600000e-01  1.000000e+00
     1
         1.000000e-01  1.000000e+00
    <BLANKLINE>
    <BLANKLINE>
    """
    _shell_chars = list("spdfghijklmnopqr")
    _data = None


    #def __new__(cls, *args, **kwargs):
    #    '''for singleton'''
    #    if '_inst' not in vars(cls):
    #        cls._inst = super(Basis2, cls).__new__(cls, *args, **kwargs)
    #    return cls._inst


    def __init__(self, basis2_path=''):
        self._initialize(basis2_path)


    def _load(self):
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
        logger.debug("load basis2 file: {}".format(self._basis2_path))
        with open(self._basis2_path) as f:
            for line in f:
                self._line_count += 1
                line = line.strip()
                if (len(line) == 0):
                    continue
                if (line[0] == '#'):
                    continue

                if (name == None):
                    logger.debug('%6d: %s' % (self._line_count, line))

                    name = line
                    if (name[0] == 'O'):
                        #basisset = copy.deepcopy(self._load_basis(f))
                        basisset = BasisSet(self._load_basis(f))

                        logger.debug(str(basisset))
                        #print('>>>> {0}'.format(name))
                        #print(basisset)
                        #print('')

                        basisset = basisset.expand()
                        basisset.name = name
                        if name in self._data['basis']:
                            logger.warning('duplicate basisset name: {0}'.format(name))
                        self._data['basis'][name] = basisset
                        name = None
                    elif (name[0] == 'A'):
                        #basisset = copy.deepcopy(self._load_basis_aux(f))
                        basisset = BasisSet(self._load_basis_aux(f))
                        basisset = basisset.expand()
                        basisset.name = name
                        self._data['basis_j'][name] = basisset

                        #basisset = copy.deepcopy(self._load_basis_aux(f))
                        basisset = BasisSet(self._load_basis_aux(f))
                        basisset = basisset.expand()
                        basisset.name = name
                        self._data['basis_xc'][name] = basisset
                        name = None
                    else:
                        logger.critical('basisset name is not found: {}'.format(name))
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

            #logger.debug('%6d: %s' % (self._line_count, line))
            if (SPDFG == None):
                SPDFG = []
                input_SPDFG = line.split()
                if len(input_SPDFG) > len(self._shell_chars):
                    logger.warning(
                        "not support shell: {} > {}".format(len(input_SPDFG),
                                                            len(self._shell_chars)))
                numOfCGTOs = 0
                for i in range(len(input_SPDFG)):
                    value = int(input_SPDFG[i])
                    SPDFG.append(value)
                    numOfCGTOs += value
                SPDFG.extend([0] * (len(self._shell_chars) - len(input_SPDFG)))
                assert(len(SPDFG) == len(self._shell_chars))
                basisset = BasisSet("", numOfCGTOs)
            elif (numOfPGTOs == None):
                numOfPGTOs = int(line)
                shell_type = ''

                order = 0
                for i in range(len(SPDFG)):
                    if (order <= SPDFG_order) and (SPDFG_order < order + SPDFG[i]):
                        shell_type = self._shell_chars[i]
                        break
                    order += SPDFG[i]

                # logger.debug("create CGTO SPD=%d shell_type=%s size=%d" % (SPD_order, shell_type, numOfPGTOs))
                cgto = ContractedGTO(shell_type, numOfPGTOs)
                assert(cgto.shell_type == shell_type)
                assert(len(cgto) == numOfPGTOs)
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

            #logger.debug('%6d: %s' % (self._line_count, line))
            if (SPDFG == None):
                SPDFG = []
                input_SPDFG = line.split()
                if len(input_SPDFG) > len(self._shell_chars):
                    logger.warning(
                        "not support shell: {} > {}".format(len(input_SPDFG),
                                                            len(self._shell_chars)))
                numOfCGTOs = 0
                for i in range(len(input_SPDFG)):
                    value = int(input_SPDFG[i])
                    SPDFG.append(value)
                    numOfCGTOs += value
                basisset = BasisSet("", numOfCGTOs)
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
                cgto = ContractedGTO(shell_type, numOfPGTOs)
                assert(cgto.shell_type == shell_type)
                assert(len(cgto) == numOfPGTOs)
                PGTO_index = 0
            else:
                values = line.split()
                assert(len(values) > 0)
                exponent = float(values[0])
                coefficient = 1.0
                #print("PGTO_index=%d" % (PGTO_index))
                cgto[PGTO_index] = PrimitiveGTO(exponent, coefficient)
                PGTO_index += 1
                if (PGTO_index >= numOfPGTOs):
                    # print('>>>> end of pGTO: spd=%d' % (SPD_order))
                    # end of PGTO list
                    # basisset[SPDFG_order] = copy.deepcopy(cgto)
                    basisset[SPDFG_order] = ContractedGTO(cgto)
                    SPDFG_order += 1
                    numOfPGTOs = None
                if (SPDFG_order >= numOfCGTOs):
                    # end of basis set
                    break

        return basisset

    def get_basisset(self, name):
        return self._data['basis'].get(name, BasisSet())

    def get_basisset_j(self, name):
        return self._data['basis_j'].get(name, BasisSet())

    def get_basisset_xc(self, name):
        return self._data['basis_xc'].get(name, BasisSet())

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
            #output += '# {}\n'.format(name)
            output += str(basis)
            output += '\n'

        # basis for J, XC
        for name, basis_j in self._data['basis_j'].items():
            if name in self._data['basis_xc']:
                basis_xc = self._data['basis_xc'][name]
                output += str(basis_j)
                output += str(basis_xc)
                output += '\n'

        return output


    def _initialize(self, basis2_path):
        if len(basis2_path) == 0:
            from .pdfcommon import pdf_home
            self._basis2_path = '%s/data/basis2' % (pdf_home())
        else:
            self._basis2_path = basis2_path

        if self._data == None:
            self._load()


if __name__ == "__main__":
    import doctest
    doctest.testmod()

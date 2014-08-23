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
import sys
import argparse
import logging
import math
try:
    import msgpack
except:
    import msgpack_pure as msgpack

import bridge
import pdf


def neutralize(atomgroup):
    modeler = bridge.Modeling()
    
    name = atomgroup.name
    name = name.upper()
    if name == 'GLU':
        counters = modeler.neutralize_GLU(atomgroup)
        atomgroup |= counters
    elif name == 'ASP':
        counters = modeler.neutralize_ASP(atomgroup)
        atomgroup |= counters
    elif name == 'LYS':
        counters = modeler.neutralize_LYS(atomgroup)
        atomgroup |= counters
    elif name == 'ARG':
        counters = modeler.neutralize_ARG(atomgroup)
        atomgroup |= counters

    # N-term
    if atomgroup.has_atom('N'):
        N = atomgroup['N']
        H_selector = bridge.Select_Atom('H')
        range_selector = bridge.Select_Range(N.xyz, 1.2)
        H_atoms = atomgroup.select(H_selector).select(range_selector)
        if H_atoms.get_number_of_atoms() == 3:
            counters = modeler.neutralize_Nterm(atomgroup)
            atomgroup |= counters
    
    # C-term
    if atomgroup.has_atom('OXT'):
        counters = modeler.neutralize_Cterm(atomgroup)
        atomgroup |= counters
        
    for key, subgroup in atomgroup.groups():
        subgroup = neutralize(subgroup)
        atomgroup[key] = subgroup

    return atomgroup

def main():
    parser = argparse.ArgumentParser(description='neutralize')
    parser.add_argument('input_brdfile',
                        nargs=1,
                        help='input Bridge file')
    parser.add_argument('output_brdfile',
                        nargs=1,
                        help='output Bridge file')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    input_brd_path = args.input_brdfile[0]
    output_brd_path = args.output_brdfile[0]

    # load
    if verbose:
        print("reading: {}".format(input_brd_path))
    infile = open(input_brd_path, 'rb')
    brddata = msgpack.unpackb(infile.read())
    infile.close()

    models = bridge.AtomGroup(brddata)
    models = neutralize(models)

    #print(models)

    # save
    if verbose:
        print('writing %s.\n'.format(outout_brd_path))
    outfile = open(output_brd_path, 'wb')
    outdata = models.get_dict_data()
    outfile.write(msgpack.packb(outdata))
    outfile.close()
    

if __name__ == '__main__':
    main()

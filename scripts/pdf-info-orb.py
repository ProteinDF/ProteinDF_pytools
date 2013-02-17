#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
output info of orb
"""

import os
import sys
import argparse
try:
    import msgpack
except:
    import msgpack_pure as msgpack
    
import bridge
import pdf

def main():
    # parse args
    parser = argparse.ArgumentParser(description='output info of orb')
    parser.add_argument('orb_index',
                        nargs='+',
                        type=int,
                        help='orbital index (start=0)')
    parser.add_argument('-d', '--db',
                        nargs=1,
                        action='store',
                        default='pdfresults.db',
                        help='ProteinDF results file')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    parser.add_argument('-D', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
        
    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    orb_index = args.orb_index
    
    entry = pdf.PdfArchive(args.db)
    orb_info = pdf.OrbInfo(entry)

    num_of_orbitals = orb_info.get_num_of_orbitals()
    for i in orb_index:
        if i < num_of_orbitals:
            print('%dth: %d %s %s' % (i,
                                      orb_info.get_atom_id(i),
                                      orb_info.get_shell_type(i),
                                      orb_info.get_basis_type(i)))
    
    
if __name__ == '__main__':
    main()
    

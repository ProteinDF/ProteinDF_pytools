#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import proteindf_bridge as bridge
import proteindf_tools as pdf

def find_atom_index(pdfparam, atomgroup):
    answer = set()
    for subgrp_key, subgrp in atomgroup.groups():
        subgrp_atoms = find_atom_index(pdfparam, subgrp)
        answer = answer | subgrp_atoms
    for atom_key, atom in atomgroup.atoms():
        index = pdfparam.find_atom_index(atom)
        answer.add(index)

    return answer

def main():
    # parse args
    parser = argparse.ArgumentParser(description='output atom index')

    parser.add_argument('brd_path',
                        nargs=1,
                        help='find atom list')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-d', '--db',
                       nargs='?',
                       action='store',
                       const='pdfresults.db',
                       help='ProteinDF results file')
    group.add_argument('-p', '--param',
                       nargs='?',
                       action='store',
                       const='pdfparam.mpac',
                       help='ProteinDF parameter file')

    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=[""],
                        help='output as msgpack format')

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    # setting
    brd_path = args.brd_path[0]
    output_path = args.output[0]
    verbose = args.verbose

    if args.db:
        entry = pdf.PdfArchive(args.db)
        atomgroup = entry.get_molecule()
    elif args.param:
        pdfparam = pdf.load_pdfparam(args.param)
        atomgroup = pdfparam.molecule

    # load brd file
    find_atoms = bridge.load_atomgroup(brd_path)
    find_atoms *= pdf.ANG2AU # angstrom -> a.u.
    # print(find_atoms)

    output_data = []
    for subgrp_key, subgrp in find_atoms.groups():
        index_list = list(find_atom_index(pdfparam, subgrp))
        output_data.append(index_list)
    if len(output_path) > 0:
        bridge.save_msgpack(output_data, output_path)
    else:
        # output to stdout
        #for subgrp_key, subgrp in find_atoms.groups():
        #    index_list = find_atom_index(pdfparam, subgrp)
        #    print(index_list)
        print(output_data)


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import math


import proteindf_bridge as bridge
import proteindf_tools as pdf


def main():
    # parse args
    parser = argparse.ArgumentParser(description='calc cell size')
    parser.add_argument('FILE',
                        nargs=1,
                        help='bridge file')
    parser.add_argument('basisset_name',
                        nargs=1,
                        type=str,
                        default="DZVP2",
                        help='basisset name')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    # setting
    mpac_file_path = args.FILE[0]
    basisset_name = args.basisset_name[0]
    verbose = args.verbose

    # reading
    if (verbose == True):
        print("reading: {}}".format(mpac_file_path))
        print("basisset: {}".format(basisset_name))

    # prepare atomgroup
    atom_group = bridge.load_atomgroup(mpac_file_path)

    basis2 = pdf.Basis2()
    number_of_AOs = get_number_of_AOs(basis2, basisset_name, atom_group)

    # output
    print(number_of_AOs)


def get_number_of_AOs(basis2, basisset_name, atom_group):
    number_of_AOs = 0
    for key, subgrp in atom_group.groups():
        number_of_AOs += get_number_of_AOs(basis2, basisset_name, subgrp)
    for key, atom in atom_group.atoms():
        bs_name = "O-{}.{}".format(basisset_name, atom.symbol)
        bs = basis2.get_basisset(bs_name)
        number_of_AOs += bs.get_number_of_AOs()
    return number_of_AOs


if __name__ == '__main__':
    main()

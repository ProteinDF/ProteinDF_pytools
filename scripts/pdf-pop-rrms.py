#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
# import pandas as pd
import csv
try:
    import msgpack
except:
    import msgpack_pure as msgpack
import math
import copy

import pdfbridge as bridge
import pdfpytools as pdf

def load_charges(csv_path):
    print(csv_path)
    with open(csv_path) as f:
        reader = csv.reader(f)
        # header = next(reader)

        charges = []
        for row in reader:
            charges.append(float(row[0]))
    
    return charges

def load_ESPs(mpac_path):
    # print('load: {}'.format(mpac_path))
    data = None
    with open(mpac_path, 'rb') as f:
        mpac_data = f.read()
        data = msgpack.unpackb(mpac_data)
        data = bridge.Utils.to_unicode_dict(data)

        grids = []
        for p in data['grids']:
            pos = bridge.Position(p[0], p[1], p[2])
            pos *= (1.0 / 0.5291772108) # Angstroam to a.u.
            grids.append(pos)
            
        ESPs = []
        for v in data['ESP']:
            ESPs.append(v)

        assert(len(grids) == len(ESPs))
    return (grids, ESPs)


def calc_RRMS(grids, ESPs, atoms):
    sum_delta2 = 0.0
    sum_v2 = 0.0
    num_of_grids = len(grids)
    # print('# of grids: {}'.format(num_of_grids))
    assert(len(ESPs) == num_of_grids)
    for i in range(num_of_grids):
        estimate_esp = calc_esp(grids[i], atoms)
        exact_esp = ESPs[i]
        delta = estimate_esp - exact_esp
        delta2 = delta * delta
        sum_delta2 += delta2
        sum_v2 += exact_esp * exact_esp
        # print("grid={} est={: 8.3e} exact={: 8.3e} >> delta2={: 8.3e}".format(grids[i], estimate_esp, exact_esp, delta2))
    rrms2 = sum_delta2 / sum_v2
    rrms = math.sqrt(rrms2)
    # print('delta2={} sum_v2={} RRMS2={} RRMS={}'.format(sum_delta2, sum_v2, rrms2, rrms))

    return rrms


def calc_esp(pos, atoms):
    esp = 0.0
    num_of_atoms = len(atoms)
    for i in range(num_of_atoms):
        d = pos.distance_from(atoms[i].xyz)
        esp += atoms[i].charge / d;
        
    return esp


def main():
    # parse args
    parser = argparse.ArgumentParser(description='calc rrms')

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

    parser.add_argument("-e", "--espdat",
                        nargs=1,
                        action='store',
                        default=['grid-esp.mpac'],
                        help="ESP values on grids by msgpack format")

    parser.add_argument("-c", "--charge_csv",
                        nargs=1,
                        action='store',
                        default=['charges.csv'],
                        help="charges by CSV format")
    
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    args = parser.parse_args()
        
    # setting
    verbose = args.verbose
    esp_data_path = args.espdat[0]
    charge_csv_path = args.charge_csv[0]

    # load grids & ESPs
    grids, ESPs = load_ESPs(esp_data_path)
    
    # setup atom list
    atomlist = []
    if args.db:
        entry = pdf.PdfArchive(args.db)
    elif args.param:
        pdfparam = pdf.load_pdfparam(args.param)
        atomgroup = pdfparam.molecule.get_atomlist()
        for k, atom in atomgroup.atoms():
            if atom.symbol != 'X':
                atomlist.append(atom)

    # charges
    charges = load_charges(charge_csv_path)
    assert(len(charges) >= len(atomlist))
    for i in range(len(atomlist)):
        atomlist[i].charge = charges[i]
                
    # calc rrms
    rrms = calc_RRMS(grids, ESPs, atomlist)
    print("RRMS={:.3f}".format(rrms))


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import csv
try:
    import msgpack
except:
    import msgpack_pure as msgpack
import math
import copy

import proteindf_bridge as bridge
import proteindf_tools as pdf

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
    parser = argparse.ArgumentParser(description='set charges to pdfparam')

    parser.add_argument('-i', '--input_path',
                        nargs=1,
                        action='store',
                        default=['pdfparam.mpac'],
                        help='ProteinDF parameter file')

    parser.add_argument("-o", "--output_path",
                        nargs=1,
                        action='store',
                        default=['new_pdfparam.mpac'],
                        help="output pdfparam")

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
    input_path = args.input_path[0]
    output_path = args.output_path[0]
    charge_csv_path = args.charge_csv[0]

    if verbose:
        sys.stderr.write("input: {}\n".format(input_path))
        sys.stderr.write("output: {}\n".format(output_path))
        sys.stderr.write("charges(CSV): {}\n".format(charge_csv_path))


    # charges
    charges = load_charges(charge_csv_path)

    # atom
    pdfparam = pdf.load_pdfparam(input_path)
    atoms = pdfparam.molecule.get_atom_list()
    assert(len(charges) >= len(atoms))
    for atom in atoms:
        if atom.symbol != 'X':
            atom.charge = charges[count]
    assert(len(charges) >= count)

    # output
    if verbose:
        sys.stderr.write("save param: {}\n".format(output_path))
    pdf.save_pdfparam(pdfparam, output_path)


if __name__ == '__main__':
    main()

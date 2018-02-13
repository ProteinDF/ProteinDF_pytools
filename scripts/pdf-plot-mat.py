#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from types import *
import argparse
import logging
try:
    import msgpack
except:
    try:
        import umsgpack as msgpack
    except:
        import msgpack_pure as msgpack

import proteindf_bridge as bridge
import proteindf_tools as pdf

def main():
    # parse args
    parser = argparse.ArgumentParser(description='plot matrix')
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
    parser.add_argument('matrix',
                        nargs=1,
                        help='matrix path')
    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=['mat.png'],
                        help='output graph path')
    parser.add_argument("-2", "--diverging",
                        action="store_true",
                        default=False)
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
    mat_path = args.matrix[0]
    output_path = args.output[0]
    is_diverging = args.diverging

    # setup orbinfo
    #if verbose:
    #    sys.stderr.write("prepare orbital information\n")
    #orb_info = None
    #molecule = None
    #if args.db:
    #    entry = pdf.PdfArchive(args.db)
    #    orb_info = pdf.OrbInfo(entry)
    #    molecule = entry.get_molecule()
    #elif args.param:
    #    pdfparam = pdf.load_pdfparam(args.param)
    #    orb_info = pdf.OrbInfo(pdfparam)
    #    molecule = pdfparam.molecule

    # molecule
    #if verbose:
    #    sys.stderr.write("prepare molecule information\n")
    #distance_mat = distance_matrix(molecule)

    # setup matrix
    mat = get_matrix(mat_path, verbose)

    # matrix elements
    mat_plot = pdf.DfMatrixGraph(mat, is_diverging)
    mat_plot.save(output_path)


    # distance v.s. matrix elements
    #dist_vs_matvalue(orb_info, distance_mat, mat, 'mat.dat')
    #distgraph = pdf.DfDistanceVsElementGraph2()
    #distgraph.load('mat.dat')
    #distgraph.save("dist.png")


def get_matrix(path, verbose = False):
    if verbose:
        sys.stderr.write("load matrix: {}\n".format(path))
    mat = None

    if pdf.SymmetricMatrix.is_loadable(path):
        mat = pdf.SymmetricMatrix()
        mat.load(path)
        mat = mat.get_general_matrix()
    elif pdf.Matrix.is_loadable(path):
        mat = pdf.Matrix()
        mat.load(path)
    else:
        print("matrix file cannnot load.")
        raise

    return mat


def distance_matrix(atomgroup, verbose = False):
    '''
    create distance matrix

    support flat-atomgroup only
    '''
    atomlist = atomgroup.get_atom_list()
    num_of_atoms = len(atomlist)
    if verbose:
        sys.stderr.write("num of atoms: {}n".format(num_of_atoms))
    mat = bridge.SymmetricMatrix(num_of_atoms)

    for i in range(num_of_atoms):
        for j in range(i):
            distance = atomlist[i].xyz.distance_from(atomlist[j].xyz)
            mat.set(i, j, distance)

    return mat


def dist_vs_matvalue(orb_info, distance_mat, mat, path, verbose = False):
    if verbose:
        sys.stderr.write("matrix size: {} x {}\n".format(mat.rows, mat.cols))

    f = open(path, 'w')
    for i in range(mat.rows):
        atom_i = orb_info.get_atom_id(i)
        for j in range(i):
            atom_j = orb_info.get_atom_id(j)

            d = distance_mat.get(atom_i, atom_j)
            v = mat.get(i, j)
            data_str = "{}, {}\n".format(d, v)
            f.write(data_str)
    f.close()


if __name__ == '__main__':
    main()

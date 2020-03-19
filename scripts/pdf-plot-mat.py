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
    parser.add_argument('matrix',
                        nargs=1,
                        help='matrix path')
    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=['mat.png'],
                        help='output graph path')
    parser.add_argument("--vmax",
                        type=float)
    parser.add_argument("--vmin",
                        type=float)
    parser.add_argument('--title',
                        nargs=1,
                        default=['matrix value'],
                        help='graph title')
    parser.add_argument('--cmap',
                        nargs=1,
                        default=['bwr'],
                        help='color map name')
    parser.add_argument("--write-values",
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
    if verbose:
        print("matrix path: {}".format(mat_path))

    output_path = args.output[0]
    if verbose:
        print("output path: {}".format(output_path))

    vmax = None
    if args.vmax != None:
        vmax = args.vmax
        if verbose:
            print("vmax: {}".format(vmax))
    vmin = None
    if args.vmin != None:
        vmin = args.vmin
        if verbose:
            print("vmin: {}".format(vmin))

    title = args.title[0]

    cmap = args.cmap[0]
    if verbose:
        print("cmap: {}".format(cmap))

    write_values = args.write_values
    if verbose:
        print("write values: {}".format(write_values))

    # setup matrix
    mat = get_matrix(mat_path, verbose)
    rows = mat.rows
    cols = mat.cols
    xticks = [x for x in range(0, rows)]
    yticks = [y for y in range(0, cols)]
    xticklabels = ["{}".format(x) for x in range(1, rows + 1)]
    yticklabels = ["{}".format(y) for y in range(1, cols + 1)]

    # matrix elements
    mat_plot = pdf.DfMatrixGraph(mat, figsize=(20, 20))
    if vmax:
        mat_plot.vmax = vmax
    if vmin:
        mat_plot.vmin = vmin
    mat_plot.xticks = xticks
    mat_plot.yticks = yticks
    mat_plot.xticklabels = xticklabels
    mat_plot.yticklabels = yticklabels
    mat_plot.title = title
    mat_plot.cmap = cmap
    mat_plot.should_write_values = write_values

    mat_plot.save(output_path)


def get_matrix(path, verbose=False):
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


def distance_matrix(atomgroup, verbose=False):
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


def dist_vs_matvalue(orb_info, distance_mat, mat, path, verbose=False):
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

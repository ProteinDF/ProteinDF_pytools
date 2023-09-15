#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from types import *
import argparse

import proteindf_bridge as bridge
import proteindf_tools as pdf

import logging


def main():
    # parse args
    parser = argparse.ArgumentParser(description="plot vector")
    parser.add_argument("vector", nargs=1, help="vector path")
    parser.add_argument("-o", "--output", nargs=1, default=["vct.png"], help="output graph path")
    parser.add_argument("--abs", action="store_true", default=False, help="plot absolute values")
    parser.add_argument("--draw-grid", action="store_true", default=False, help="draw grid")
    parser.add_argument("--yscale", type=str)
    parser.add_argument("--ymax", type=float)
    parser.add_argument("--ymin", type=float)
    # parser.add_argument("--title", nargs=1, default=["matrix value"], help="graph title")
    # parser.add_argument("--no-x-tick-labels", action="store_true", default=False, help="draw x-tick labels")
    # parser.add_argument("--no-y-tick-labels", action="store_true", default=False, help="draw y-tick labels")
    # parser.add_argument("--cmap", nargs=1, default=["bwr"], help="color map name")
    # parser.add_argument("--write-values", action="store_true", default=False)
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("-D", "--debug", action="store_true", default=False)
    args = parser.parse_args()

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    vct_path = args.vector[0]
    if verbose:
        print("vector path: {}".format(vct_path))

    output_path = args.output[0]
    if verbose:
        print("output path: {}".format(output_path))

    plot_abs_values = args.abs

    is_draw_grid = args.draw_grid

    yscale = None
    if args.yscale != None:
        yscale = args.yscale
        if verbose:
            print("yscale: {}".format(yscale))
    ymax = None
    if args.ymax != None:
        ymax = args.ymax
        if verbose:
            print("ymax: {}".format(ymax))
    ymin = None
    if args.ymin != None:
        ymin = args.ymin
        if verbose:
            print("ymin: {}".format(ymin))

    # title = args.title[0]
    # x_tick_labels = not args.no_x_tick_labels
    # y_tick_labels = not args.no_y_tick_labels
    # if verbose:
    #     print("draw x tick labels: {}".format(x_tick_labels))
    #     print("draw y tick labels: {}".format(y_tick_labels))

    # cmap = args.cmap[0]
    # if verbose:
    #     print("cmap: {}".format(cmap))

    # write_values = args.write_values
    # if verbose:
    #     print("write values: {}".format(write_values))

    # setup vector
    vct = get_vector(vct_path, verbose)
    if plot_abs_values:
        vct = vct.abs()

    # vector elements
    vct_plot = pdf.DfVectorGraph(vct, figsize=(20, 20))

    vct_plot.is_draw_grid = is_draw_grid    
    vct_plot.yscale = yscale
    if ymax:
        vct_plot.ymax = ymax
    if ymin:
        vct_plot.ymin = ymin

    # mat_plot.xticks = []
    # mat_plot.xticklabels = []
    # if x_tick_labels:
    #     xticks = [x for x in range(0, rows)]
    #     xticklabels = ["{}".format(x) for x in range(1, rows + 1)]
    #     mat_plot.xticks = xticks
    #     mat_plot.xticklabels = xticklabels

    # mat_plot.yticks = []
    # mat_plot.yticklabels = []
    # if y_tick_labels:
    #     yticks = [y for y in range(0, cols)]
    #     yticklabels = ["{}".format(y) for y in range(1, cols + 1)]
    #     mat_plot.yticks = yticks
    #     mat_plot.yticklabels = yticklabels
    # mat_plot.title = title
    # mat_plot.cmap = cmap
    # mat_plot.should_write_values = write_values

    vct_plot.save(output_path)


def get_vector(path, verbose=False):
    if verbose:
        sys.stderr.write("load matrix: {}\n".format(path))
    vct = None

    if pdf.Vector.is_loadable(path):
        vct = pdf.Vector()
        vct.load(path)
    else:
        print("vector file cannnot load.")
        raise

    return vct


# def distance_matrix(atomgroup, verbose=False):
#     """
#     create distance matrix

#     support flat-atomgroup only
#     """
#     atomlist = atomgroup.get_atom_list()
#     num_of_atoms = len(atomlist)
#     if verbose:
#         sys.stderr.write("num of atoms: {}n".format(num_of_atoms))
#     mat = bridge.SymmetricMatrix(num_of_atoms)

#     for i in range(num_of_atoms):
#         for j in range(i):
#             distance = atomlist[i].xyz.distance_from(atomlist[j].xyz)
#             mat.set(i, j, distance)

#     return mat


# def dist_vs_matvalue(orb_info, distance_mat, mat, path, verbose=False):
#     if verbose:
#         sys.stderr.write("matrix size: {} x {}\n".format(mat.rows, mat.cols))

#     f = open(path, "w")
#     for i in range(mat.rows):
#         atom_i = orb_info.get_atom_id(i)
#         for j in range(i):
#             atom_j = orb_info.get_atom_id(j)

#             d = distance_mat.get(atom_i, atom_j)
#             v = mat.get(i, j)
#             data_str = "{}, {}\n".format(d, v)
#             f.write(data_str)
#     f.close()


if __name__ == "__main__":
    main()

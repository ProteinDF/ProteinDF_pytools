#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import csv
import logging
import math
import os
import sys
from types import *

import numpy
import proteindf_bridge as bridge
import proteindf_tools as pdf
from matplotlib import use


def output_csv(csv_data, path):
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(csv_data)


def main():
    # parse args
    parser = argparse.ArgumentParser(description='plot G-values')
    parser.add_argument('-d', '--db',
                        nargs=1,
                        default=['pdfresults.h5'],
                        help='ProteinDF results file')

    parser.add_argument('-i', '--iteration',
                        nargs=1,
                        default=[0],
                        help='SCF iteration')

    parser.add_argument('--occ',
                        action="store_true",
                        default=False,
                        help='occ only')

    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=['g-values.png'],
                        help='output graph path')
    parser.add_argument('--csv',
                        nargs=1,
                        default=['dos.csv'],
                        help='output csv path')

    parser.add_argument("--frontier",
                        action="store_true",
                        default=False)

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    parser.add_argument('-D', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    # print(args)

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    pdfresults_db = args.db[0]

    iteration = int(args.iteration[0])
    output_path = args.output[0]
    csv_path = args.csv[0]
    is_frontier = args.frontier

    AU2EV = 27.2116

    # load DB
    pdfparam = pdf.PdfParam_H5()
    pdfparam.open(pdfresults_db)

    # eLevels
    method = pdfparam.method
    run_type = "rks"
    HOMO_level = pdfparam.get_HOMO_level(run_type)

    if iteration == 0:
        iteration = pdfparam.iterations
    print("iteration: {}".format(iteration))
    eigvals = pdfparam.get_energy_level(method, iteration)
    eigvals *= AU2EV

    # make CSC and g-values
    CSC_mat_path = "CSC.mat"
    g_vtr_path = "g-values.vtr"
    args1 = [
        "component",
        "-v",
        "-i", iteration,
        "-S", CSC_mat_path,
        "-g", g_vtr_path
    ]
    # print(args)
    pdf.run_pdf(args1)

    vtr = pdf.Vector()
    vtr.load(g_vtr_path)

    num_of_MOs = len(vtr)
    print("HOMO: {}".format(HOMO_level))

    data_x_occ = [ eigvals[i] for i in range(HOMO_level +1)]
    data_y_occ = [ vtr.get(i) for i in range(HOMO_level +1)]
    data_x_vir = [ eigvals[i] for i in range(HOMO_level +1, len(vtr))]
    data_y_vir = [ vtr.get(i) for i in range(HOMO_level +1, len(vtr))]

    # plot data
    graph = pdf.DfGvalsGraph()
    if is_frontier:
        graph.xmin = -20.0
        graph.xmax = 5.0

    graph.set_data(data_x_occ, data_y_occ, data_x_vir, data_y_vir)
    print(output_path)
    graph.save(output_path)


if __name__ == '__main__':
    main()

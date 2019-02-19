#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import logging

import proteindf_bridge as bridge
import proteindf_tools as pdf


def make_pop_data(pop_vtr, output_path):
    num_of_items = len(pop_vtr)
    with open(output_path, "w") as f:
        for index in range(num_of_items):
            charge = pop_vtr[index]
            f.write("{}, {}\n".format(index +1, charge))


def make_pop_data_diff(pop_vtr1, pop_vtr2, output_path):
    num_of_items = len(pop_vtr1)
    assert(num_of_items == len(pop_vtr2))

    with open(output_path, "w") as f:
        for index in range(num_of_items):
            charge = pop_vtr1[index] - pop_vtr2[index]
            f.write("{}, {}\n".format(index +1, charge))


def main():
    # parse args
    parser = argparse.ArgumentParser(description='plot population')
    group = parser.add_mutually_exclusive_group(required=True)
    #group.add_argument('-e', '--elevel_vector',
    #                   nargs=1,
    #                   help='ProteinDF energy level vector file')
    group.add_argument('--h5',
                       nargs='?',
                       action='store',
                       const='pdfresults.h5',
                       help='ProteinDF results file (HDF5 format)')

    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=['pop.png'],
                        help='output graph path')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    args = parser.parse_args()


    # setting
    verbose = args.verbose
    pdfparam_h5_path = args.h5
    output_path = args.output[0]

    #
    if verbose:
        print("H5 path: {}".format(pdfparam_h5_path))
        print("output path: {}".format(output_path))

    #
    pdfparam = pdf.PdfParam_H5()
    pdfparam.open(pdfparam_h5_path)

    run_type = "rks"
    itr = pdfparam.iterations

    pop_vtr_path1 = pdfparam.get_pop_mulliken_path(run_type, itr)
    if os.path.exists(pop_vtr_path1):
        pop_vtr1 = pdf.Vector()
        pop_vtr1.load(pop_vtr_path1)

    pop_vtr_path2 = pdfparam.get_pop_mulliken_path(run_type, 1)
    if os.path.exists(pop_vtr_path2):
        pop_vtr2 = pdf.Vector()
        pop_vtr2.load(pop_vtr_path2)


    data_path = "pop.mulliken.txt"
    # make_pop_data(pop_vtr, data_path)
    make_pop_data_diff(pop_vtr1, pop_vtr2, data_path)

    # draw
    pop_graph = pdf.DfPopulationGraph()
    pop_graph.load_data(data_path)
    pop_graph.save(output_path)


if __name__ == '__main__':
    main()

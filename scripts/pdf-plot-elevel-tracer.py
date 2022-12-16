#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import logging

import proteindf_tools as pdf

AU2EV = 27.2116


def find_pdfparam_path(pdf_path, pdfparam_filename="pdfparam.mpac"):
    pdfparam_path = None
    path = os.path.join(pdf_path, pdfparam_filename)
    if os.path.isfile(path):
        pdfparam_path = path

    return pdfparam_path


def find_db(pdf_path, h5_filename="pdfresults.h5"):
    h5_path = None
    path = os.path.join(pdf_path, h5_filename)
    if os.path.isfile(path):
        h5_path = path

    return h5_path


def get_entry(pdfparam_path, h5_path, verbose=False):
    entry = None
    if h5_path:
        if verbose:
            print("load: {}".format(h5_path))
        entry = pdf.PdfParam_H5(h5_path)
    elif pdfparam_path:
        if verbose:
            print("load: {}".format(pdfparam_path))
        entry = pdf.load_pdfparam(pdfparam_path)
    else:
        print("cannot load parameter files")

    return entry


def make_corresponding_orbital_matrix(path1, path2, pdfparam_path1, pdfparam_path2, CSC_mat_path):
    pdfparam1 = pdf.load_pdfparam(pdfparam_path1)
    cmat_path1 = pdfparam1.get_c_mat_path()

    pdfparam2 = pdf.load_pdfparam(pdfparam_path2)
    cmat_path2 = pdfparam2.get_c_mat_path()

    args = [
        "info-corresponding-orbital",
        "-v",
        "-o",
        CSC_mat_path,
        pdfparam_path1,
        pdfparam_path2,
        os.path.join(path1, cmat_path1),
        os.path.join(path2, cmat_path2),
    ]
    print(args)

    pdf.run_pdf(args)


def make_pair(CSC_mat, pickup_ratio=0.1):
    pickup_ratio = float(pickup_ratio)
    print("pickup ratio: {}".format(pickup_ratio))
    pairs = []

    rows = CSC_mat.rows
    cols = CSC_mat.cols
    for r in range(rows):
        # print("make pair for [{}]".format(r))
        col_vec = CSC_mat.get_row_vector(r)
        norm_col_vec = pdf.Vector(cols)
        norm = 0.0
        for c in range(cols):
            v = abs(col_vec[c])
            norm_col_vec[c] = v
            norm += v
        norm_col_vec *= 1.0 / norm
        col_index = norm_col_vec.argsort().flip()
        assert len(col_index) == cols

        # pickup by ratio
        rank = 0
        ratio = 0.0
        while ratio < 0.8:
            index = col_index[rank]
            value = norm_col_vec[index]

            if value > pickup_ratio:
                pair = [r, index, value]
                pairs.append(pair)
                print("pickup[{}]: {} -> {:5.2f}%".format(rank, index, value * 100.0))

            rank += 1
            ratio += value

    return pairs


def main():
    # parse args
    parser = argparse.ArgumentParser(description="plot energy levels")

    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("-D", "--debug", action="store_true", default=False)

    parser.add_argument("-c", "--corresponding-orbital-matrix-path", default="CSC.mat")

    parser.add_argument("--pickup-ratio", default=0.1, type=float)

    parser.add_argument("pdf_path1", help="ProteinDF path1")
    parser.add_argument("pdf_path2", help="ProteinDF path2")
    args = parser.parse_args()
    # print(args)

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    pickup_ratio = args.pickup_ratio

    path1 = args.pdf_path1
    path2 = args.pdf_path2
    # print(path1, path2)
    pdfparam_path1 = find_pdfparam_path(path1)
    pdfparam_path2 = find_pdfparam_path(path2)

    # entry1 = get_entry(pdfparam_path1, None, verbose)
    # entry2 = get_entry(pdfparam_path2, None, verbose)

    # itr1 = entry1.iterations
    # itr2 = entry2.iterations
    # print(itr1, itr2)

    corresponding_orbital_matrix_path = args.corresponding_orbital_matrix_path
    # print(CSC_matrix_path)

    make_corresponding_orbital_matrix(path1, path2, pdfparam_path1, pdfparam_path2, corresponding_orbital_matrix_path)

    db_path1 = find_db(path1)
    db_path2 = find_db(path2)
    output_path = "elevel-tracer.png"

    # load DBs
    db1 = pdf.PdfParam_H5(db_path1)
    db2 = pdf.PdfParam_H5(db_path2)

    runtype = "rks"
    itr1 = db1.iterations
    if verbose:
        print("iteration1: {}".format(itr1))
    elevel1 = db1.get_energy_level(runtype, itr1)
    elevel1 *= AU2EV
    if verbose:
        print("size of eigvals1: {}".format(len(elevel1)))

    itr2 = db2.iterations
    if verbose:
        print("iteration2: {}".format(itr2))
    elevel2 = db2.get_energy_level(runtype, itr2)
    elevel2 *= AU2EV
    if verbose:
        print("size of eigvals2: {}".format(len(elevel2)))

    # corresponding orbital
    CSC_mat = pdf.Matrix()
    CSC_mat.load(corresponding_orbital_matrix_path)
    print("CSC matrix: {} x {}".format(CSC_mat.rows, CSC_mat.cols))

    # make pair
    print("> make pair")
    pairs = make_pair(CSC_mat, pickup_ratio)

    # graph
    print("> make graph")
    graph = pdf.DfEnergyLevelTraceGraph()
    graph.ymin = -30
    graph.ymax = 10
    graph.set_data(elevel1.to_list(), elevel2.to_list(), pairs)
    graph.save(output_path)


if __name__ == "__main__":
    main()

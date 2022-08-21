#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging

import proteindf_tools as pdf

AU2EV = 27.2116


def make_pair(CO_mat):
    pairs = []

    rows = CO_mat.rows
    cols = CO_mat.cols
    for r in range(rows):
        # print("make pair for [{}]".format(r))
        col_vec = CO_mat.get_row_vector(r)
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

            if value > 0.3:
                pair = [r, index, value]
                pairs.append(pair)
                # print("pickup[{}]: {} -> {:5.2f}%".format(rank, index, value * 100.0))

            rank += 1
            ratio += value

    return pairs


def main():
    # parse args
    parser = argparse.ArgumentParser(description="plot energy levels")

    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("-D", "--debug", action="store_true", default=False)

    parser.add_argument("-c", "--corresponding-orbital-matrix-path")

    parser.add_argument("db_path1", help="pdfresults.h5 path(left)")
    parser.add_argument("db_path2", help="pdfresults.h5 path(right)")
    args = parser.parse_args()
    # print(args)

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    corresponding_orbital_matrix_path = ""
    # print(args.corresponding_orbital_matrix_path)
    if args.corresponding_orbital_matrix_path:
        corresponding_orbital_matrix_path = args.corresponding_orbital_matrix_path

    db_path1 = args.db_path1
    db_path2 = args.db_path2
    output_path = "elevel-tracer.png"

    if verbose:
        print("DB1: {}".format(db_path1))
        print("DB2: {}".format(db_path2))

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
    CO_mat = pdf.Matrix()
    CO_mat.load(corresponding_orbital_matrix_path)
    print("CO matrix: {} x {}".format(CO_mat.rows, CO_mat.cols))

    # make pair
    print("> make pair")
    pairs = make_pair(CO_mat)

    # graph
    print("> make graph")
    graph = pdf.DfEnergyLevelTraceGraph()
    graph.set_data(elevel1.to_list(), elevel2.to_list(), pairs)
    graph.save(output_path)


if __name__ == "__main__":
    main()

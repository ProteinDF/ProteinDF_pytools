#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse

import proteindf_bridge as bridge
import proteindf_tools as pdf

import logging
logger = logging.getLogger(__name__)


def main():
    # parse args
    parser = argparse.ArgumentParser(description='show matrix info')
    parser.add_argument("path",
                        nargs=1,
                        help='matrix path')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    parser.add_argument('-D', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    # setting
    mat_path = args.path[0]
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    if verbose:
        print("loading: {}".format(mat_path))

    stdmat = pdf.Matrix()
    symmat = pdf.SymmetricMatrix()

    if stdmat.is_loadable(mat_path):
        stdmat.load(mat_path)
        print("standard dens matrix")
        print("size: {row} x {col}".format(row=stdmat.rows,
                                           col=stdmat.cols))
        # print(stdmat)
    elif symmat.is_loadable(mat_path):
        symmat.load(mat_path)
        print("symmetric dens matrix")
        print("size: {row} x {col}".format(row=symmat.rows,
                                           col=symmat.cols))
        # print(symmat)
    else:
        print("cannot load file.")


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2014 The ProteinDF development team.
# see also AUTHORS and README if provided.
#
# This file is a part of the ProteinDF software package.
#
# The ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
from types import *
import argparse
import logging
try:
    import msgpack
except:
    import msgpack_pure as msgpack

import proteindf_bridge as bridge
import proteindf_tools as pdf

def main():
    # parse args
    parser = argparse.ArgumentParser(description='plot DOS)')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-e', '--elevel_vector',
                       nargs=1,
                       help='ProteinDF energy level vector file')
    group.add_argument('-a', '--alpha',
                       nargs=1,
                       type=float,
                       default=[1.0],
                       help='alpha parameter')
    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=['dos.png'],
                        help='output graph path')
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
    output_path = args.output[0]
    HOMO_level = int(args.level[0])

    eigval_path = args.elevel_vector[0]
    eigvals = pdf.Vector()
    alpha = args.alpha[0]
    if verbose:
        print('load: %s' % (eigval_path))
        print('output: %s' % (output_path))
        print('alpha: %d' % (alpha))

    x_min = -20.0
    x_max = 10.0
    x_pitch = 0.01
    x = [x for x in range(x_min, x_max, x_pitch)]
    y = [dos(alpha, vct, i) for i in x]

    eiggvals.load(eigval_path)
    elevel2dat(itr, eigvals, data_path)
    plot_elevel_single(data_path, itr, HOMO_level, output_path)

def elevel2dat(itr, eigvals, data_path):
    """
    エネルギー準位ベクトルファイルを
    グラフ描画用データに変換する

    itr: iteraton(int)
    eigvals: energy level(bridge.Vector object)
    data_path: output path of data
    """
    assert(isinstance(itr, IntType))
    assert(isinstance(eigvals, bridge.Vector))
    assert(isinstance(data_path, StringType))

    dat = open(data_path, 'w')
    for level, e in enumerate(eigvals):
        e *= 27.2116
        dat.write('%d, %d, % 16.10f\n' % (itr, level, e))
    dat.close()

def func_gauss(alpha, x):
    alpha = float(alpha)
    x = float(x)
    answer = exp(- alpha * x * x)

    return answer

def func_dos(alpha, vct, x):
    assert (isinstance(vct, bridge.Vector))

    answer = 0.0
    vct_size = vct.size()
    for i in range(vct_size):
        v = vct.get(i)
        answer += func_gauss(alpha, v)
    
    return answer

if __name__ == '__main__':
    main()

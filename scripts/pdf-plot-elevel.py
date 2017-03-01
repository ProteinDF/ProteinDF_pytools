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
    
import pdfbridge
import pdfpytools as pdf

def main():
    # parse args
    parser = argparse.ArgumentParser(description='plot energy level(single; vertical)')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-e', '--elevel_vector',
                       nargs=1,
                       help='ProteinDF energy level vector file')
    group.add_argument('-d', '--db',
                       nargs='?',
                       action='store',
                       const='pdfresults.db',
                       help='ProteinDF results file')
    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=['elevel.png'],
                        help='output graph path')
    parser.add_argument('-l', '--level',
                        nargs=1,
                        default=['0'],
                        help='HOMO level')
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

    if args.elevel_vector:
        eigval_path = args.elevel_vector[0]

        itr = 1
        data_path = 'tmp.dat'

        eigvals = pdf.Vector()
        if verbose:
            print('load: %s' % (eigval_path))
            print('output: %s' % (output_path))
            print('HOMO: %d' % (HOMO_level))
        eigvals.load(eigval_path)
        elevel2dat(itr, eigvals, data_path)
        plot_elevel_single(data_path, itr, HOMO_level, output_path)
    else:
        entry = pdf.PdfArchive(args.db)

        itr = entry.iterations
        HOMO_level = entry.get_HOMO_level('ALPHA') # TODO
        data_path = 'tmp.dat'

        eigvals = entry.get_energylevel('ALPHA', itr) # TODO
        if verbose:
            print('load: %s' % (eigval_path))
            print('output: %s' % (output_path))
            print('HOMO: %d' % (HOMO_level))
        eigvals.load(eigval_path)
        elevel2dat(itr, eigvals, data_path)
        plot_elevel_single(data_path, itr, HOMO_level, output_path)
        
    
def elevel2dat(itr, eigvals, data_path):
    """
    エネルギー準位ベクトルファイルを
    グラフ描画用データに変換する

    itr: iteraton(int)
    eigvals: energy level(pdfbridge.Vector object)
    data_path: output path of data
    """
    assert(isinstance(itr, IntType))
    assert(isinstance(eigvals, pdfbridge.Vector))
    assert(isinstance(data_path, StringType))
    
    dat = open(data_path, 'w')
    for level, e in enumerate(eigvals):
        e *= 27.2116
        dat.write('%d, %d, % 16.10f\n' % (itr, level, e))
    dat.close()
    
def plot_elevel_single(data_path, itr, HOMO_level, output_path):
    """
    横向きのエネルギー準位グラフを描画する

    data_path: グラフ描画用データファイルのパス
    HOMO: 最高占有軌道(0からスタート)
    itr: iteration
    output_path: 
    """
    graphV = pdf.DfEnergyLevelHistoryGraphV()
    graphV.set_HOMO_level(HOMO_level) # option base 0
    graphV.load_data(data_path)
    graphV.select_iterations([itr])
    graphV.is_draw_grid = False
    graphV.yticks = []
    graphV.ylabel = ''
    graphV.ymin = 0.4
    graphV.ymax = 1.6
    graphV.aspect = 3.0
    
    graphV.save(output_path)

if __name__ == '__main__':
    main()
    

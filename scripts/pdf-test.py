#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2002-2014 The ProteinDF project
# see also AUTHORS and README.
# 
# This file is part of ProteinDF.
# 
# ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

"""
compare ProteinDF results 
"""

import sys
import os.path
import argparse
import logging
try:
    import msgpack
except:
    import msgpack_pure as msgpack

import pdf

def main():
    parser = argparse.ArgumentParser(description='compare ProteinDF results')
    parser.add_argument('FILE1',
                        nargs=1,
                        help='ProteinDF parameter file1')
    parser.add_argument('FILE2',
                        nargs=1,
                        help='ProteinDF parameter file2')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    verbose = args.verbose
    logging.basicConfig()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)

    path1 = args.FILE1[0]
    path2 = args.FILE2[0]

    if not os.path.isfile(path1):
        sys.exit('file not found: %s' % (path1))
    if not os.path.isfile(path2):
        sys.exit('file not found: %s' % (path2))

    data1 = pdf.PdfArchive(path1)
    data2 = pdf.PdfArchive(path2)
    
    if data1 == data2:
        logging.debug('ProteinDF results are OK.')
        sys.exit(0)
    else:
        logging.error('ProteinDF results are not consistent.')
        sys.exit(1)

        
if __name__ == '__main__':
    main()


#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
compare ProteinDF results 
"""

import sys
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
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    path1 = args.FILE1[0]
    path2 = args.FILE2[0]

    data1 = pdf.PdfArchive(path1)
    data2 = pdf.PdfArchive(path2)

    if data1 == data2:
        logging.debug('ProteinDF results are OK.')
        sys.exit(0)
    else:
        logging.warning('ProteinDF results are not consistent.')
        sys.exit(1)
    
if __name__ == '__main__':
    main()


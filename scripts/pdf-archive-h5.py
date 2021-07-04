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

"""
archive ProteinDF results.
"""

import os
import sys
import argparse
import traceback
#import pprint

import proteindf_bridge as bridge
import proteindf_tools as pdf

import logging
import logging.config


def main():
    parser = argparse.ArgumentParser(description='archive ProteinDF results')
    parser.add_argument('-p', '--pdfparam',
                        nargs=1,
                        dest='pdfparam_path',
                        default='pdfparam.mpac',
                        help='ProteinDF parameter file')
    parser.add_argument('-o', '--output',
                        dest='output',
                        nargs='?',
                        default='pdfresults.h5',
                        help='ProteinDF results file')

    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    # setting from command-line
    pdfparam_path = args.pdfparam_path
    output = args.output
    verbose = args.verbose
    #is_little_endian = True

    try:
        # read ProteinDF parameters
        pdfparam = pdf.load_pdfparam(pdfparam_path)

        # setup DB
        entry = pdf.PdfParam_H5(pdfparam)
        entry.save_standard(output)

    except:
        print('-'*60)
        traceback.print_exc(file=sys.stdout)
        # print(traceback.format_exc())
        print('-'*60)


if __name__ == '__main__':
    logconfig_file = 'logconfig.ini'
    if os.path.exists(logconfig_file):
        logging.config.fileConfig(logconfig_file)
    main()

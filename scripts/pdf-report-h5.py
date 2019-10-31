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
import argparse
import logging
import logging.config
logger = logging.getLogger(__name__)

try:
    import msgpack
except:
    import msgpack_pure as msgpack

import proteindf_tools as pdf

def main():
    # parse args
    parser = argparse.ArgumentParser(description='make ProteiDF report')
    parser.add_argument('pdfresults_db',
                        nargs='?',
                        default='pdfresults.h5',
                        help='ProteinDF results file')
    parser.add_argument('output_dir',
                        nargs='?',
                        default='report',
                        help='output directory')

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        default=False)
    parser.add_argument("-p", "--profile",
                        nargs='?',
                        const='pdf-report.stat',
                        help="do profiling")

    args = parser.parse_args()

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    do_profile = args.profile

    output_dir=args.output_dir
    pdfparam = pdf.PdfParam_H5()
    pdfparam.open(args.pdfresults_db)

    # run
    report = pdf.PdfReport(pdfparam=pdfparam, workdir=output_dir)
    if do_profile != None:
       import cProfile
       cProfile.runctx("report()",
                       globals(), locals(),
                       report.report())
    else:
       report.report()


if __name__ == '__main__':
    if os.path.exists("config.ini"):
        logging.config.fileConfig("config.ini",
                                  disable_existing_loggers=False)
    main()

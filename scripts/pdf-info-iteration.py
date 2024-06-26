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
output geometory in XYZ format
"""

import os
import sys
import argparse

import proteindf_bridge as bridge
import proteindf_tools as pdf


def main():
    # parse args
    parser = argparse.ArgumentParser(description="output geometory in XYZ format")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("-d", "--db", nargs="?", action="store", const="pdfresults.h5", help="ProteinDF results file")
    group.add_argument(
        "-p", "--param", nargs="?", action="store", const="pdfparam.mpac", help="ProteinDF parameter file"
    )
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("-D", "--debug", action="store_true", default=False)
    args = parser.parse_args()

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    entry = None
    if args.db:
        entry = pdf.PdfParam_H5()
        entry.open(args.db)
    elif args.param:
        entry = pdf.load_pdfparam(args.param)

    iterations = entry.iterations
    print(iterations)


if __name__ == "__main__":
    main()

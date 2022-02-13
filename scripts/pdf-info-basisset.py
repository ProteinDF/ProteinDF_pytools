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
compare ProteinDF results 
"""

import sys
import os.path
import argparse

import proteindf_bridge as bridge
import proteindf_tools as pdf

import logging


def main():
    parser = argparse.ArgumentParser(
        description='show ProteinDF basis information')
    parser.add_argument('basisset_name',
                        nargs=1,
                        help='ProteinDF basisset name')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    verbose = args.verbose

    basisset_name = args.basisset_name[0]

    basis2 = pdf.Basis2()
    basisset = basis2.get_basisset(basisset_name)

    print(basisset)

    sys.exit(0)


if __name__ == '__main__':
    main()

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
import logging
try:
    import msgpack
except:
    import msgpack_pure as msgpack

import proteindf_tools as pdf


def get_OK_NG(yn):
    if yn:
        return "OK"
    else:
        return "NG"


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
    parser.add_argument('--compare-info',
                        action='store_true',
                        help='check info',
                        default=False)
    parser.add_argument('--compare-energy',
                        action='store_true',
                        help='check TE',
                        default=False)
    parser.add_argument('--compare-pop',
                        action='store_true',
                        help='check population',
                        default=False)
    parser.add_argument('--compare-grad',
                        action='store_true',
                        help='check gradient',
                        default=False)
    args = parser.parse_args()

    verbose = args.verbose
    debug = args.debug

    compare_info = args.compare_info
    compare_energy = args.compare_energy
    compare_pop = args.compare_pop
    compare_grad = args.compare_grad

    # 全てoffならenergyをon
    if (compare_info | compare_energy | compare_pop | compare_grad) == 0:
        compare_info = True
        compare_energy = True

    # setup logger
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.INFO)
    logger = logging.getLogger(__name__)

    logger.debug('check items')
    logger.debug(' info: {}'.format(compare_info))
    logger.debug(' TE: {}'.format(compare_energy))
    logger.debug(' pop: {}'.format(compare_pop))
    logger.debug(' grad: {}'.format(compare_grad))

    path1 = args.FILE1[0]
    path2 = args.FILE2[0]

    if not os.path.isfile(path1):
        sys.exit('file not found: %s' % (path1))
    if not os.path.isfile(path2):
        sys.exit('file not found: %s' % (path2))

    data1 = pdf.PdfArchive(path1)
    data2 = pdf.PdfArchive(path2)


    # check ------------------------------------------------------------
    answer = True
    if compare_info:
        judge = data1.compare_info(data2)
        logger.info('check info: {}'.format(get_OK_NG(judge)))
        answer = answer & judge

    if compare_energy:
        judge = data1.compare_energy(data2)
        logger.info('check TE: {}'.format(get_OK_NG(judge)))
        answer = answer & judge

    if compare_pop:
        judge = data1.compare_pop(data2)
        logger.info('check pop: {}'.format(get_OK_NG(judge)))
        answer = answer & judge

    if compare_grad:
        judge = data1.compare_grad(data2)
        logger.info('check grad: {}'.format(get_OK_NG(judge)))
        answer = answer & judge

    # output -----------------------------------------------------------
    if answer:
        logger.info('check OK')
        sys.exit(0)
    else:
        logger.error('ProteinDF results are not consistent.')
        sys.exit(1)

if __name__ == '__main__':
    if os.path.exists("config.ini"):
        logging.config.fileConfig("config.ini",
                                  disable_existing_loggers=False)
    main()

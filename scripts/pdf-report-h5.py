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
import yaml
import rich.logging

import proteindf_bridge as bridge
import proteindf_tools as pdf

import logging
import logging.config
logger = logging.getLogger(__name__)


def setup_logger(is_debug=False, output_file='pdf-report-h5.log'):
    handlers = []
    level = logging.INFO
    if is_debug:
        level = logging.DEBUG

    rich_handler = rich.logging.RichHandler(rich_tracebacks=True)
    rich_handler.setLevel(level)
    rich_handler.setFormatter(logging.Formatter("%(message)s"))
    handlers.append(rich_handler)

    if len(output_file) > 0:
        file_handler = logging.FileHandler(output_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(logging.Formatter("%(asctime)s@%(name)s[%(levelname)s] %(message)s"))
        handlers.append(file_handler)

    logging.basicConfig(level=logging.NOTSET,
                        handlers=handlers)


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

    parser.add_argument("-i", "--iteration",
                        nargs=1,
                        default=[0],
                        help='SCF iteration')

    # parser.add_argument("-v", "--verbose",
    #                     action="store_true",
    #                     default=False)
    parser.add_argument('--debug',
                        action='store_true',
                        default=False)
    parser.add_argument("-p", "--profile",
                        nargs='?',
                        const='pdf-report.stat',
                        help="do profiling")
    parser.add_argument('--yaml',
                        nargs=1,
                        default=None,
                        help='output yaml path')

    args = parser.parse_args()
    do_profile = args.profile
    is_debug = args.debug
    yaml_path = args.yaml[0] if args.yaml else None

    setup_logger(is_debug)
    # verbose = args.verbose
    # if verbose:
    #     logger.setLevel(logging.DEBUG)
    #     logger.info("info: debug on")
    # do_profile = args.profile

    # setting
    iteration = int(args.iteration[0])

    output_dir = args.output_dir
    pdfparam = pdf.PdfParam_H5()
    pdfparam.open(args.pdfresults_db)


    # run
    report = pdf.PdfReport(pdfparam=pdfparam, workdir=output_dir)
    if do_profile != None:
        import cProfile
        cProfile.runctx(
            "report.report(iteration=iteration)",
            globals(),
            locals(),
            filename=do_profile,
        )
    else:
        report.report(iteration=iteration)

    if yaml_path:
        data = report.get_report_data(iteration=iteration, source=args.pdfresults_db)
        report.save_report_yaml(data, yaml_path)


if __name__ == '__main__':
    if os.path.exists("config.ini"):
        logging.config.fileConfig("config.ini",
                                  disable_existing_loggers=False)
    main()

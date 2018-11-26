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
import logging
import logging.config
import traceback
try:
    import msgpack
except:
    import msgpack_pure as msgpack
#import pprint

import proteindf_bridge as bridge
import proteindf_tools as pdf

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
                        default='pdfresults.db',
                        help='ProteinDF results file')

    parser.add_argument('--calc-pop',
                        nargs='*',
                        choices=['mulliken', 'mk'],
                        default=argparse.SUPPRESS,
                        help='calc population')

    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    # setting from command-line
    pdfparam_path = args.pdfparam_path
    output = args.output
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    #is_little_endian = True

    do_calc_pop_types = [None]
    if 'calc_pop' in args:
        if len(args.calc_pop) == 0:
            do_calc_pop_types = ['mulliken']
        else:
            do_calc_pop_types = args.calc_pop

    try:
        # read ProteinDF parameters
        pdfparam = pdf.load_pdfparam(pdfparam_path)

        # setup DB
        entry = pdf.PdfArchive(output,
                               pdfparam=pdfparam)

        pdf2db = Pdf2Db(pdfparam, entry,
                        pop_types=do_calc_pop_types)
        pdf2db.regist()
    except:
        print('-'*60)
        traceback.print_exc(file=sys.stdout)
        #print(traceback.format_exc())
        print('-'*60)


class Pdf2Db(object):
    def __init__(self, pdfparam, entry, pop_types=[]):
        self._logger = logging.getLogger(__name__)

        #self._is_little_endian = True

        assert(isinstance(pdfparam, pdf.PdfParam))
        assert(isinstance(entry, pdf.PdfArchive))
        self._pdfparam = pdfparam
        self._entry = entry
        self._pop_types = pop_types

    def regist(self):
        # energy level
        self.set_energylevels('full')

        # OCC
        self.set_occ()

        # LCAO matrix
        self.set_lcao('last')

        # calc population
        if self._pdfparam.scf_converged:
            self.set_pop(self._pop_types)

    #def _get_is_little_endian(self):
    #    return self._is_little_endian

    #def _set_is_little_endian(self, value):
    #    assert(isinstance(value, bool))
    #    self._is_little_endian = value

    #is_little_endian = property(_get_is_little_endian, _set_is_little_endian)

    def set_energylevels(self, mode ='last'):
        """
        エネルギー準位をDBに格納する

        TODO: modeの実装
        """
        pdfparam = self._pdfparam
        entry = self._entry
        iterations = pdfparam.iterations

        for itr in range(1, iterations +1):
            v = pdf.Vector()
            # energy level
            for runtype in pdfparam.runtypes():
                path = pdfparam.get_energy_level_path(itr, runtype)
                self._logger.debug('load %s' % (path))
                if v.is_loadable(path):
                    v.load(path)
                    entry.set_energylevel(runtype, itr, v)

    def set_occ(self):
        """
        占有軌道情報をDBに格納する
        """
        pdfparam = self._pdfparam
        entry = self._entry

        v = pdf.Vector()
        for runtype in pdfparam.runtypes():
            path = pdfparam.get_occ_path(runtype)
            self._logger.debug('load %s' % (path))
            if v.is_loadable(path):
                v.load(path)
                entry.set_occupations(runtype, v)
            else:
                self._logger.warning('could not load %s' % (path))

    def set_lcao(self, mode ='last'):
        """
        LCAO行列をDBに格納する

        TODO: modeの実装
        """
        pdfparam = self._pdfparam
        entry = self._entry
        iterations = pdfparam.iterations

        start_itr = 1
        if mode == 'last':
            start_itr = iterations

        m = pdf.Matrix()
        for itr in range(start_itr, iterations +1):
            for runtype in pdfparam.runtypes():
                path = pdfparam.get_cmat_path(runtype, itr)
                self._logger.debug('load %s' % (path))
                if m.is_loadable(path):
                    m.load(path)
                    entry.set_lcao(runtype, itr, m)

    def set_pop(self, pop_types = [], iteration = 0):
        for pop_type in pop_types:
            if pop_type == 'mulliken':
                self._logger.info('calc pop mulliken')
                if iteration == 0:
                    iteration = self._pdfparam.iterations
                pop_matrix_path = 'mulliken_{}.mat'.format(iteration)
                pdf.run_pdf(['pop-mulliken', '-i', iteration, '-s', pop_matrix_path])
                m = pdf.Matrix()
                m.load(pop_matrix_path)
                self._entry.set_population('mulliken', iteration, m)
            else:
                self._logger.warning('unknown pop type: {}'.format(pop_type))


if __name__ == '__main__':
    logconfig_file = 'logconfig.ini'
    if os.path.exists(logconfig_file):
        logging.config.fileConfig(logconfig_file)
    main()

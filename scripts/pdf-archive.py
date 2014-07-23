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
archive ProteinDF results.
"""

import sys
import argparse
import logging
import traceback
try:
    import msgpack
except:
    import msgpack_pure as msgpack
#import pprint
    
import bridge
import pdf

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
    is_little_endian = True

    try:
        # read ProteinDF parameters
        pdfparam = pdf.load_pdfparam(pdfparam_path)
    
        # setup DB
        entry = pdf.PdfArchive(output,
                               pdfparam=pdfparam)

        pdf2db = Pdf2Db(pdfparam, entry)
        pdf2db.regist()
    except:
        print('-'*60)
        traceback.print_exc(file=sys.stdout)
        #print(traceback.format_exc())
        print('-'*60)
    

class Pdf2Db(object):
    def __init__(self, pdfparam, entry):
        nullHandler = bridge.NullHandler()
        self._logger = logging.getLogger(__name__)
        self._logger.addHandler(nullHandler)

        self._is_little_endian = True

        assert(isinstance(pdfparam, pdf.PdfParam))
        assert(isinstance(entry, pdf.PdfArchive))
        self._pdfparam = pdfparam
        self._entry = entry
        
    def regist(self):
        # energy level
        self.set_energylevels('full')

        # OCC
        self.set_occ()
        
        # LCAO matrix
        self.set_lcao('last')

    def _get_is_little_endian(self):
        return self._is_little_endian
    
    def _set_is_little_endian(self, value):
        assert(isinstance(value, bool))
        self._is_little_endian = value

    is_little_endian = property(_get_is_little_endian, _set_is_little_endian)
        
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
                if v.is_loadable(path, self.is_little_endian):
                    v.load(path, self.is_little_endian)
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
            if v.is_loadable(path, self.is_little_endian):
                v.load(path, self.is_little_endian)
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
                if m.is_loadable(path, self.is_little_endian):
                    m.load(path, self.is_little_endian)
                    entry.set_lcao(runtype, itr, m)
                
    
if __name__ == '__main__':
    main()


#!/usr/bin/env python

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

import sys, os
import shutil
from distutils.core import setup

shutil.copyfile('scripts/pdfcmd.py', 'scripts/pdf')

setup(name='ProteinDF_pytools',
      version='0.1',
      description='python scripts for ProteinDF package',
      author='Toshiyuki HIRANO',
      author_email='hiracchi@gmail.com',
      url='http://proteindf.github.io/',
      packages=['pdfpytools'],
      scripts=[
          'scripts/pdf',
          'scripts/pdf-archive.py',
          'scripts/pdf-env.py',
#          'scripts/pdf-info.py',
	  'scripts/pdf-info-geom.py',
          'scripts/pdf-info-orb.py',
	  'scripts/pdf-info-xyz.py',
          'scripts/pdf-make-basis2.py',
          'scripts/pdf-plot-basisset.py',
          'scripts/pdf-plot-elevel.py',
          'scripts/pdf-plot-mat.py',
          'scripts/pdf-plot-decaymat.py',
          'scripts/pdf-report.py',
	  'scripts/pdf-reg-harris.py',
          'scripts/pdf-test.py',
          'scripts/pdf-test-eri.py',
          'scripts/xyz2gau.py',
          'scripts/xyz2pdf.py',
          'scripts/g0xeri2mpac.py',
          'scripts/pdf-pop-lasso.py',
          'scripts/pdf-pop-ridge.py',
          'scripts/pdf-pop-elasticnet.py',
          'scripts/pdf-pop-resp.py',
          'scripts/pdf-pop-rrms.py',
          'scripts/pdf-set-charges.py'
      ],
      data_files=[
          ('data', ['data/basis2.cache'])
      ]
)


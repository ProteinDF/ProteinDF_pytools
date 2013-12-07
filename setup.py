#!/usr/bin/env python

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

import sys, os
import shutil
from distutils.core import setup

#if not os.path.exists('scripts'):
#    os.makedirs('scripts')
shutil.copyfile('scripts/pdfcmd.py', 'scripts/pdf')

setup(name='pdfutils',
      version='0.1',
      description='ProteinDF python scripts',
      author='Toshiyuki HIRANO',
      author_email='t-hirano@iis.u-tokyo.ac.jp',
      url='http://satolab.iis.u-tokyo.ac.jp/',
      packages=['pdf'],
      scripts=[
        'scripts/mpac2yml.py',
        'scripts/yml2mpac.py',
        'scripts/pdf',
        'scripts/pdf-archive.py',
        'scripts/pdf-env.py',
	'scripts/pdf-info-geom.py',
        'scripts/pdf-info-orb.py',
	'scripts/pdf-info-xyz.py',
        'scripts/pdf-plot-elevel.py',
        'scripts/pdf-report.py',
	'scripts/pdf-reg-harris.py',
        'scripts/pdf-test.py',
        'scripts/pdf-test-eri.py',
        'scripts/g0xeri2mpac.py',
        ]
      )


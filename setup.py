#!/usr/bin/env python

import sys, os
from distutils.core import setup

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
        'scripts/pdfcmd.py',
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


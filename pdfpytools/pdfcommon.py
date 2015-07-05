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

import io
import sys
import os
import subprocess
import shlex
import tempfile
import logging
import traceback

try:
    import msgpack
except:
    import msgpack_pure as msgpack

import pdfbridge as bridge
import pdfpytools as pdf

epsilon = 1.0E-10 # 計算機イプシロン
error = 1.0E-5 # 許容誤差

module_logger = logging.getLogger(__name__)
module_logger.addHandler(logging.NullHandler())
module_logger.setLevel(logging.INFO)

#class NullHandler(logging.Handler):
#    def emit(self, record):
#        pass

def pdf_home():
    """
    return PDF_HOME environment parameter value
    """
    answer = os.environ.get('PDF_HOME', '')
    return answer

def get_default_pdfparam():
    """
    defaultのpdfparamを返す
    """
    # make temp dir & filepath
    tempfile_fd, tempfile_path = tempfile.mkstemp()
    os.close(tempfile_fd)
    
    # 一時ファイルの初期化情報を読取る
    pdf.run_pdf(['init-param', '-v', '-o', tempfile_path])
    f = open(tempfile_path, "rb")
    tempdata = msgpack.unpackb(f.read())
    tempdata = bridge.Utils.byte2str(tempdata)
    f.close()

    # remove temp
    os.remove(tempfile_path)
    
    pdfparam = pdf.PdfParam(tempdata)
    pdfparam.step_control = 'create integral guess scf'

    pdfparam.guess = 'harris'
    pdfparam.orbital_independence_threshold = 0.007
    pdfparam.orbital_independence_threshold_canonical = 0.007
    pdfparam.orbital_independence_threshold_lowdin = 0.007
    pdfparam.scf_acceleration = 'damping'
    pdfparam.scf_acceleration_damping_factor = 0.85
    pdfparam.convergence_threshold_energy = 1.0E-4
    pdfparam.convergence_threshold = 1.0E-3
    pdfparam.scf_acceleration_damping_damping_type = 'density_matrix'
    pdfparam.xc_functional = "b3lyp"
    pdfparam.j_engine = "CD"
    pdfparam.k_engine = "CD"
    pdfparam.xc_engine = "grid"
    pdfparam.gridfree_orthogonalize_method = "canonical"

    return pdfparam

def set_basisset(pdfparam,
                 basisset_name_ao = "DZVP2",
                 basisset_name_rij = "DZVP2",
                 basisset_name_rixc = "DZVP2",
                 basisset_name_gridfree = "cc-pVDZ-SP"):
    """
    pdfparamにbasissetを設定する
    """
    assert(isinstance(pdfparam, pdf.PdfParam))

    basis2 = Basis2()
    atoms = ['C', 'H', 'N', 'O', 'S']
    basisset = {}

    if basisset_name_ao == basisset_name_gridfree:
        pdfparam.gridfree_dedicated_basis = False
    else:
        pdfparam.gridfree_dedicated_basis = True

    for atom in atoms:
        basisset_ao = basis2.get_basisset('O-{}.{}'.format(basisset_name_ao, atom))
        basisset_j = basis2.get_basisset_j('A-{}.{}'.format(basisset_name_rij, atom))
        basisset_xc = basis2.get_basisset_xc('A-{}.{}'.format(basisset_name_rixc, atom))
        basisset_gf = basis2.get_basisset('O-{}.{}'.format(basisset_name_gridfree, atom))

        pdfparam.set_basisset(atom, basisset_ao)
        pdfparam.set_basisset_j(atom, basisset_j)
        pdfparam.set_basisset_xc(atom, basisset_xc)
        pdfparam.set_basisset_gridfree(atom, basisset_gf)

    return pdfparam

def run_pdf(subcmd):
    """
    run ProteinDF command
    """
    module_logger.info("run_pdf({})".format(subcmd))

    try:
        if isinstance(subcmd, list):
            subcmd_tmp = [str(x) for x in subcmd ]
            subcmd = " ".join(subcmd_tmp)
    except:
        print(subcmd)
        raise

    cmd = os.path.join(pdf_home(), "bin", "pdf") + " " + subcmd
    module_logger.debug("run: {0}".format(cmd))

    p = pdf.Process()
    return_code = p.cmd(cmd).commit()
    
    module_logger.debug('return code={}'.format(return_code))
    if return_code != 0:
        sys.stderr.write('Failed to execute command: %s' % cmd)
        sys.stderr.write('return code = {}'.format(return_code))
        module_logger.critical('Failed to execute command: %s' % cmd)
        module_logger.critical('return code = {}'.format(return_code))
        raise
    
    
def mpac2py(path):
    """
    load message pack binary file to python dictionary data
    """
    assert(isinstance(path, str) == True)

    f = open(path, "rb")
    contents = f.read()
    data = bridge.Utils.byte2str(msgpack.unpackb(contents))
    f.close()

    return data

def load_pdfparam(pdfparam_path='pdfparam.mpac'):
    data = mpac2py(pdfparam_path)
    param = pdf.PdfParam(data)

    return param

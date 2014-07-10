#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import subprocess
import logging
import traceback

try:
    import msgpack
except:
    import msgpack_pure as msgpack

import bridge
import pdf

epsilon = 1.0E-10 # 計算機イプシロン
error = 1.0E-5 # 許容誤差

class NullHandler(logging.Handler):
    def emit(self, record):
        pass

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
    pdfparam = pdf.PdfParam()
    pdfparam.step_control = 'create integral guess scf'

    pdfparam.guess = 'harris'
    pdfparam.orbital_independence_threshold = 0.0
    pdfparam.orbital_independence_threshold_canonical = 0.0
    pdfparam.orbital_independence_threshold_lowdin = 0.0
    pdfparam.scf_acceleration = 'damping'
    pdfparam.scf_acceleration_damping_factor = 0.85
    pdfparam.convergence_threshold_energy = 1.0E-4
    pdfparam.convergence_threshold = 1.0E-3
    pdfparam.convergence_target = 'density_matrix'
    pdfparam.xc_functional = "b3lyp"
    pdfparam.j_engine = "CD"
    pdfparam.k_engine = "CD"
    pdfparam.xc_engine = "gridfree_CD"

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
    if isinstance(subcmd, str):
        subcmd = [subcmd]

    nullHandler = NullHandler()
    logger = logging.getLogger(__name__)
    logger.addHandler(nullHandler)

    cmd = os.path.join(pdf_home(), "bin", "pdf")
    cmdlist = [cmd]
    cmdlist.extend(subcmd)
    logger.debug("run: {0}".format(cmdlist))

    try:
        subprocess.check_call(cmdlist)
    except:
        print('-'*60)
        traceback.print_exc(file=sys.stdout)
        #print(traceback.format_exc())
        print('-'*60)

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

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import subprocess
import tempfile
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

module_logger = logging.getLogger(__name__)

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
    print('get_default_pdfparam()')
    # 一時ファイル名を取得する
    tmpfile_fd, tmpfile_path = tempfile.mkstemp()
    os.close(tmpfile_fd)
    # 一時ファイルの初期化情報を読取る
    pdf.run_pdf(['init-param', '-o', tmpfile_path])
    f = open(tmpfile_path, "rb")
    tmpdata = msgpack.unpackb(f.read())
    tmpdata = bridge.Utils.byte2str(tmpdata)
    f.close()
    # remove
    os.remove(tmpfile_path)
    
    pdfparam = pdf.PdfParam(tmpdata)
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
    pdfparam.xc_engine = "gridfree_CD"
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
    #handler = logging.StreamHandler()
    #handler.setLevel(logging.DEBUG)
    #module_logger.setLevel(logging.DEBUG)
    #logger.addHandler(handler)

    module_logger.info("pdf.run_pdf >> {}".format(subcmd))

    if isinstance(subcmd, str):
        subcmd = [subcmd]

    for i, c in enumerate(subcmd):
        subcmd[i] = str(c)
        
    cmd = os.path.join(pdf_home(), "bin", "pdf")
    cmdlist = [cmd]
    cmdlist.extend(subcmd)
    module_logger.critical("run: {0}".format(cmdlist))
    
    #try:
    #    p = subprocess.Popen(cmdlist,
    #                         stdout=subprocess.PIPE,
    #                         stderr=subprocess.PIPE)
    #    stdout, stderr = p.communicate()
    #    if stdout:
    #        sys.stdout.write(stdout)
    #        logger.info(stdout)
    #    if stderr:
    #        sys.stderr.write(stderr)
    #        logger.error(stderr)
    #except:
    #    print('-'*60)
    #    traceback.print_exc(file=sys.stdout)
    #    #print(traceback.format_exc())
    #    print('-'*60)

    subproc_args = {'stdin': None,
                    'stdout': subprocess.PIPE,
                    'stderr': subprocess.STDOUT,
                    'shell': False,
                    'universal_newlines': True
                    }
    try:
        proc = subprocess.Popen(cmdlist, **subproc_args)
    except OSError as e:
        module_logger.critical('Failed to execute command: %s' % cmdlist)
        module_logger.critical(str(errno.errorcode[e.errno]))
        module_logger.critical(str(os.strerror(e.errno)))
        raise e

    if proc.stdout != None:
        while True:
            line = proc.stdout.readline()
            if not line:
                break
            module_logger.info(line)
        
    return_code = proc.wait()
    (stdouterr, stdin) = (proc.stdout, proc.stdin)
    module_logger.info('return code={}'.format(return_code))
    if return_code != 0:
        print('Failed to execute command: %s' % cmdlist)
        print('return code = '.format(return_code))
        module_logger.critical('Failed to execute command: %s' % cmdlist)
        module_logger.critical('return code = '.format(return_code))
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

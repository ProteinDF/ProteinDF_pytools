#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import subprocess
import logging

try:
    import msgpack
except:
    import msgpack_pure as msgpack

from pdfparam import PdfParam

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

def pdf_setup(workdir = "."):
    """
    setup to run ProteinDF
    """
    dirs = ['fl_Work']
    for d in dirs:
        path = os.path.join(workdir, d)
        if not os.path.exists(path):
            os.mkdir(path)

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
    print(cmdlist)
    subprocess.check_call(cmdlist)
    
    
def mpac2py(path):
    """
    load message pack binary file to python dictionary data
    """
    assert(isinstance(path, str) == True)
    
    f = open(path, "rb")
    contents = f.read()
    data = msgpack.unpackb(contents)
    f.close()

    return data

def load_pdfparam(pdfparam_path='pdfparam.mpac'):
    data = mpac2py(pdfparam_path)
    param = PdfParam(data)
    
    return param

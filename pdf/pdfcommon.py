#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    import msgpack
except:
    import msgpack_pure as msgpack

from pdfparam import PdfParam

epsilon = 1.0E-10 # 計算機イプシロン

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

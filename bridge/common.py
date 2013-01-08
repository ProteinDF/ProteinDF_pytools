#!/usr/bin/env python
# -*- coding: utf-8 -*-

try:
    import msgpack
except:
    import msgpack_pure as msgpack

    
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

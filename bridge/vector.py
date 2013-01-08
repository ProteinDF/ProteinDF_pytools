#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import struct
import copy
import numpy
from types import *

class Vector(object):
    """
    >>> a = Vector(10)
    >>> len(a)
    10
    >>> a.resize(20)
    >>> len(a)
    20
    >>> a.resize(5)
    >>> len(a)
    5
    >>> a.set(1, 3.14)
    >>> (a.get(0) - 0.00 < 1.0E-10)
    True
    >>> (a.get(1) - 3.14 < 1.0E-10)
    True
    >>> (a[1] - 3.14 < 1.0E-10)
    True
    >>> a[2] = 1.73
    >>> (a[2] - 1.73 < 1.0E-10)
    True
    """

    __header_struct_little_endian = "<i"
    __header_struct_big_endian = ">i" 
    __body_struct_little_endian = "<d"
    __body_struct_big_endian = ">d"
    
    def __init__(self, obj =[]):
        if isinstance(obj, int):
            size = obj
            self._data = numpy.array([0.0 for x in range(size)])
        elif isinstance(obj, list):
            self._data = numpy.array(obj)
        elif isinstance(obj, dict):
            size = raw_data.get('size', 0)
            buf = raw_data.get('data', None)
            if (buf != None):
                self._data = numpy.array([0.0 for x in range(size)])
        elif isinstance(obj, Vector):
            self._data = obj._data
        else:
            raise TypeError
            
    def size(self):
        return self.__len__()
                
    def resize(self, new_size):
        new_data = numpy.array([0.0 for x in range(new_size)])
        for i in range(min(self.size(), new_size)):
            new_data[i] = self._data[i]
        self._data = new_data

    def get(self, index):
        return self._data[index]
        
    def set(self, index, value):
        self._data[index] = value

    def to_list(self):
        return self._data.tolist()

    def get_buffer(self):
        return buffer(self._data.tostring())

    def set_buffer(self, buf):
        self._data = numpy.fromstring(buf)
        
    def __get_header_struct(self, is_little_endian):
        if is_little_endian:
            return self.__header_struct_little_endian
        else:
            return self.__header_struct_big_endian

    def __get_body_struct(self, is_little_endian):
        if is_little_endian:
            return self.__body_struct_little_endian
        else:
            return self.__body_struct_big_endian
        
    def __str__(self):
        output = ''
        for order in range(0, len(self), 10):
            output += '\n'
            for j in range(order, min(order +10, len(self))):
                output += '   %5d th' % (j +1)
            output += '\n'
            for j in range(order, min(order +10, len(self))):
                output += '-----------'
            output += '----\n\n'
            for j in range(order, min(order +10, len(self))):
                output += ' %10.6lf' % (self[j])
            output += '\n'
        return output
  
    def __add__(self, other):
        assert isinstance(other, Vector)
        assert (len(self) == len(other))

        answer = self.copy()
        answer += other
        return answer

    def __iadd__(self, other):
        assert isinstance(other, Vector)
        assert (len(self) == len(other))

        self._data += other._data
        return self

    def __sub__(self, other):
        assert isinstance(other, Vector)

        answer = self.copy()
        answer -= other
        return answer

    def __isub__(self, other):
        assert isinstance(other, Vector)
        
        self._data -= other._data
        return self

    def __mul__(self, other):
        answer = copy.deepcopy(self)
        answer *= other
        return answer
    
    def __imul__(self, other):
        if type(other) is FloatType:
            for i in range(len(self)):
                self[i] *= other
        else:
            raise
        return self

    def __len__(self):
        return len(self._data)

    def __getitem__(self, key):
        return self._data[key]

    def __setitem__(self, key, value):
        self._data[key] = value
        
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()


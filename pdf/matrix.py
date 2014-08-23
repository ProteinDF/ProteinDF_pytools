#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import struct
import logging

import bridge

class Matrix(bridge.Matrix):
    """
    """
    __header_struct_little_endian = "<iii"
    __header_struct_big_endian = ">iii"
    __body_struct_little_endian = "<d"
    __body_struct_big_endian = ">d"

    def __init__(self, *args, **kwargs):
        self._logger = logging.getLogger(__name__)
        super(Matrix, self).__init__(self, *args, **kwargs)

    @classmethod
    def __get_header_struct(cls, is_little_endian):
        if is_little_endian:
            return cls.__header_struct_little_endian
        else:
            return cls.__header_struct_big_endian

    @classmethod
    def __get_body_struct(cls, is_little_endian):
        if is_little_endian:
            return cls.__body_struct_little_endian
        else:
            return cls.__body_struct_big_endian

    @classmethod
    def is_loadable(cls, file_path, is_little_endian = True):
        answer = False
        row = None
        col = None
        if os.path.isfile(file_path):
            fin = open(file_path, "rb")
            header_struct = cls.__get_header_struct(is_little_endian)
            size_of_header = struct.calcsize(header_struct)
            data = fin.read(size_of_header)
            fin.close()
            header = struct.unpack(header_struct, data[0: size_of_header])
            matrix_type = header[0]
            row = header[1]  
            col = header[2]
            if (matrix_type == 0):
                answer = True
        return (answer, row, col)

    def load(self, file_path, is_little_endian = True):
        if os.path.isfile(file_path):
            fin = open(file_path, "rb")
            # read header
            header_struct = self.__get_header_struct(is_little_endian)
            size_of_header = struct.calcsize(header_struct);
            header_bin = fin.read(size_of_header)
            header = struct.unpack(header_struct, header_bin)
            matrix_type = header[0]
            row = header[1]
            col = header[2]
            assert(matrix_type == 0)

            # read contents
            body_struct = self.__get_body_struct(is_little_endian)
            self.resize(row, col)
            size_of_double = struct.calcsize(body_struct);
            for r in range(row):
                for c in range(col):
                    value = struct.unpack(body_struct, fin.read(size_of_double))
                    self.set(r, c, value[0])
            fin.close()
        else:
            self._logger.error("file not found: %s" % (file_path))
            
    def save(self, file_path, is_little_endian = True):
        row = self.rows
        col = self.cols
        matrix_type = 0

        fout = open(file_path, "wb")
        # write header
        header_struct = self.__get_header_struct(is_little_endian)
        header = struct.pack(header_struct,
                             matrix_type, row, col)
        fout.write(header)
        
        # write elements
        body_struct = self.__get_body_struct(is_little_endian)
        for r in range(row):
            for c in range(col):
                value = self.get(r, c)
                value_str = struct.pack(body_struct, value)
                fout.write(value_str)
        
        fout.close()


class SymmetricMatrix(bridge.SymmetricMatrix):
    """
    """
    __header_struct_little_endian = "<iii"
    __header_struct_big_endian = ">iii"
    __body_struct_little_endian = "<d"
    __body_struct_big_endian = ">d"

    def __init__(self, *args, **kwargs):
        self._logger = logging.getLogger(__name__)
        super(SymmetricMatrix, self).__init__(self, *args, **kwargs)

    @classmethod
    def __get_header_struct(cls, is_little_endian):
        if is_little_endian:
            return cls.__header_struct_little_endian
        else:
            return cls.__header_struct_big_endian

    @classmethod
    def __get_body_struct(cls, is_little_endian):
        if is_little_endian:
            return cls.__body_struct_little_endian
        else:
            return cls.__body_struct_big_endian

    @classmethod
    def is_loadable(cls, file_path, is_little_endian = True):
        answer = False
        row = None
        col = None
        if os.path.isfile(file_path):
            fin = open(file_path, "rb")
            header_struct = cls.__get_header_struct(is_little_endian)
            size_of_header = struct.calcsize(header_struct);
            data = fin.read(size_of_header)
            fin.close()
            header = struct.unpack(header_struct, data[0: size_of_header])
            matrix_type = header[0]
            row = header[1]  
            col = header[2]
            if (matrix_type == 2):
                answer = True
        return (answer, row, col)

    def load(self, file_path, is_little_endian = True):
        if os.path.isfile(file_path):
            fin = open(file_path, "rb")
            # read header
            header_struct = self.__get_header_struct(is_little_endian)
            size_of_header = struct.calcsize(header_struct);
            header = struct.unpack(header_struct,
                                   fin.read(size_of_header))
            matrix_type = header[0]
            row = header[1]  
            col = header[2]
            dim = row
            assert(matrix_type == 2)
            assert(row == col)
            self.resize(dim)

            # body
            body_struct = self.__get_body_struct(is_little_endian)
            size_of_double = struct.calcsize(body_struct);
            for r in range(dim):
                for c in range(r +1):
                    value = struct.unpack(body_struct, fin.read(size_of_double))
                    self.set(r, c, value[0])
            fin.close()
        else:
            self._logger.error("file not found: %s" % (file_path))
            
    def save(self, file_path, is_little_endian = True):
        dim = self.rows
        assert(dim == self.cols)
        matrix_type = 2

        fout = open(file_path, "wb")
        # write header
        header_struct = self.__get_header_struct(is_little_endian)
        header = struct.pack(header_struct,
                             matrix_type, dim, dim)
        fout.write(header)
        
        # write elements
        body_struct = self.__get_body_struct(is_little_endian)
        for r in range(dim):
            for c in range(r +1):
                value = self.get(r, c)
                value_str = struct.pack(body_struct, value)
                fout.write(value_str)
        
        fout.close()
        
        

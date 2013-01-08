#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import struct

import bridge

class Matrix(bridge.Matrix):
    """
    """
    def __init__(self, *args, **kwargs):
        bridge.Matrix.__init__(self, *args, **kwargs)

    def is_loadable(self, file_path, is_little_endian = True):
        fin = open(file_path, "rb")
        header_struct = self.__get_header_struct(is_little_endian)
        size_of_header = struct.calcsize(header_struct)
        data = fin.read(size_of_header)
        fin.close()
        header = struct.unpack(header_struct, data[0: size_of_header])
        matrix_type = header[0]
        row = header[1]  
        col = header[2]
        if (matrix_type == 0):
            return True
        else:
            return False
        
    def load(self, file_path, is_little_endian = True):
        if (os.path.exists(file_path) == True):
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
            print("file not found: %s" % (file_path))

    def save(self, file_path, is_little_endian = True):
        row = self.get_num_of_rows()
        col = self.get_num_of_cols()
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
    def __init__(self, *args, **kwargs):
        bridge.SymmetricMatrix.__init__(self, *args, **kwargs)

    def is_loadable(self, file_path, is_little_endian = True):
        fin = open(file_path, "rb")
        header_struct = self.__get_header_struct(is_little_endian)
        size_of_header = struct.calcsize(header_struct);
        data = fin.read(size_of_header)
        fin.close()
        header = struct.unpack(header_struct, data[0: size_of_header])
        matrix_type = header[0]
        row = header[1]  
        col = header[2]
        if (matrix_type == 2):
            return True
        else:
            return False

    def load(self, file_path, is_little_endian = True):
        if (os.path.exists(file_path) == True):
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
            self.get_num_of_rows()

            # body
            body_struct = self.__get_body_struct(is_little_endian)
            size_of_double = struct.calcsize(body_struct);
            for r in xrange(dim):
                for c in xrange(r +1):
                    value = struct.unpack(body_struct, fin.read(size_of_double))
                    self.set(r, c, value[0])
            fin.close()
        else:
            print("file not found: %s" % (file_path))
            
    def save(self, file_path, is_little_endian = True):
        dim = self.get_num_of_rows()
        assert(dim == self.get_num_of_cols())
        matrix_type = 2

        fout = open(file_path, "wb")
        # write header
        header_struct = self.__get_header_struct(is_little_endian)
        header = struct.pack(header_struct,
                             matrix_type, dim, dim)
        fout.write(header)
        
        # write elements
        body_struct = self.__get_body_struct(is_little_endian)
        for r in xrange(dim):
            for c in xrange(r +1):
                value = self.get(r, c)
                value_str = struct.pack(body_struct, value)
                fout.write(value_str)
        
        fout.close()
        
        

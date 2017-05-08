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

import os
import struct
import logging

import pdfbridge as bridge

class Matrix(bridge.Matrix):
    """
    """
    __header_struct_little_endian = "<iii"
    __header_struct_big_endian = ">iii"
    __body_struct_little_endian = "<d"
    __body_struct_big_endian = ">d"

    def __init__(self, *args, **kwargs):
        self._logger = logging.getLogger(__name__)
        if kwargs.get('debug'):
            self._logger.addHandler(logging.StreamHandler())
            self._logger.setLevel(logging.DEBUG)
        else:
            self._logger.addHandler(logging.NullHandler())
            self._logger.setLevel(logging.INFO)

        super(Matrix, self).__init__(*args, **kwargs)

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
        return answer

    
    @classmethod
    def get_size(cls, file_path, is_little_endian = True):
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
            if (matrix_type == 0):
                row = header[1]  
                col = header[2]
        return (row, col)

    
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
        super(SymmetricMatrix, self).__init__(*args, **kwargs)

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
        return answer

    
    @classmethod
    def get_size(cls, file_path, is_little_endian = True):
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
            if (matrix_type == 2):
                row = header[1]  
                col = header[2]
        return (row, col)

    
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
        
        

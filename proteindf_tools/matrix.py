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

import proteindf_bridge as bridge
import os
import struct
import logging
logger = logging.getLogger(__name__)


class Matrix(bridge.Matrix):
    """
    """
    __header_struct_little_endian = "<iii"
    __header_struct_big_endian = ">iii"
    __body_struct_little_endian = "<d"
    __body_struct_big_endian = ">d"

    def __init__(self, *args, **kwargs):
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
    def find_header_struct(cls, file_path):
        def get_header(file_path, header_struct="<iii"):
            fin = open(file_path, "rb")
            size_of_header = struct.calcsize(header_struct)
            data = fin.read(size_of_header)
            fin.close()

            header = struct.unpack(header_struct, data[0: size_of_header])
            matrix_type = header[0]
            row = header[1]
            col = header[2]
            size_of_data = struct.calcsize("d") * row * col
            size_of_file = size_of_header + size_of_data
            return (matrix_type, row, col, size_of_file)

        answer = False
        endian = None
        matrix_type = None
        row = None
        col = None
        header_struct = None
        if os.path.isfile(file_path):
            filesize = os.path.getsize(file_path)
            endian_list = ["=", "<", ">"]
            header_struct_list = ["iii", "bii", "ill", "bll"]
            for endian in endian_list:
                for header_struct in header_struct_list:
                    (matrix_type, row, col, chk_filesize) = get_header(
                        file_path, endian + header_struct)
                    answer = (filesize == chk_filesize)
                    if answer:
                        break
                if answer:
                    break

        return (endian, header_struct, matrix_type, row, col)

    @classmethod
    def is_loadable(cls, file_path):
        answer = False
        (endian, header_struct, matrix_type, row,
         col) = cls.find_header_struct(file_path)
        if endian != None:
            answer = True

        return answer

    @classmethod
    def get_size(cls, file_path):
        (endian, header_struct, matrix_type, row,
         col) = cls.find_header_struct(file_path)
        return (row, col)

    def load(self, file_path):
        """Load the matrix. True is returned if the reading is successful.
        """
        if os.path.isfile(file_path):
            (endian, header_struct, matrix_type, row,
             col) = self.find_header_struct(file_path)

            with open(file_path, "rb") as fin:
                # read header
                endian_header_struct = endian + header_struct
                size_of_header = struct.calcsize(endian_header_struct)
                header_bin = fin.read(size_of_header)

                # read contents
                body_struct = endian + "d"
                self.resize(row, col)
                size_of_double = struct.calcsize(body_struct)
                if matrix_type == 0:
                    # RSFD type
                    for r in range(row):
                        for c in range(col):
                            value = struct.unpack(
                                body_struct, fin.read(size_of_double))
                            self.set(r, c, value[0])
                elif matrix_type == 1:
                    # CSFD type
                    for c in range(col):
                        for r in range(row):
                            value = struct.unpack(
                                body_struct, fin.read(size_of_double))
                            self.set(r, c, value[0])
                else:
                    logger.critical(
                        "file type mismatch: type={}".format(matrix_type))
                    raise

            return True
        else:
            logger.error("file not found: %s" % (file_path))
            return False

    def save(self, file_path, matrix_type=1):
        row = self.rows
        col = self.cols

        fout = open(file_path, "wb")
        # write header
        header_struct = "=bii"
        header = struct.pack(header_struct,
                             matrix_type, row, col)
        fout.write(header)

        # write elements
        body_struct = "=d"
        if matrix_type == 0:
            for r in range(row):
                for c in range(col):
                    value = self.get(r, c)
                    value_str = struct.pack(body_struct, value)
                    fout.write(value_str)
        elif matrix_type == 1:
            for c in range(col):
                for r in range(row):
                    value = self.get(r, c)
                    value_str = struct.pack(body_struct, value)
                    fout.write(value_str)
        else:
            logger.error("unknown matrix type: {}".format(matrix_type))

        fout.close()


class SymmetricMatrix(bridge.SymmetricMatrix):
    """
    """
    __header_struct_little_endian = "<iii"
    __header_struct_big_endian = ">iii"
    __body_struct_little_endian = "<d"
    __body_struct_big_endian = ">d"

    def __init__(self, *args, **kwargs):
        logger = logging.getLogger(__name__)
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
    def find_header_struct(cls, file_path):
        def get_header(file_path, header_struct="<iii"):
            fin = open(file_path, "rb")
            size_of_header = struct.calcsize(header_struct)
            data = fin.read(size_of_header)
            fin.close()

            header = struct.unpack(header_struct, data[0: size_of_header])
            matrix_type = header[0]
            row = header[1]
            col = header[2]
            size_of_data = struct.calcsize("d") * row * (row + 1) / 2
            size_of_file = size_of_header + size_of_data
            return (matrix_type, row, col, size_of_file)

        answer = False
        endian = None
        matrix_type = None
        row = None
        col = None
        header_struct = None
        if os.path.isfile(file_path):
            filesize = os.path.getsize(file_path)
            endian_list = ["=", "<", ">"]
            header_struct_list = ["iii", "bii", "ill", "bll"]
            for endian in endian_list:
                for header_struct in header_struct_list:
                    (matrix_type, row, col, chk_filesize) = get_header(
                        file_path, endian + header_struct)
                    answer = (filesize == chk_filesize)
                    if answer:
                        break
                if answer:
                    break

        return (endian, header_struct, matrix_type, row, col)

    @classmethod
    def is_loadable(cls, file_path):
        answer = False
        (endian, header_struct, matrix_type, row,
         col) = cls.find_header_struct(file_path)
        if endian != None:
            answer = True

        return answer

    @classmethod
    def get_size(cls, file_path):
        (endian, header_struct, matrix_type, row,
         col) = cls.find_header_struct(file_path)
        return (row, col)

    def load(self, file_path):
        if os.path.isfile(file_path):
            (endian, header_struct, matrix_type, row,
             col) = self.find_header_struct(file_path)

            fin = open(file_path, "rb")
            # read header
            endian_header_struct = endian + header_struct
            size_of_header = struct.calcsize(endian_header_struct)
            header_bin = fin.read(size_of_header)
            dim = row
            assert(row == col)
            self.resize(dim)

            # body
            body_struct = endian + "d"
            size_of_double = struct.calcsize(body_struct)
            for r in range(dim):
                for c in range(r + 1):
                    value = struct.unpack(
                        body_struct, fin.read(size_of_double))
                    self.set(r, c, value[0])
            fin.close()
        else:
            logger.error("file not found: %s" % (file_path))

    def save(self, file_path, is_little_endian=True):
        dim = self.rows
        assert(dim == self.cols)
        matrix_type = 2

        fout = open(file_path, "wb")
        # write header
        header_struct = "=bii"
        header = struct.pack(header_struct,
                             matrix_type, dim, dim)
        fout.write(header)

        # write elements
        body_struct = "=d"
        for r in range(dim):
            for c in range(r + 1):
                value = self.get(r, c)
                value_str = struct.pack(body_struct, value)
                fout.write(value_str)

        fout.close()

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

import proteindf_bridge as bridge


class Vector(bridge.Vector):
    """
    """
    def __init__(self, obj = []):
        self._logger = logging.getLogger(__name__)
        super(Vector, self).__init__(obj)

    @classmethod
    def find_header_struct(cls, file_path):
        def get_header(file_path, header_struct="=i"):
            fin = open(file_path, "rb")
            size_of_header = struct.calcsize(header_struct)
            data = fin.read(size_of_header)
            fin.close()

            header = struct.unpack(header_struct, data[0: size_of_header])
            size = header[0]
            size_of_data = struct.calcsize("d") * size;
            size_of_file = size_of_header + size_of_data
            return (size, size_of_file)

        answer = False
        endian = None
        size = None
        header_struct = None
        if os.path.isfile(file_path):
            filesize = os.path.getsize(file_path)
            endian_list = ["=", "<", ">"]
            header_struct_list = ["i", "l"]
            for endian in endian_list:
                for header_struct in header_struct_list:
                    (size, chk_filesize) = get_header(file_path, endian + header_struct)
                    answer = (filesize == chk_filesize)
                    if answer:
                        break
                if answer:
                    break

        return (endian, header_struct, size)

    @classmethod
    def is_loadable(cls, file_path):
        answer = False
        (endian, header_struct, size) = cls.find_header_struct(file_path)
        if endian != None:
            answer = True

        return answer

    def load(self, file_path):
        if os.path.isfile(file_path):
            (endian, header_struct, size) = self.find_header_struct(file_path)

            data = open(file_path, "rb").read()
            # read header
            header_struct = endian + header_struct
            start = 0
            size_of_header = struct.calcsize(header_struct);
            header = struct.unpack(header_struct, data[start: start + size_of_header])
            start += size_of_header
            size = header[0]
            assert(size >= 0)

            # read contents
            body_struct = endian + "d"
            self.resize(size)
            size_of_double = struct.calcsize(body_struct);
            for i in range(len(self)):
                value = struct.unpack(body_struct, data[start: start + size_of_double])
                self[i] = float(value[0])
                start += size_of_double
        else:
            self._logger.error("file not found: %s" % (file_path))

    def save(self, file_path):
        size = len(self)

        fout = open(file_path, "wb")
        # write header
        header_struct = "l"
        header = struct.pack(header_struct, size)
        fout.write(header)

        # write elements
        body_struct = "d"
        for i in range(size):
            value_str = struct.pack(body_struct, self[i])
            fout.write(value_str)

        fout.close()

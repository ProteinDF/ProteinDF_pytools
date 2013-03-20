#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import struct
import logging

import bridge

class Vector(bridge.Vector):
    """
    """
    def __init__(self, obj = []):
        self._logger = logging.getLogger(__name__)
        bridge.Vector.__init__(self, obj)

    def is_loadable(self, file_path, is_little_endian = True):
        answer = False
        if os.path.isfile(file_path):
            fin = open(file_path, "rb")
            header_struct = self.__get_header_struct(is_little_endian)
            size_of_header = struct.calcsize(header_struct);
            data = fin.read(size_of_header)
            fin.close()
            header = struct.unpack(header_struct, data[0: size_of_header])
            size = header[0]
            if (size >= 0):
                answer = True
        return answer

    def load(self, file_path, is_little_endian = True):
        if os.path.isfile(file_path):
            data = open(file_path, "rb").read()
            # read header
            header_struct = self.__get_header_struct(is_little_endian)
            start = 0
            size_of_header = struct.calcsize(header_struct);
            header = struct.unpack(header_struct, data[start: start + size_of_header])
            start += size_of_header
            size = header[0]
            assert size >= 0

            # read contents
            body_struct = self.__get_body_struct(is_little_endian)
            self.resize(size)
            size_of_double = struct.calcsize(body_struct);
            for i in range(len(self)):
                value = struct.unpack(body_struct, data[start: start + size_of_double])
                self[i] = float(value[0])
                start += size_of_double
        else:
            self._logger.error("file not found: %s" % (file_path))

    def save(self, file_path, is_little_endian = True):
        size = self._num_of_size

        fout = open(file_path, "wb")
        # write header
        header_struct = self.__get_header_struct(is_little_endian)
        header = struct.pack(header_struct, size)
        fout.write(header)
        
        # write elements
        body_struct = self.__get_header_struct(is_little_endian)
        for i in range(size):
            value_str = struct.pack(body_struct, self[i])
            fout.write(value_str)
        
        fout.close()

                


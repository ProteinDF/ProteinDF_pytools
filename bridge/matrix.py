#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import struct
import copy
import math
import numpy

from vector import Vector

"""
Matrix (for general matrix) and SymmetricMatrix (for symmetric matrix) class
using numpy.ndarray
"""

class Matrix(object):
    """
    >>> a = Matrix()
    >>> a.rows
    0
    >>> a.cols
    0
    >>> B = Matrix(3, 5)
    >>> B.rows
    3
    >>> B.cols
    5
    >>> B.set(0, 1, 1.0)
    >>> math.fabs(B.get(0, 1) - 1.0) < 1.0E-5
    True
    >>> C = Matrix([[7, 4, -1], [3, 0, 5]])
    >>> math.fabs(C.get(0, 1) - 4.0) < 1.0E-5
    True
    >>> D = Matrix([[8, 4, 2], [1, 3, -6], [-7, 0, 5]])
    >>> CD = C * D
    >>> (CD == Matrix([[67, 40, -15], [-11, 12, 31]]))
    True
    """
    __header_struct_little_endian = "<iii"
    __header_struct_big_endian = ">iii"
    __body_struct_little_endian = "<d"
    __body_struct_big_endian = ">d"
    _type = 'GE'
    
    def __init__(self, *args, **kwargs):
        self._rows = 0
        self._cols = 0
        self._type = 'GE'
        self._data = None

        size_of_args = len(args)
        if size_of_args == 1:
            if args[0] is None:
                return
            elif isinstance(args[0], list):
                self._data = numpy.array(args[0], numpy.float)
                assert(self._data.ndim == 2)
                self._rows, self._cols = self._data.shape
                return
            else:
                raise
        elif size_of_args == 2:
            if isinstance(args[0], int) and isinstance(args[1], int):
                self._rows = args[0]
                self._cols = args[1]
                self._data = numpy.array(
                    [[0.0 for c in range(self.cols)] for r in range(self.rows)],
                    numpy.float)
                return
            else:
                raise
            
        if kwargs:
            self._rows = kwargs.get('row', 0)
            self._cols = kwargs.get('col', 0)
            matrix_type = kwargs.get('type', None)
            if (matrix_type == 'GE'):
                data = kwargs.get('data', None)
                if data:
                    self._data = numpy.array(
                        [ [0.0 for c in range(self.cols)] for r in range(self.rows)],
                        numpy.float)
                    index = 0
                    for r in range(self.rows):
                        for c in range(self.cols):
                            self.set(r , c, buf[index])
                            index += 1
                    return
                else:
                    raise
            else:
                raise
            
    def copy(self):
        answer = copy.deepcopy(self)
        return answer

    def get_symmetric_matrix(self):
        answer = None
        dim = self.rows
        if (dim == self.col):
            answer = SymmetricMatrix(dim)
            for r in range(dim):
                for c in range(r):
                    v1 = self.get(r, c)
                    v2 = self.get(c, r)
                    if (math.fabs(v1 - v2) > 1.0E-5):
                        print("warning: %f(%d, %d) != %f(%d, %d)"
                              % (v1, r, c, v2, c, r))
                    answer.set(r, c, v1)
                answer.set(r, r, self.get(r, r))
        else:
            raise

        return answer
    
    def clear(self):
        self._data = None

    def resize(self, new_rows, new_cols):
        new_data = numpy.array([ [0.0 for c in range(new_cols)] for r in range(new_rows) ])
        for r in range(min(self.rows, new_rows)):
            for c in range(min(self.cols, new_cols)):
                new_data[r, c] = self._data[r, c]
        self._rows = new_rows
        self._cols = new_cols
        self._data = new_data

    # --------------------------------------------------------------------------
    @property
    def rows(self):
        return self._rows

    @property
    def cols(self):
        return self._cols

    @property
    def type(self):
        return self._type
    
    # --------------------------------------------------------------------------
    def get(self, row, col):
        assert((0 <= row) and (row < self.rows))
        assert((0 <= col) and (col < self.cols))
        return self._data[row, col]
    
    def set(self, row, col, value):
        assert((0 <= row) and (row < self.rows))
        assert((0 <= col) and (col < self.cols))
        self._data[row, col] = value

    def add(self, row, col, value):
        assert((0 <= row) and (row < self.rows))
        assert((0 <= col) and (col < self.cols))
        self._data[row, col] += value

    def transpose(self):
        numpy.transpose(self._data)
        return self
        
    def select(self, start_row, start_col, end_row, end_col):
        """
        select matrix sub-block
        """
        assert((0 <= start_row) and (start_row < self.rows))
        assert((0 <= start_col) and (start_col < self.cols))
        assert((0 < end_row) and (end_row <= self.rows))
        assert((0 < end_col) and (end_col <= self.cols))
        new_row_size = end_row - start_row
        new_col_size = end_col - start_col
        assert(new_row_size > 0)
        assert(new_col_size > 0)
        answer = Matrix(new_row_size, new_col_size)
        for r in range(start_row, end_row):
            for c in range(start_col, end_col):
                answer.set(r - start_row, c - start_col, self.get(r, c))
        return answer

    def get_row_vector(self, row):
        assert(0 <= row)
        assert(row < self.rows)

        cols = self.col
        v = Vector(cols)
        for i in range(cols):
            v[i] = self.get(row, i)
        return v
        
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
        answer = ''
        for order in range(0, self.col, 10):
            answer += '       '
            for j in range(order, min(order +10, self.col)):
                answer += '   %5d th' % (j +1)
            answer += '\n  ---- '

            for j in range(order, min(order +10, self.col)):
                answer += '-----------'
            answer += '\n'

            for i in range(0, self.get_num_of_rows()):
                answer += ' %5d ' % (i +1)
                for j in range(order, min(order +10, self.col)):
                    answer += ' % 10.6f' % (self.get(i, j))
                answer += '\n'
            answer += '\n\n'
        return answer

    def get_raw_data(self):
        raw = {}
        raw['row'] = self.rows
        raw['col'] = self.cols
        raw['type'] = self.type

        # setup data
        data = [ 0.0 for x in range(self.rows * self.col) ]
        index = 0
        for r in range(self.rows):
            for c in range(self.col):
                data[index] = self.get(r, c)
                index += 1
        raw['data'] = data
        
        return raw

    def get_buffer(self):
        return buffer(self._data.tostring())
    
    def __add__(self, other):
        assert isinstance(other, Matrix)
        assert (self.rows == other.rows)
        assert (self.col == other.cols)

        answer = self.copy()
        answer += other
        return answer

    def __iadd__(self, other):
        assert isinstance(other, Matrix)
        assert (self.rows == other.rows)
        assert (self.cols == other.cols)
            
        self.__data += other.__data
        return self

    def __sub__(self, other):
        assert isinstance(other, Matrix)
        assert (self.rows == other.rows)
        assert (self.cols == other.cols)

        answer = self.copy()
        answer -= other
        return answer

    def __isub__(self, other):
        assert isinstance(other, Matrix)
        assert (self.rows == other.rows)
        assert (self.cols == other.cols)
            
        self.__data -= other.__data
        return self

    def __mul__(self, other):
        if isinstance(other, Matrix):
            # matrix * matrix
            assert(self.col == other.row)
            
            A = numpy.matrix(self._data)
            B = numpy.matrix(other._data)
            C = A * B
            
            answer = Matrix()
            answer._data = C.getA()
            return answer
        elif isinstance(other, Vector):
            # matrix * (coulmn)vector
            assert(self.col == other.size())

            A = numpy.matrix(self._data)
            B = numpy.matrix([other._data])
            B = B.getT()
            C = A * B

            # TODO: to be simply!
            C = C.getT()
            a = C.tolist()
            answer = Vector(a[0])
            return answer

    def __eq__(self, other):
        assert(isinstance(other, Matrix))
        answer = False
        if ((self.rows == other.rows) and
            (self.cols == other.cols)):
            answer = True
            for r in range(self.rows):
                for c in range(self.cols):
                    if math.fabs(self.get(r, c) - other.get(r, c)) > 1.0E-5:
                        answer = False
                        break
        return answer

    def __ne__(self, other):
        return not self.__eq__(other)
    
########################################################################
# 

class SymmetricMatrix(Matrix):
    """
    >>> A = SymmetricMatrix()
    >>> A.rows
    0
    >>> A.cols
    0
    >>> B = SymmetricMatrix(5)
    >>> B.rows
    5
    >>> B.cols
    5
    >>> B.set(0, 1, 1.0)
    >>> B.get(0, 1)
    1.0
    >>> B.get(1, 0)
    1.0
    """
    __header_struct_little_endian = "<iii" 
    __header_struct_big_endian = ">iii"
    __body_struct_little_endian = "<d" 
    __body_struct_big_endian = ">d"
    
    def __init__(self, *args, **kwargs):
        Matrix.__init__(self, None)
        self._type = 'SY'

        size_of_args = len(args)
        if size_of_args == 1:
            if isinstance(args[0], int):
                dim = args[0]
                self._rows = dim
                self._cols = dim
                self._data = numpy.array(
                    [[0.0 for c in range(dim)] for r in range(dim)],
                    numpy.float)
                return
            elif isinstance(args[0], list):
                self._data = numpy.array(args[0], numpy.float)
                assert(self._data.ndim == 2)
                self._rows, self._cols = self._data.shape
                assert(self.rows == self.cols)
                return
            else:
                raise
                
        if kwargs:
            self._rows = kwargs.get('row', 0)
            self._cols = kwargs.get('col', 0)
            assert(self.rows == self.cols)
            matrix_type = kwargs.get('type', None)
            if matrix_type == 'SP':
                data = kwargs.get('data', None)
                if data:
                    self._data = numpy.array(
                        [ [0.0 for c in range(self.cols)] for r in range(self.rows)],
                        numpy.float)
                    index = 0
                    for r in range(self.rows):
                        for c in range(r +1):
                            self.set(r , c, data[index])
                            index += 1
                    return
                else:
                    raise
            elif matrix_type == 'SY':
                data = kwargs.get('data', None)
                if data:
                    self._data = numpy.array(
                        [ [0.0 for c in range(self.cols)] for r in range(self.rows)],
                        numpy.float)
                    index = 0
                    for r in range(self.rows):
                        for c in range(self.col):
                            if r >= c:
                                self.set(r , c, data[index])
                            index += 1
                    return
                else:
                    raise
            else:
                raise

    @property
    def dim(self):
        assert(self.rows == self.cols)
        return self.rows
            
    def get_general_matrix(self):
        answer = Matrix(self.rows, self.cols)
        for r in range(self.rows):
            for c in range(r):
                v = self.get(r, c)
                answer.set(r, c, v)
                answer.set(c, r, v)
            answer.set(r, r, self.get(r, r))
        return answer

    def resize(self, new_dim):
        new_data = numpy.array([ [0.0 for c in range(new_dim)] for r in range(new_dim) ])
        for r in range(min(self.rows, new_dim)):
            for c in range(r +1):
                new_matrix[r, c] = self._matrix[r, c]
        self._rows = new_dim
        self._cols = new_dim
        self._data = new_data

    def get(self, row, col):
        if (row < col):
            row, col = col, row
        return Matrix.get(self, row, col)


    def set(self, row, col, value):
        if (row < col):
            row, col = col, row
        Matrix.set(self, row, col, value)

    def eig(self):
        """
        return the eigenvalues and eigenvectors.
        """
        w = Vector()
        v = Matrix()
        if (self.dim > 1):
            w._data, v._data = numpy.linalg.eigh(self._data, 'L')
            v._data = numpy.transpose(v._data) # to treat column vector as eigenvectors
        return w, v
    
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
        answer = ''
        dim = self.dim
        for order in range(0, dim, 10):
            answer += '       '
            for j in range(order, min(order +10, dim)):
                answer += '  %6d th' % (j +1)
            answer += '\n ======'

            for j in range(order, min(order +10, dim)):
                answer += '==========='
            answer += '\n'

            for i in range(0, dim):
                answer += ' %6d ' % (i +1)
                for j in range(order, min(order +10, dim)):
                    if (j > i):
                        answer += '   ------- '
                    else:
                        answer += ' % 10.6f' % (self.get(i, j))
                answer += '\n'
            answer += '\n\n'
        return answer

    def get_raw_data(self):
        dim = self.dim
        
        raw = {}
        raw['row'] = dim
        raw['col'] = dim
        raw['type'] = 'SP'

        # setup data
        data = [ 0.0 for x in range(dim * (dim +1) / 2) ]
        index = 0
        for r in range(dim):
            for c in range(r +1):
                data[index] = self.get(r, c)
                index += 1
        raw['data'] = data

        return raw

    def __mul__(self, other):
        A = self.get_general_matrix()
        B = copy.deepcopy(other)
        if isinstance(other, SymmericMatrix):
            B = other.get_general_matrix()
        return A * B

    def __eq__(self, other):
        assert(isinstance(other, Matrix))
        answer = False
        if ((self.rows == other.rows) and
            (self.cols == other.cols)):
            answer = True
            for r in range(self.rows):
                for c in range(r +1):
                    if math.fabs(self.get(r, c) - other.get(r, c)) > 1.0E-5:
                        answer = False
                        break
        return answer


if __name__ == '__main__':
    import doctest
    doctest.testmod()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import optparse
import math
import copy

from vector import Vector

class Position(object):
    """
    >>> p = Position([5.0, 3.0, -1.2])
    >>> p.x
    5.0
    >>> p.y
    3.0
    >>> p.z
    -1.2

    >>> p.x = 0.0
    >>> p.y = 1.0
    >>> p.z = 2.0
    >>> p == Position([0, 1, 2])
    True

    >>> abs(abs(p) - 2.23606) < 1.0E-5
    True

    >>> p.norm()
    >>> p == Position([0, 1/math.sqrt(5), 2/math.sqrt(5)])
    True
    
    >>> p.move_to([3, 4, 5])
    >>> p == Position([3.0, 4.0, 5.0])
    True

    >>> tmp = p * 2.0
    >>> tmp == Position([6.0, 8.0, 10.0])
    True
    
    >>> tmp = (-1.0 * p)
    >>> tmp == Position([-3.0, -4.0, -5.0])
    True
    
    >>> a = Position([1, 2, 3])
    >>> b = Position([2, 3, 4])
    >>> a * b
    20.0

    >>> p = a + b
    >>> p == Position([3, 5, 7])
    True

    >>> a += b
    >>> a == Position([3, 5, 7])
    True

    >>> p = a - b
    >>> p == Position([1, 2, 3])
    True

    >>> a -= b
    >>> a == Position([1, 2, 3])
    True
    
    """
    _position = [0.0, 0.0, 0.0]
    
    def __init__(self, rhs = [0.0, 0.0, 0.0]):
        self.epsilon = 1.0E-5
        self._position = [0.0, 0.0, 0.0]
        if ((isinstance(rhs, (list, tuple)) == True) and (len(rhs) == 3)):
            self._position[0] = float(rhs[0])
            self._position[1] = float(rhs[1])
            self._position[2] = float(rhs[2])
        elif (isinstance(rhs, Position) == True):
            self._position = rhs._position
        else:
            raise InputError, "illegal input"

    # --------------------------------------------------------------------------
    @property
    def xyz(self):
        return self._position

    def __get_x(self):
        return self._position[0]

    def __set_x(self, new_x):
        self._position[0] = float(new_x)

    x = property(__get_x, __set_x)

    def __get_y(self):
        return self._position[1]

    def __set_y(self, new_y):
        self._position[1] = float(new_y)

    y = property(__get_y, __set_y)

    def __get_z(self):
        return self._position[2]

    def __set_z(self, new_z):
        self._position[2] = float(new_z)

    z = property(__get_z, __set_z)

    # --------------------------------------------------------------------------
    def get_dict_data(self):
        return self._position

    def move_to(self, position):
        tmp = Position(position)
        self._position = tmp._position

    def square_distance_from(self, other):
        d2 = 0.0
        for i in range(3):
            tmp = self._position[i] - other._position[i]
            d2 += tmp * tmp
        return d2
            
    def distanceFrom(self, other):
        d2 = self.square_distance_from(other)
        return math.sqrt(d2)

    def norm(self):
        n = self.__abs__()
        self._position = [x / n for x in self._position]

    def rotate(self, mat):
        v1 = Vector(self._position)
        v2 = mat * v1
        self._position = v2.to_list()
        
    def __str__(self):
        answer = "(% 10.6f, % 10.6f, % 10.6f)" % (self.x,
                                                  self.y,
                                                  self.z)
        return answer

    def __eq__(self, rhs):
        answer = False
        if (isinstance(rhs, Position) == True):
            if (self.distanceFrom(rhs) < self.epsilon):
                answer = True
        else:
            raise
        return answer

    def __ne__(self, rhs):
        return not(self.__eq__(rhs))

    def __neg__(self):
        return Position([-x for x in self._position])

    def __abs__(self):
        return math.sqrt(sum([x * x for x in self._position]))

    def __add__(self, rhs):
        if (isinstance(rhs, Position) == True):
            return Position([ x + y for x, y in zip(self._position, rhs._position)])
        else:
            raise InputError, "illegal input: Position is required."

    __iadd__ = __add__

    def __sub__(self, rhs):
        return self.__add__(-rhs)

    __isub__ = __sub__

    def __sub__(self, rhs):
        return self.__add__(-rhs)

    def __mul__(self, rhs2):
        rhs1 = copy.deepcopy(self)
        if (isinstance(rhs2, Position) == True):
            return sum([x * y for x, y in zip(rhs1._position, rhs2._position)])
        elif (isinstance(rhs2, float) == True):
            rhs1._position = [x * rhs2 for x in rhs1._position]
            return rhs1
        else:
            raise

    __rmul__ = __mul__

    def __imul__(self, rhs):
        v = float(rhs)
        self._position[0] *= v
        self._position[1] *= v
        self._position[2] *= v
        return self
        
    def __div__(self, rhs):
        v = float(rhs)
        return Position([x / v for x in self._position])
        
    def __idiv__(self, rhs):
        return self.__imul__(1.0 / float(rhs))
        
if __name__ == "__main__":
    import doctest
    doctest.testmod()

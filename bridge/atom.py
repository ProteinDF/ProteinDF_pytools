#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
import math

from position import Position
from periodictable import PeriodicTable

class Atom(object):
    """
    >>> a = Atom()
    >>> a.symbol = 'Fe'
    >>> a.atomic_number
    26
    >>> a.symbol
    'Fe'
    >>> a.charge = -0.2
    >>> math.fabs(a.charge - (-0.2)) < 1.0E-10
    True
    >>> b = Atom(a)
    >>> b.symbol = 'Na'
    >>> a.symbol
    'Fe'
    >>> b.symbol
    'Na'
    """
    _xyz = Position()
    _atomic_number = 0
    _name = ""
    _charge = 0.0
    _path = ""
    
    def __init__(self, rhs = None,
                 symbol = 'X', position = Position(),
                 charge = 0.0, name = "", path = ""):
        self._xyz = position
        self._atomic_number = PeriodicTable.get_atomic_number(symbol)
        self._name = name
        self._charge = charge
        self._path = path
        
        if (isinstance(rhs, Atom) == True):
            self._xyz = copy.copy(rhs._xyz)
            self._atomic_number = copy.copy(rhs._atomic_number)
            self._name = copy.copy(rhs._name)
            self._charge = copy.copy(rhs._charge)
            self._path = copy.copy(rhs._path)
        elif (isinstance(rhs, dict) == True):
            self.set_by_dict_data(rhs)

    # move ---------------------------------------------------------------------
    def move_to(self, position):
        self.position.moveTo(position)
        return self

    def shift_by(self, direction):
        self.position += direction
        return self

    # --------------------------------------------------------------------------
    def __get_xyz(self):
        return self._xyz

    def __set_xyz(self, p):
        assert(isinstance(p, Position))
        self._xyz = p

    xyz = property(__get_xyz, __set_xyz)
        
    # --------------------------------------------------------------------------
    def __get_atomic_number(self):
        return self._atomic_number

    def __set_atomic_number(self, an):
        self._atomic_number = an
        
    atomic_number = property(__get_atomic_number, __set_atomic_number)
        
    # --------------------------------------------------------------------------
    def __get_symbol(self):
        answer = PeriodicTable.get_symbol(self.atomic_number)
        return answer

    def __set_symbol(self, symbol):
        self._atomic_number = PeriodicTable.get_atomic_number(symbol)

    symbol = property(__get_symbol, __set_symbol)

    # --------------------------------------------------------------------------
    def __get_name(self):
        return self._name

    def __set_name(self, name):
        self._name = str(name)

    name = property(__get_name, __set_name)

    # --------------------------------------------------------------------------
    def __get_charge(self):
        return self._charge

    def __set_charge(self, charge):
        self._charge = float(charge)

    charge = property(__get_charge, __set_charge)
    # --------------------------------------------------------------------------
    def __get_path(self):
        return self._path

    def __set_path(self, path):
        self._path = path

    path = property(__get_path, __set_path)
    # --------------------------------------------------------------------------
    def set_by_dict_data(self, data):
        self._atomic_number = data.get('Z', 0)
        self._name = data.get('name', '')
        self._charge = data.get('Q', 0.0)
        self._xyz = Position(data.get('xyz'))
        return self

    def get_dict_data(self):
        data = {}
        data['Z'] = self._atomic_number
        data['name'] = self._name
        data['Q'] = self._charge
        data['xyz'] = self._xyz.get_dict_data()

        return data

    # --------------------------------------------------------------------------
    def __str__(self):
        answer = "%2s(%4s) (% 8.3f, % 8.3f, % 8.3f) Z=% .2f" % (
            self.symbol,
            self.name,
            self.position.x,
            self.position.y,
            self.position.z,
            self.charge)
        return answer
    
    def __eq__(self, rhs):
        answer = False
        if ((isinstance(rhs, Atom) == True) and
            (self.atomic_number == rhs.atomic_number) and
            (math.fabs(self.charge - rhs.charge) < 1.0E-10) and
            (self.position == rhs.position)):
            answer = True
        return answer

    def __ne__(self, rhs):
        return not(self.__eq__(rhs))

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()

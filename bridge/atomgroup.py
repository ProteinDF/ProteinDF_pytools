#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import copy

from utils import *
from atom import Atom
from select import Select

class AtomGroup(object):
    """
    >>> group1 = AtomGroup()
    >>> atom1 = Atom(symbol='C')
    >>> atom2 = Atom(symbol='H')
    >>> atom3 = Atom(symbol='N')
    >>> subgrp = AtomGroup()
    >>> subgrp.set_atom('C1', atom1)
    >>> subgrp.set_atom('H1', atom2)
    >>> subgrp.set_atom('N1', atom3)
    >>> group1.set_group('grp', subgrp)
    >>> group1['grp']['C1'].symbol
    'C'
    >>> group1.get_number_of_atoms()
    0
    >>> group1.get_number_of_groups()
    1
    >>> group1.get_number_of_all_atoms()
    3
    >>> group1.sum_of_atomic_number()
    14.0
    """
    _atoms = {}
    _groups = {}
    _name = ""
    _charge = 0.0
    _path = "/"

    def __init__(self, rhs = None):
        self.initialize()
        
        if (isinstance(rhs, AtomGroup) == True):
            self._atoms = copy.deepcopy(rhs._atoms)
            self._groups = copy.deepcopy(rhs._groups)
            self._name = copy.copy(rhs._name)
            self._charge = copy.copy(rhs._charge)
            self._path = copy.copy(rhs._path)
        elif (isinstance(rhs, dict) == True):
            self.set_by_dict_data(rhs)

    def initialize(self):
        self._atoms = {}
        self._groups = {}
        self._name = ""
        self._charge = 0.0
        self._path = "/"
            
    # --------------------------------------------------------------------------
    def get_number_of_groups(self):
        return len(self._groups)

    def get_number_of_atoms(self):
        return len(self._atoms)

    def get_number_of_all_atoms(self):
        answer = 0
        for key, grp in self._groups.items():
            answer += grp.get_number_of_all_atoms()
        answer += len(self._atoms)
        return answer

    # move ---------------------------------------------------------------------
    def shift_by(self, direction):
        for key, grp in self._groups.items():
            grp.shift_by(direction)
        for key, atm in self._atoms.items():
            atm.shift_by(direction)

    # --------------------------------------------------------------------------
    def sum_of_atomic_number(self):
        """
        原子数の総和を返す
        """
        answer = 0.0
        for key, grp in self._groups.items():
            answer += grp.sum_of_atomic_number()
        for key, atm in self._atoms.items():
            answer += atm.atomic_number
        return answer
            
    # --------------------------------------------------------------------------
    def groups(self):
        """
        原子団のリストを返す
        """
        keys = self._groups.keys()
        sort_nicely(keys)
        for k in keys:
            yield(k, self._groups[k])

    def get_group(self, key):
        key = str(key)
        return self._groups.get(key, None)

    def set_group(self, key, value):
        key = str(key)
        assert(isinstance(value, AtomGroup))
        self._groups[key] = copy.deepcopy(value)

    def has_group(self, key):
        key = str(key)
        return self._groups.has_key(key)

    def erase_group(self, key):
        assert(isinstance(key, str))
        self._groups.pop(key, None)

    #def get_group_list(self):
    #    return self.data['groups'].keys()

    # --------------------------------------------------------------------------
    def atoms(self):
        """
        原子のリストを返す
        """
        keys = self._atoms.keys()
        sort_nicely(keys)
        for k in keys:
            yield(k, self._atoms[k])

    def get_atom(self, key):
        key = str(key)
        return self._atoms.get(key, None)

    def set_atom(self, key, value):
        key = str(key)
        assert(isinstance(value, Atom))
        self._atoms[key] = copy.deepcopy(value)

    def has_atom(self, key):
        key = str(key)
        return self._atoms.has_key(key)

    def erase_atom(self, key):
        assert(isinstance(key, str))
        self._atoms.pop(key, None)

    #def get_atom_list(self):
    #    return self.data['atoms'].keys()

    # --------------------------------------------------------------------------
    def __get_name(self):
        return self._name

    def __set_name(self, name):
        assert(isinstance(name, str))
        self._name = name

    name = property(__get_name, __set_name)
    # --------------------------------------------------------------------------
    def __get_charge(self):
        return self._charge

    def __set_charge(self, value):
        assert(isinstance(value, float))
        self._charge = value

    charge = property(__get_charge, __set_charge)
    # --------------------------------------------------------------------------
    def __get_path(self):
        return self._path

    def __set_path(self, value):
        if (self._path != value):
            self._path = value
            if (self._path[-1] != '/'):
                self._path.append('/')
            self._update_path()

    path = property(__get_path, __set_path)
    # --------------------------------------------------------------------------
    def merge(self, rhs):
        """原子団を結合する
        """
        assert(isinstance(rhs, AtomGroup) == True)
        for key, group in rhs._groups.items():
            self.merge_group(key, group)
        for key, atom in rhs._atoms.items():
            self.set_atom(key, copy.deepcopy(atom))

    # --------------------------------------------------------------------------
    def select(self, selecter):
        assert(isinstance(selecter, Select) == True)
        self._update_path()

        answer = None
        if (selecter.is_match(self) == True):
            answer = copy.deepcopy(self)
        else:
            answer = AtomGroup()
            answer.name = self.name
            for key, group in self.groups():
                tmp = group.select(selecter)
                if ((tmp.get_number_of_groups() != 0) or
                    (tmp.get_number_of_atoms() != 0)):
                    answer.set_group(key, tmp)
            for key, atom in self.atoms():
                if (selecter.is_match(atom) == True):
                    answer.set_atom(key, atom)
            answer.path = self.path
        return answer

    # private method -----------------------------------------------------------
    def _merge_group(self, key, group):
        assert(isinstance(key, str) == True)
        assert(isinstance(group, AtomGroup) == True)
        if (self.has_group(key) == True):
            self._groups[key].merge(group)
        else:
            self.set_group(key, group)

    def _update_path(self):
        for key, group in self._groups.iteritems():
            group.path = "%s%s/" % (self._path, key)
        for key, atom in self._atoms.items():
            atom.path = "%s%s" % (self._path, key)
    
    # --------------------------------------------------------------------------
    def __iand__(self, rhs):
        """
        implement of '&=' operator
        """
        assert(isinstance(rhs, AtomGroup) == True)
        #self.update_path(self.get_path())
        #rhs.update_path(rhs.get_path())

        for key, group in rhs._groups:
            if (self.has_group(key) != True):
                self.erase_group(key)
            else:
                self._groups[key] &= rhs._groups[key]
                if ((self._groups[key].get_number_of_groups() == 0) and
                    (self._groups[key].get_number_of_atoms() == 0)):
                    self.erase_group(key)

        for key, atom in rhs._atoms:
            if (self.has_atom(key) != True):
                self.erase_atom(key)

        return self

    def __ior__(self, rhs):
        """
        implement of '|=' operator
        """
        assert(isinstance(rhs, AtomGroup) == True)
        #self.update_path(self.get_path())
        #rhs.update_path(rhs.get_path())

        self.merge(rhs)
        return self

    def __ixor__(self, rhs):
        """
        implement of '^=' operator
        """
        assert(isinstance(rhs, AtomGroup) == True)
        #self.update_path(self.get_path())
        #rhs.update_path(rhs.get_path())

        for key, group in rhs._groups:
            if (self.has_group(key) == True):
                self._groups.__ixor__(group)
                if ((self._groups[key].get_number_of_groups() == 0) and
                    (self._groups[key].get_number_of_atoms() == 0)):
                    self.erase_group(key)
            else:
                self.set_group(key, group)

        for key, atom in rhs._atoms.items():
            if (self.has_atom(key) == True):
                self.erase_atom(key)
            else:
                self.set_atom(key, atom)

        return self

    # --------------------------------------------------------------------------
    def set_by_dict_data(self, data):
        assert(isinstance(data, dict) == True)

        if (data.has_key('groups') == True):
            for key in data['groups'].keys():
                #print("key=", key)
                atomgroup = AtomGroup(data['groups'][key])
                #print("atomgroup=", atomgroup)
                self.set_group(key, atomgroup)
        if (data.has_key('atoms') == True):
            for key in data['atoms'].keys():
                #print("key(atom)=", key)
                atom = Atom(data['atoms'][key])
                #print("atom=", atom)
                self.set_atom(key, atom)
        self.name = data.get('name', '')
        self.charge = data.get('charge', 0.0)

        self._update_path()
        return self

    def get_dict_data(self):
        self._update_path()
        data = {}
        if (len(self._groups) > 0):
            data.setdefault('groups', {})
            for key in self._groups.keys():
                data['groups'][key] = self._groups[key].get_dict_data()
        if (len(self._atoms) > 0):
            data.setdefault('atoms', {})
            for key in self._atoms.keys():
                data['atoms'][key] = self._atoms[key].get_dict_data()
        data['name'] = self._name
        data['charge'] = self._charge
        return data

    def __str__(self):
        self._update_path()
        
        answer = ""
        for key, atomgroup in self.groups():
            #answer += "<atom_group(%s):path=%s>\n" % (self.name, self.path)
            answer += str(atomgroup)
        for key, atom in self.atoms():
            answer += "%s%s: %s\n" % (self.path, key,
                                      str(atom))
        return answer

    def __getitem__(self, key):
        """operator[] for getter"""
        key = str(key)
        if (self.has_group(key) == True):
            return self._groups[key]
        elif (self.has_atom(key) == True):
            return self._atoms[key]
        else:
            raise KeyError, key

    def __setitem__(self, key, value):
        """operator[] for setter"""
        assert(isinstance(key, str))
        if (isinstance(value, AtomGroup) == True):
            self._groups.__setitem__(key, value)
        elif (isinstance(value, Atom) == True):
            self._atoms.__setitem__(key, value)
        else:
            raise ValueError

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()

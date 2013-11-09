#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re
import array
import copy
import math
import pprint

class PrimitiveGTO(object):
    """
    >>> pgto = PrimitiveGTO(2.80806400E+03, 2.01783000E-03)
    >>> math.fabs(pgto.exp - 2.80806400E+03) < 1.0E-15
    True
    >>> math.fabs(pgto.coef - 2.01783000E-03) < 1.0E-15
    True
    """
    def __init__(self, *args, **kwargs):
        self.exp = kwargs.get('exp', 0.0)
        self.coef = kwargs.get('coef', 1.0)
        if (len(args) == 2):
            self.exp = args[0]
            self.coef = args[1]
        
    def __get_exp(self):
        return self._exp

    def __set_exp(self, e):
        self._exp = float(e)

    exp = property(__get_exp, __set_exp)

    def __get_coef(self):
        return self._coef

    def __set_coef(self, coef):
        self._coef = float(coef)

    coef = property(__get_coef, __set_coef)

    def __str__(self):
        output = "    {0: e} {1: e}\n".format(self.exp, self.coef)
        return output

    #def __getstate__(self):
    #    odict = {}
    #    odict['exp'] = self.exp
    #    odict['coef'] = self.coef
    #    return odict
    
    #def _init_by_dict(self, odict):
    #    assert(isinstance(odict, dict) == True)
    #    self.exp  = odict.get('exp', 0.0)
    #    self.coef = odict.get('coef', 1.0)
    
        
class ContractedGTO(list):
    """
    >>> cgto = ContractedGTO('p', 3)
    >>> cgto[0] = PrimitiveGTO(2.80806400E+03, 2.01783000E-03)
    >>> cgto[1] = PrimitiveGTO(4.21138300E+02, 1.54332000E-02)
    >>> cgto[2] = PrimitiveGTO(9.55866200E+01, 7.55815500E-02)
    >>> cgto.shell_type
    'p'
    >>> len(cgto)
    3
    """
    _shell_types = ['s', 'p', 'd', 'spd']
    _orb_types = ['s', 'px', 'py', 'pz', 'dxy', 'dyz', 'dzx', 'dxx-yy', 'dzz']

    def __init__(self, *args, **kwargs):
        shell_type = 's'
        size = 0
        
        if len(args) == 2:
            shell_type = args[0]
            size = int(args[1])
        if kwargs.has_key('shell_type'):
            shell_type = kwargs.get('shell_type')
        if kwargs.has_key('size'):
            size = kwargs.get('size')

        self.shell_type = shell_type
        list.__init__(self, [PrimitiveGTO() for x in range(size)])

        if kwargs.has_key('pGTOs'):
            pgtos = kwargs.get('pGTOs', [])
            self.__init__(size = len(pgtos))
            for i in range(len(pgtos)):
                self[i] = PrimitiveGTO(**(pgtos[i]))
            self.shell_type = kwargs.get('shell_type', 's')
            self.scale_factor = kwargs.get('scale_factor', 1.0)
    
    # shell_type ---------------------------------------------------------------
    def __get_shell_type_id(self):
        return self._shell_type_id

    shell_type_id = property(__get_shell_type_id)

    def __get_shell_type(self):
        if not '_shell_type' in self.__dict__:
            self._shell_type = 's'
        return self._shell_types[self._shell_type_id]

    def __set_shell_type(self, shell_type):
        self._shell_type_id = None
        for index, st in enumerate(self._shell_types):
            if st == shell_type:
                self._shell_type_id = index
                break

    shell_type = property(__get_shell_type, __set_shell_type)

    # scale factor -------------------------------------------------------------
    def __get_scale_factor(self):
        if not '_scale_factor' in self.__dict__:
            self._scale_factor = 1.0
        return self._scale_factor

    def __set_scale_factor(self, value):
        self._scale_factor = value

    scale_factor = property(__get_scale_factor, __set_scale_factor)

    @classmethod
    def get_supported_shell_types(cls):
        return cls._shell_types

    @classmethod
    def get_shell_type_id(cls, shell_type):
        """
        shell_typeに対応するidを返す

        s: 0, p: 1, d: 2
        """
        answer = -1
        for i, st in enumerate(cls._shell_types):
            if shell_type == st:
                answer = i
                break
        return answer

    @classmethod
    def get_basis_type(cls, shell_type_id, basis_id):
        """
        軌道種類の文字列を返す
        shell_type_id=0, basis_id=0: s
        shell_type_id=1, basis_id=0: px
        shell_type_id=1, basis_id=1: py
        shell_type_id=1, basis_id=2: pz
        shell_type_id=2, basis_id=0: dxy
        shell_type_id=2, basis_id=1: dyz
        shell_type_id=2, basis_id=2: dzx
        shell_type_id=2, basis_id=3: dxx-yy
        shell_type_id=2, basis_id=4: dzz
        """
        tbl = [0, 1, 4]
        index = tbl[shell_type_id] + basis_id
        return cls._orb_types[index]
    
    def expand(self):
        """
        shell_typeが'spd'などの場合に分解する
        """
        answer = None
        if (self.shell_type == 'spd'):
            cgto_s = copy.deepcopy(self)
            cgto_s.shell_type = 's'
            cgto_p = copy.deepcopy(self)
            cgto_p.shell_type = 'p'
            cgto_d = copy.deepcopy(self)
            cgto_d.shell_type = 'd'
            answer = [cgto_s, cgto_p, cgto_d]
        else:
            answer = [self]
        return answer
    
    def __str__(self):
        output = ""
        #output += " %s %d\n".format(self.shell_type, len(self))
        output += " %d\n" % (len(self))
        for pgto in self:
            output += str(pgto)
        return output

    #def __getstate__(self):
    #    odict = {}
    #    odict['shell_type'] = self.shell_type
    #    odict['scale_factor'] = self.scale_factor
    #    odict['pGTOs'] = []
    #    for pgto in self:
    #        odict['pGTOs'].append(pgto.__getstate__())
    #    assert(len(odict['pGTOs']) == len(self))
    #    return odict
        
    #def _init_by_dict(self, odict):
    #    assert(isinstance(odict, dict) == True)
    #    pgtos = odict.get('pGTOs', [])
    #    self.__init__(size = len(pgtos))
    #    for i in range(len(pgtos)):
    #        pgto = PrimitiveGTO()
    #        pgto.__setstate__(pgtos[i])
    #        self[i] = pgto
    #    self.shell_type = odict.get('shell_type', 's')
    #    self.scale_factor = odict.get('scale_factor', 1.0)
    
class BasisSet(list):
    """
    >>> bs = BasisSet('sample', 3)
    >>> bs.name
    'sample'
    >>> len(bs)
    3
    >>> bs[0] = ContractedGTO('p', 3)
    >>> bs[0][0] = PrimitiveGTO(2.80806400E+03, 2.01783000E-03)
    >>> bs[0][1] = PrimitiveGTO(4.21138300E+02, 1.54332000E-02)
    >>> bs[0][2] = PrimitiveGTO(9.55866200E+01, 7.55815500E-02)
    """
    def __init__(self, *args, **kwargs):
        self._name = ''
        size = 0
        if (len(args) > 0):
            self.name = args[0]
            if (len(args) > 1):
                size = args[1]
        list.__init__(self, [ContractedGTO() for x in range(size)])

        if kwargs.has_key('name'):
            self.name = kwargs.get('name')
        if kwargs.has_key('cGTOs'):
            cgtos = kwargs.get('cGTOs', [])
            size = len(cgtos)
            list.__init__(self, [ContractedGTO() for x in range(size)])
            for i in range(len(cgtos)):
                self[i] = ContractedGTO(**(cgtos[i]))
                
    def _get_name(self):
        return self._name

    def _set_name(self, name):
        self._name = name

    name = property(_get_name, _set_name)

    def get_num_of_CGTOs(self, shell_type):
        answer = 0
        for i in range(len(self)):
            if (self[i].shell_type == shell_type):
                answer += 1
        return answer

    def expand(self):
        """
        shell_typeが'spd'などの場合に分解する
        """
        tmp = BasisSet(self.name)
        for i in self:
            tmp.extend(i.expand())
        del self[:]
        return tmp

    def sort(self):
        CGTOs = {}
        for shell_type in ContractedGTO.get_supported_shell_types():
            CGTOs.setdefault(shell_type, [])
        for cgto in self:
            shell_type = cgto.shell_type
            CGTOs[shell_type].append(cgto)
        del self[:] # 全削除
        for shell_type in ContractedGTO.get_supported_shell_types():
            self.extend(CGTOs[shell_type])

    def get_raw_data(self):
        return self.__getstate__()

    #def setstate(self, odict):
    #    self.__setstate__(odict)
            
    def __str__(self):
        self.sort()
        
        output = ""
        output += "name:%s\n" % (self.name)

        for shell_type in ContractedGTO.get_supported_shell_types():
            num_of_CGTOs = self.get_num_of_CGTOs(shell_type)
            if ((shell_type == 'spd') and (num_of_CGTOs == 0)):
                continue
            output += " %d" % (num_of_CGTOs)
        output += "\n"
        
        for cgto in self:
            output += str(cgto)
        
        return output

    #def __getstate__(self):
    #    odict = {}
    #    odict['name'] = self.name
    #    odict['cGTOs'] = []
    #    for i in range(len(self)):
    #        odict['cGTOs'].append(self[i].__getstate__())
    #    assert(len(odict['cGTOs']) == len(self))
    #    return odict
    
    #def _init_by_dict(self, odict):
    #    assert(isinstance(odict, dict) == True)
    #    cgtos = odict.get('cGTOs', [])
    #    self.__init__(size = len(cgtos))
    #    for i in range(len(cgtos)):
    #        cgto = ContractedGTO()
    #        cgto.__setstate__(cgtos[i])
    #        self[i] = cgto
    #    self.name = odict.get('name', '')

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()

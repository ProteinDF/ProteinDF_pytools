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

import re
import array
import copy
import math

import pdfpytools as pdf

class PrimitiveGTO(object):
    """
    >>> pgto = PrimitiveGTO(2.80806400E+03, 2.01783000E-03)
    >>> math.fabs(pgto.exp - 2.80806400E+03) < 1.0E-15
    True
    >>> math.fabs(pgto.coef - 2.01783000E-03) < 1.0E-15
    True
    """
    def __init__(self, *args, **kwargs):
        if (len(args) == 1):
            if isinstance(args[0], PrimitiveGTO):
                rhs = args[0]
                self.exp = rhs.exp
                self.coef = rhs.coef
                return
            elif isinstance(args[0], dict):
                self.set_by_raw_data(args[0])
                return
            
        self.exp = kwargs.get('exp', 0.0)
        self.coef = kwargs.get('coef', 1.0)
        if (len(args) == 2):
            self.exp = args[0]
            self.coef = args[1]

    # exp --------------------------------------------------------------
    def __get_exp(self):
        return self._exp

    def __set_exp(self, e):
        self._exp = float(e)

    exp = property(__get_exp, __set_exp)

    # coef -------------------------------------------------------------
    def __get_coef(self):
        return self._coef

    def __set_coef(self, coef):
        self._coef = float(coef)

    coef = property(__get_coef, __set_coef)

    # normalize --------------------------------------------------------
    def normalize(self, shell_type):
        shell_type = shell_type.upper()
        max_angular = 0
        l = m = n = 0
        
        if shell_type == 'S':
            l = m = n = 0
            max_angular = 0
        elif shell_type == 'P':
            l = 1
            m = n = 0
            max_angular = 1
        elif shell_type == 'D':
            l = m = 1
            n = 0
            max_angular = 2
        elif shell_type == 'F':
            l = m = n = 1
            max_angular = 3
        elif shell_type == 'G':
            l = 2
            m = n = 1
            max_angular = 4
        else:
            sys.stderr.write('not support: {}\n'.format(shell_type))

        pwr = float(l + m + n)
        answer = math.pow(2.0, pwr)
        answer *= math.pow(pdf.Math.dbfact(2*l-1) *
                           pdf.Math.dbfact(2*m-1) *
                           pdf.Math.dbfact(2*n-1), -1.0/2.0);
        answer *= math.pow(2.0 / math.pi, 3.0 / 4.0);
        answer *= math.pow(self.exp, (pwr + 3.0/2.0) / 2.0);

        return answer
        
    # ==================================================================
    # raw data
    # ==================================================================
    def set_by_raw_data(self, odict):
        self.coef = odict.get('coef', 0.0)
        self.exp = odict.get('exp', 0.0)

    def get_raw_data(self):
        odict = {}
        odict['coef'] = self.coef
        odict['exp'] = self.exp
        return odict
    
    # ==================================================================
    # debug
    # ==================================================================
    def __str__(self):
        output = "    {0: e} {1: e}\n".format(self.exp, self.coef)
        return output

        
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
    _shell_types = ['s', 'p', 'd', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm']
    _orb_types = ['s', 'px', 'py', 'pz', 'dxy', 'dyz', 'dzx', 'dxx-yy', 'dzz',
                  'z3', 'xz2', 'yz2', '3x2y-y3', 'x3-3xy2', 'xyz', 'x2z-y2z']

    def __init__(self, *args, **kwargs):
        shell_type = 's'
        size = 0
        
        if (len(args) == 1):
            if isinstance(args[0], ContractedGTO):
                self._copy_constructer(args[0])
                return
            elif isinstance(args[0], dict):
                self.set_by_raw_data(args[0])
                return 
        
        if len(args) == 2:
            shell_type = args[0]
            size = int(args[1])
        if 'shell_type' in kwargs:
            shell_type = kwargs.get('shell_type')
        if 'size' in kwargs:
            size = kwargs.get('size')

        self.shell_type = shell_type
        list.__init__(self, [PrimitiveGTO() for x in range(size)])

        if 'pGTOs' in kwargs:
            pgtos = kwargs.get('pGTOs', [])
            self.__init__(size = len(pgtos))
            for i in range(len(pgtos)):
                self[i] = PrimitiveGTO(**(pgtos[i]))
            self.shell_type = kwargs.get('shell_type', 's')
            self.scale_factor = kwargs.get('scale_factor', 1.0)

    def _copy_constructer(self, rhs):
        self.shell_type = rhs.shell_type
        self.scale_factor = rhs.scale_factor
        list.__init__(self, [PrimitiveGTO() for x in range(len(rhs))])
        for i, pgto in enumerate(rhs):
            self[i] = PrimitiveGTO(pgto)
        
    # shell_type_id ------------------------------------------------------------
    def _get_shell_type_id(self):
        if not '_shell_type_id' in self.__dict__:
            raise
        return self._shell_type_id

    def _set_shell_type_id(self, id):
        self._shell_type_id = id
    
    shell_type_id = property(_get_shell_type_id, _set_shell_type_id)

    # shell_type ---------------------------------------------------------------
    def _get_shell_type(self):
        return self.get_shell_type(self.shell_type_id)

    def _set_shell_type(self, shell_type):
        self.shell_type_id = self.get_shell_type_id(shell_type)

    shell_type = property(_get_shell_type, _set_shell_type)

    # scale factor -------------------------------------------------------------
    def _get_scale_factor(self):
        if not '_scale_factor' in self.__dict__:
            self._scale_factor = 1.0
        return self._scale_factor

    def _set_scale_factor(self, value):
        self._scale_factor = value

    scale_factor = property(_get_scale_factor, _set_scale_factor)

    # normalize ----------------------------------------------------------------
    def normalize(self):
        shell_type = self.shell_type.upper()
        max_angular = 0
        l = m = n = 0

        if shell_type == 'S':
            l = m = n = 0
            max_angular = 0
        elif shell_type == 'P':
            l = 1
            m = n = 0
            max_angular = 1
        elif shell_type == 'D':
            l = m = 1
            n = 0
            max_angular = 2
        elif shell_type == 'F':
            l = m = n = 1
            max_angular = 3
        elif shell_type == 'G':
            l = 2
            m = n = 1
            max_angular = 4
        else:
            sys.stderr.write('not support: {}\n'.format(shell_type))

        pwr = float(l + m + n) + 3.0 / 2.0
        answer = 0.0
        for a in range(len(self)):
            coef_a = self[a].coef
            norm_a = self[a].normalize(self.shell_type)
            exp_a = self[a].exp

            for b in range(len(self)):
                coef_b = self[b].coef
                norm_b = self[b].normalize(self.shell_type)
                exp_b = self[b].exp

                trm = coef_a * coef_b
                trm *= norm_a * norm_b
                trm *= math.pow(exp_a + exp_b, -1.0 * pwr)
                answer += trm

        answer *= pdf.Math.dbfact(2*l-1) * pdf.Math.dbfact(2*m-1) * pdf.Math.dbfact(2*n-1)
        answer *= math.pow(2.0, -1.0 * float(l + m + n))
        answer *= math.pow(math.pi,  3.0 / 2.0)
        answer = math.sqrt(1.0 / answer)

        return answer
    
    # --------------------------------------------------------------------------
    @classmethod
    def get_supported_shell_types(cls):
        return cls._shell_types

    @classmethod
    def get_shell_type_id(cls, shell_type):
        """
        shell_typeに対応するidを返す

        s: 0, p: 1, d: 2
        """
        answer = None
        for i, st in enumerate(cls._shell_types):
            if shell_type == st:
                answer = i
                break
        if answer == None:
            print(shell_type)
            raise
        return answer

    @classmethod
    def get_shell_type(cls, id):
        return cls._shell_types[id]
        
    
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
    
    # ==================================================================
    # raw data
    # ==================================================================
    def set_by_raw_data(self, odict):
        list.__init__(self, [])

        self.shell_type = odict.get('shell_type')
        self.scale_factor = odict.get('scale_factor')
        pGTOs = odict.get('pGTOs', [])
        for pGTO in pGTOs:
            self.append(PrimitiveGTO(pGTO))
        
    def get_raw_data(self):
        odict = {}

        odict['shell_type'] = self.shell_type
        odict['scale_factor'] = self.scale_factor
        odict.setdefault('pGTOs', [])
        for pgto in self:
            odict['pGTOs'].append(pgto.get_raw_data())
        
        return odict
    
    # ==================================================================
    # debug
    # ==================================================================
    def __str__(self):
        output = ""
        #output += " %s %d\n".format(self.shell_type, len(self))
        output += " %d\n" % (len(self))
        for pgto in self:
            output += str(pgto)
        return output

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

        if (len(args) == 1):
            if isinstance(args[0], BasisSet):
                self._copy_constructer(args[0])
                return
            elif isinstance(args[0], dict):
                self.set_by_raw_data(args[0])
                return
        
        if (len(args) > 0):
            self.name = args[0]
            if (len(args) > 1):
                size = args[1]
        list.__init__(self, [ContractedGTO() for x in range(size)])

        if 'name' in kwargs:
            self.name = kwargs.get('name')
        if 'cGTOs' in kwargs:
            cgtos = kwargs.get('cGTOs', [])
            size = len(cgtos)
            list.__init__(self, [ContractedGTO() for x in range(size)])
            for i in range(len(cgtos)):
                self[i] = ContractedGTO(**(cgtos[i]))

    def _copy_constructer(self, rhs):
        self.name = rhs.name
        list.__init__(self, [ContractedGTO() for x in range(len(rhs))])
        for i, cgto in enumerate(rhs):
            self[i] = ContractedGTO(cgto)

    # name -------------------------------------------------------------
    def _get_name(self):
        return self._name

    def _set_name(self, name):
        self._name = str(name)

    name = property(_get_name, _set_name)

    # max_shell_type_id ------------------------------------------------
    def _get_max_shell_type_id(self):
        max_shell_type_id = 0
        for cgto in self:
            max_shell_type_id = max(max_shell_type_id, cgto.shell_type_id)
        return max_shell_type_id

    max_shell_type_id = property(_get_max_shell_type_id)

    # max_shell_type ---------------------------------------------------
    def _get_max_shell_type(self):
        return ContractedGTO.get_shell_type(self.max_shell_type_id)

    max_shell_type = property(_get_max_shell_type)
    
    # ------------------------------------------------------------------
    
    def get_number_of_AOs(self):
        answer = 0
        for cgto in self:
            st_id = cgto.shell_type_id
            answer += st_id * 2 + 1
        return answer
    
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

    # ==================================================================
    # raw data
    # ==================================================================
    def set_by_raw_data(self, odict):
        list.__init__(self, [])

        self.name = odict.get('name', '')
        cGTOs = odict.get('cGTOs', [])
        for cGTO in cGTOs:
            self.append(ContractedGTO(cGTO))
        
    def get_raw_data(self):
        odict = {}

        odict['name'] = self.name
        odict['cGTOs'] = []
        for cGTO in self:
            odict['cGTOs'].append(cGTO.get_raw_data())

        return odict
        
    # ==================================================================
    # debug
    # ==================================================================
    def get_basis2(self):
        self.sort()
        
        output = ""
        output += "%s\n" % (self.name)

        for shell_type in ContractedGTO.get_supported_shell_types():
            num_of_CGTOs = self.get_num_of_CGTOs(shell_type)
            if ((shell_type == 'spd') and (num_of_CGTOs == 0)):
                continue
            output += " %d" % (num_of_CGTOs)
        output += "\n"
        
        for cgto in self:
            output += str(cgto)
        
        return output

    def __str__(self):
        return self.get_basis2()
    
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import copy
from pdfarchive import PdfArchive
from pdfparam import PdfParam
from basisset import BasisSet, ContractedGTO, PrimitiveGTO

class OrbInfo(object):
    """
    軌道情報を管理する
    """
    def __init__(self, obj = None):
        self._orb_info = []
        self._atoms = []
        self._basissets = {}
        if isinstance(obj, PdfArchive):
            self._setup_by_db(obj)
        elif isinstance(obj, PdfParam):
            self._setup_by_param(obj)
        
    def _setup_by_db(self, db):
        assert(isinstance(db, PdfArchive))

        mol = db.get_molecule()
        # make atom and basisset list
        atom_index = 0
        for key, atom in mol.atoms():
            atom_label = self._get_atom_label(atom)
            basisset_name = db.get_basisset_name(atom_label)
            if not self._basissets.has_key(basisset_name):
                basisset = db.get_basisset(basisset_name)
                self._basissets[basisset_name] = basisset

            basisset = self._basissets[basisset_name]
            num_of_CGTOs = len(basisset)
            for CGTO_index in range(num_of_CGTOs):
                CGTO = basisset[CGTO_index]
                shell_type = CGTO.shell_type
                shell_type_id = ContractedGTO.get_shell_type_id(shell_type)
                num_of_basis_type = shell_type_id * 2 + 1
                for basis_type in range(num_of_basis_type):
                    data = {'atom_index': atom_index,
                            'basisset_name': basisset_name,
                            'CGTO_index': CGTO_index,
                            'basis_type': basis_type}
                    self._orb_info.append(data)
            self._atoms.append(copy.deepcopy(atom))
            atom_index += 1

    def _setup_by_param(self, param):
        assert(isinstance(param, PdfParam))

        mol = param.molecule
        # make atom and basisset list
        atom_index = 0
        for key, atom in mol.atoms():
            atom_label = self._get_atom_label(atom)
            basisset = param.get_basisset(atom_label)
            basisset_name = basisset.name
            if not self._basissets.has_key(basisset_name):
                self._basissets[basisset_name] = basisset

            num_of_CGTOs = len(basisset)
            for CGTO_index in range(num_of_CGTOs):
                CGTO = basisset[CGTO_index]
                shell_type = CGTO.shell_type
                shell_type_id = ContractedGTO.get_shell_type_id(shell_type)
                num_of_basis_type = shell_type_id * 2 + 1
                for basis_type in range(num_of_basis_type):
                    data = {'atom_index': atom_index,
                            'basisset_name': basisset_name,
                            'CGTO_index': CGTO_index,
                            'basis_type': basis_type}
                    self._orb_info.append(data)
            self._atoms.append(copy.deepcopy(atom))
            atom_index += 1
            
    def _get_atom_label(self, atom):
        atom_label = atom.symbol
        if len(atom.name) > 0:
            atom_label += '@' + atom.name
        return atom_label

    def get_num_of_orbitals(self):
        return len(self._orb_info)
    
    def get_atom_id(self, orb_index):
        answer = None
        if orb_index < self.get_num_of_orbitals():
            answer = self._orb_info[orb_index]['atom_index']
        return answer

    def get_atom(self, orb_index):
        atom = None
        atom_id = self.get_atom_id(orb_index)
        if atom_id != None:
            atom = self._atoms[atom_id]
        return atom
    
    def get_shell_type(self, orb_index):
        answer = None
        if orb_index < self.get_num_of_orbitals():
            basisset_name = self._orb_info[orb_index]['basisset_name']
            basisset = self._basissets[basisset_name]
            CGTO_index = self._orb_info[orb_index]['CGTO_index']
            answer = basisset[CGTO_index].shell_type
        return answer

    def get_basis_type(self, orb_index):
        answer = None
        if orb_index < self.get_num_of_orbitals():
            basisset_name = self._orb_info[orb_index]['basisset_name']
            basisset = self._basissets[basisset_name]
            CGTO_index = self._orb_info[orb_index]['CGTO_index']
            shell_type = basisset[CGTO_index].shell_type_id
            basis_type = self._orb_info[orb_index]['basis_type']
            answer = ContractedGTO.get_basis_type(shell_type, basis_type)
        return answer
        
    
if __name__ == "__main__":
    import doctest
    doctest.testmod()
        

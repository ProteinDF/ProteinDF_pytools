#!/usr/bin/env python
# -*- coding: utf-8 -*-

import types
import hashlib
import pickle
import types

import bridge
from basisset import BasisSet

class PdfParam(object):
    """
    ProteinDF parameter file (pdfparam)を扱うクラス
    """
    def __init__(self, rhs):
        self._data = {}
        if isinstance(rhs, types.DictType):
            self.__setstate__(rhs)

    def digest(self):
        md5obj = hashlib.md5()
        md5obj.update(pickle.dumps(self._data))
        return md5obj.hexdigest()

    #
    def runtypes(self):
        """
        RKS, ROKS計算なら['ALPHA']を返す。
        UKS計算なら['ALPHA', 'BETA']を返す。
        """
        runtypes = ['ALPHA']
        if self.method == 'UKS':
            runtypes = ['ALPHA', 'BETA']
        return runtypes
    
    # path ---------------------------------------------------------------------
    def _get_work_path(self):
        return self._data.get('work_path', './fl_Work')

    def _set_work_path(self, value):
        self._data['work_path'] = str(value)

    work_path = property(_get_work_path, _set_work_path)

    def get_energy_level_path(self, itr, runtype):
        """
        TODO: 'file_base_name'から取得すること
        """
        #if self._data.has_key('file_base_name'):
        #    base_name = self._data['file_base_name'].get('eigenvalues', '')
        base_name = 'eigenvalues.rks%s.vtr'
        file_name = base_name % (itr)
        return self.work_path + '/' + file_name

    def get_occ_path(self, runtype):
        """
        占有軌道情報ベクトルファイルのパスを返す
        TODO: 現在はRKS専用
        """
        file_name = ''
        if self._data.has_key('file_base_name'):
            file_name = self._data['file_base_name'].get('occupation_vtr', None)
        return self.work_path + '/' + file_name

    def get_cmat_path(self, itr, runtype):
        """
        TODO: 'file_base_name'から取得すること
        """
        #if self._data.has_key('file_base_name'):
        #    base_name = self._data['file_base_name'].get('C_matrix', '')
        base_name = 'C.rks%s.mat'
        file_name = base_name % (itr)
        return self.work_path + '/' + file_name
        
    # property -----------------------------------------------------------------
    # comment
    def _get_comment(self):
        return self._data.get('comment', '')

    def _set_comment(self, value):
        self._data['comment'] = str(value)

    comment = property(_get_comment, _set_comment)

    # method
    def _get_method(self):
        return self._data.get('method', None)

    def _set_method(self, value):
        value = value.upper()
        assert(value == 'RKS' or value == 'UKS' or value == 'ROKS')
        self._data['method'] = value

    method = property(_get_method, _set_method)

    # guess
    def _get_guess(self):
        return self._data.get('guess', None)

    def _set_guess(self, value):
        value = value.upper()
        self._data['guess'] = value

    guess = property(_get_guess, _set_guess)
    
    # xc_functional
    def _get_xc_functional(self):
        return self._data.get('xc_functional', None)

    def _set_xc_functional(self, value):
        value = value.upper()
        self._data['xc_functional'] = value

    xc_functional = property(_get_xc_functional, _set_xc_functional)

    # num_of_atoms
    def _get_num_of_atoms(self):
        return self._data.get('num_of_atoms', None)

    def _set_num_of_atoms(self, value):
        self._data['num_of_atoms'] = int(value)

    num_of_atoms = property(_get_num_of_atoms, _set_num_of_atoms)
    
    # num_of_AOs
    def _get_num_of_AOs(self):
        return self._data.get('num_of_AOs', None)

    def _set_num_of_AOs(self, value):
        self._data['num_of_AOs'] = int(value)

    num_of_AOs = property(_get_num_of_AOs, _set_num_of_AOs)

    # num_of_MOs
    def _get_num_of_MOs(self):
        return self._data.get('num_of_MOs', None)

    def _set_num_of_MOs(self, value):
        if value == None or len(value) == 0:
            value = 0
        self._data['num_of_MOs'] = int(value)

    num_of_MOs = property(_get_num_of_MOs, _set_num_of_MOs)

    # iterations
    def _get_iterations(self):
        return self._data.get('iterations', None)

    def _set_iterations(self, value):
        if value == None or len(value) == 0:
            value = 0
        self._data['iterations'] = int(value)

    iterations = property(_get_iterations, _set_iterations)
        
    # max_iterations
    def _get_max_iterations(self):
        return self._data.get('max_iterations', None)

    def _set_max_iterations(self, value):
        self._data['max_iterations'] = int(value)

    max_iterations = property(_get_max_iterations, _set_max_iterations)

    # SCF converged
    @property
    def scf_converged(self):
        return self._data.get('scf_converged', False)
    
    # molecule
    def _get_molecule(self):
        return self._data.get('molecule', bridge.AtomGroup())

    def _set_moleucule(self, value):
        assert(isinstance(value, bridge.AtomGroup))
        self._data['molecule'] = value

    molecule = property(_get_molecule, _set_moleucule)

    # basisset
    def set_basisset(self, atom_label, basisset):
        assert(isinstance(atom_label, str))
        assert(isinstance(basisset, BasisSet))
        if not self._data.has_key('basisset'):
            self._data.setdefault('basisset', {})
        self._data['basisset'][atom_label] = basisset

    def get_basisset_atomlabels(self):
        """
        原子(ラベル)名のリストを返す
        """
        self._data.setdefault('basisset', {})
        for i in self._data['basisset'].keys():
            yield i

    def get_basisset(self, atom_label):
        """
        原子(ラベル)名のBasisSetがあれば、そのBasisSetオブジェクトを返す
        """
        assert(isinstance(atom_label, types.StringType))
        answer = BasisSet()
        if self._data.has_key('basisset'):
            answer = self._data['basisset'].get(atom_label, BasisSet())
        return answer
        
    # TEs
    def _get_TEs(self):
        return self._data.get('TEs', None)

    def _set_TEs(self, value):
        self._data['TEs'] = value

    TEs = property(_get_TEs, _set_TEs)
        
    # --------------------------------------------------------------------------
    def __setstate__(self, rhs):
        assert(isinstance(rhs, types.DictType))
        rhs = self._alias_conversion(rhs)
    
        self.method = rhs.get('method', None)
        self.guess = rhs.get('guess', None)
        self.xc_functional = rhs.get('xc_functional', None)
        self.convergence_target = rhs.get('convergence_target', None)
        self.num_of_atoms = rhs.get('num_of_atoms', None)
        self.num_of_AOs = rhs.get('num_of_AOs', None)
        self.num_of_MOs = rhs.get('num_of_MOs', None)

        self.max_iteration = rhs.get('max_iteration', None)
        self.iterations = rhs.get('iterations', None)

        # basis set
        self._basissets = {}
        for atom_label, basisset_state in rhs['basis_sets'].items():
            basisset = BasisSet()
            basisset.__setstate__(basisset_state)
            self.set_basisset(atom_label, basisset)
            
        # coordinates
        index = 0
        molecule = bridge.AtomGroup()
        for atomgroup_label, atomgroup_data in rhs['coordinates'].items():
            for atom_data in atomgroup_data:
                symbol = atom_data['symbol']
                xyz = bridge.Position(atom_data['xyz'])
                z = atom_data['charge']
                label = atom_data['label']
                atom = bridge.Atom(symbol = symbol,
                                   position = xyz,
                                   charge = z,
                                   name = label)
                molecule.set_atom(index, atom)
                index += 1
        self.molecule = molecule

        # total energy
        self.TEs = rhs.get('TEs', {})

        # control
        self._data['file_base_name'] = {}
        control = rhs.get('control', None)
        if isinstance(control, types.DictType):
            file_base_name = control.get('file_base_name', None)
            if file_base_name is not None:
                for key, value in file_base_name.items():
                    self._data['file_base_name'][key] = value
            self._data['scf_converged'] = control.get('scf_converged', False)
        
    def __getstate__(self):
        return self._data

    def _alias_conversion(self, rhs):
        """
        キーおよび値に別名(alias)が使用されていた場合、処理すべきキーに変換する
        """
        for k, v in rhs.items():
            if k == 'method':
                v = v.upper()
                if v == 'NSP':
                    v = 'RKS'
                elif v == 'SP':
                    v = 'UKS'
                rhs[k] = v

            if k == 'scf-start-guess':
                rhs['guess'] = v
                rhs.pop(k)

            if k == 'xc-potential':
                rhs['xc_functional'] = v
                rhs.pop(k)
                
            if k == 'max-iteration':
                rhs['max_iteration'] = v
                rhs.pop(k)

            if k == 'num_of_iterations':
                rhs['iterations'] = v
                rhs.pop(k)

            if k == 'TE':
                rhs['TEs'] = v
                rhs.pop(k)

        return rhs

    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
    

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import types
import hashlib
import pickle
import types
import math

import bridge
from basisset import BasisSet

class PdfParam(object):
    """
    ProteinDF parameter file (pdfparam)を扱うクラス
    """
    def __init__(self, rhs = None):
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
        """
        runtypes = []
        if self.method == 'rks':
            runtypes = ['rks']
        elif self.method == 'uks':
            runtypes = ['uks_alpha', 'uks_beta']
        elif self.method == 'roks':
            runtypes = ['roks_close', 'roks_open']
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
        TODO: 現在はrks専用
        """
        file_name = ''
        if self._data.has_key('file_base_name'):
            file_name = self._data['file_base_name'].get('occupation_vtr', None)
        file_name = file_name % (runtype)
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

    # step_control
    def _get_step_control(self):
        return self._data.get('step_control', '')

    def _set_step_control(self, values):
        self._data['step_control'] = str(values)

    step_control = property(_get_step_control, _set_step_control)

    # comment
    def _get_comment(self):
        return self._data.get('comment', '')

    def _set_comment(self, value):
        self._data['comment'] = str(value)

    comment = property(_get_comment, _set_comment)

    # method
    def _get_method(self):
        return self._data.get('method', 'rks')

    def _set_method(self, value):
        value = value.lower()
        assert(value == 'rks' or value == 'uks' or value == 'roks')
        self._data['method'] = value

    method = property(_get_method, _set_method)

    # guess
    def _get_guess(self):
        return self._data.get('guess', None)

    def _set_guess(self, value):
        value = str(value)
        self._data['guess'] = value

    guess = property(_get_guess, _set_guess)

    # orbital_independence_threshold
    def _get_orbital_independence_threshold(self):
        return self._data.get('orbital_independence_threshold', 0.007)

    def _set_orbital_independence_threshold(self, value):
        value = float(value)
        self._data['orbital_independence_threshold'] = value

    orbital_independence_threshold = property(_get_orbital_independence_threshold,
                                              _set_orbital_independence_threshold)

    # scf_acceleration
    def _get_scf_acceleration(self):
        return self._data.get('scf_acceleration', 'damping')

    def _set_scf_acceleration(self, value):
        value = str(value)
        self._data['scf_acceleration'] = value

    scf_acceleration = property(_get_scf_acceleration,
                                _set_scf_acceleration)

    # scf_acceleration/damping_factor
    def _get_scf_acceleration_damping_factor(self):
        return self._data.get('scf_acceleration_damping_factor', 0.85)

    def _set_scf_acceleration_damping_factor(self, value):
        value = float(value)
        self._data['scf_acceleration_damping_factor'] = value

    scf_acceleration_damping_factor = property(_get_scf_acceleration_damping_factor,
                                               _set_scf_acceleration_damping_factor)
        
    # xc_functional
    def _get_xc_functional(self):
        return self._data.get('xc_functional', 'SVWN')

    def _set_xc_functional(self, value):
        value = str(value)
        self._data['xc_functional'] = value

    xc_functional = property(_get_xc_functional, _set_xc_functional)

    def _is_fitting_xc(self):
        """
        使用している交換相関汎関数がRI法を利用するならばTrueを返す
        """
        last_char = self.xc_functional[-1]
        return (last_char == '~')
    
    # num_of_atoms
    def _get_num_of_atoms(self):
        return self._data.get('num_of_atoms', None)

    def _set_num_of_atoms(self, value):
        self._data['num_of_atoms'] = int(value)

    num_of_atoms = property(_get_num_of_atoms, _set_num_of_atoms)

    # num_of_electrons
    def _get_num_of_electrons(self):
        value = self._data.get('num_of_electrons', None)
        if (value == None):
            value = self.molecule.sum_of_atomic_number()
            self.num_of_electrons = value
        return value

    def _set_num_of_electrons(self, value):
        self._data['num_of_electrons'] = value

    num_of_electrons = property(_get_num_of_electrons, _set_num_of_electrons)
    
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
        if value == None:
            value = 0
        self._data['num_of_MOs'] = int(value)

    num_of_MOs = property(_get_num_of_MOs, _set_num_of_MOs)

    # occupation_level
    def _get_occupation_level(self):
        value = self._data.get('occupation_level', None)
        if value == None:
            value = '[1 - {0}]'.format(int(self.num_of_electrons / 2))
            self.occupation_level = value 
        return value

    def _set_occupation_level(self, level):
        self._data['occupation_level'] = level

    occupation_level = property(_get_occupation_level, _set_occupation_level)

    # iterations
    def _get_iterations(self):
        return self._data.get('iterations', None)

    def _set_iterations(self, value):
        if value == None:
            value = 0
        self._data['iterations'] = int(value)

    iterations = property(_get_iterations, _set_iterations)
        
    # max_iterations
    def _get_max_iterations(self):
        return self._data.get('max_iterations', 100)

    def _set_max_iterations(self, value):
        self._data['max_iterations'] = int(value)

    max_iterations = property(_get_max_iterations, _set_max_iterations)

    # convergence_target
    def _get_convergence_target(self):
        value = self._data.get('convergence_target', None)
        if value == None:
            if (self._is_fitting_xc() == True):
                value = 'density'
            else:
                value = 'density_matrix'
            self.convergence_target = value
        return value

    def _set_convergence_target(self, target):
        self._data['convergence_target'] = target

    convergence_target = property(_get_convergence_target, _set_convergence_target)

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
    def get_basisset_atomlabels(self):
        """
        原子(ラベル)名のリストを返す
        """
        self._data.setdefault('basisset', {})
        for i in self._data['basisset'].keys():
            yield i

    def get_basisset(self, atom_label):
        return self._get_basisset_common(atom_label, 'basisset')

    def set_basisset(self, atom_label, basisset):
        self._set_basisset_common(atom_label, basisset, 'basisset')

    def get_basisset_j(self, atom_label):
        return self._get_basisset_common(atom_label, 'basisset_j')

    def set_basisset_j(self, atom_label, basisset):
        self._set_basisset_common(atom_label, basisset, 'basisset_j')

    def get_basisset_xc(self, atom_label):
        return self._get_basisset_common(atom_label, 'basisset_xc')

    def set_basisset_xc(self, atom_label, basisset):
        self._set_basisset_common(atom_label, basisset, 'basisset_xc')

    def _get_basisset_common(self, atom_label, key):
        """
        原子(ラベル)名のBasisSetがあれば、そのBasisSetオブジェクトを返す
        """
        assert(isinstance(atom_label, types.StringType))
        answer = BasisSet()
        if self._data.has_key(key):
            answer = self._data[key].get(atom_label, BasisSet())
        return answer

    def _set_basisset_common(self, atom_label, basisset, key):
        """
        原子(ラベル)名にBasisSetオブジェクトを設定する
        """
        assert(isinstance(atom_label, str))
        assert(isinstance(basisset, BasisSet))
        self._data.setdefault(key, {})
        self._data[key][atom_label] = basisset
        
    # TEs
    def _get_TEs(self):
        return self._data.get('TEs', None)

    def _set_TEs(self, value):
        self._data['TEs'] = value

    TEs = property(_get_TEs, _set_TEs)
        
    # --------------------------------------------------------------------------
    def get_inputfile_contents(self):
        """
        ProteinDF入力ファイル(fl_Userinput)の内容を作成する
        """
        output = ""
        output += ">>>>MAIN\n"
        output += "    step-control = [%s]\n" % (self.step_control)
        output += "\n"
        output += ">>>>MODEL\n"
        output += "    scf-start-guess = %s\n" % (self.guess)
        output += "    method = %s\n" % (self.method)
        if self.method == 'rks':
            output += "    method/nsp/electron-number = %d\n" % (
                self.num_of_electrons)
            output += "    method/nsp/occlevel = %s\n" % (
                self.occupation_level)
        else:
            raise
        output += "    max_iteration = {0}\n".format(self.max_iterations)
        output += "    orbital_independence_threshold  = {0}\n".format(self.orbital_independence_threshold)
        output += "    scf-acceleration/damping/damping-type = {0}\n".format(self.convergence_target)
        output += "    scf_acceleration = {0}\n".format(self.scf_acceleration)
        output += "    scf_acceleration/damping/damping_factor = {0}\n".format(self.scf_acceleration_damping_factor)
        output += "    xc_potential = {0}\n".format(self.xc_functional)
        output += "\n"
        output += ">>>>MOLECULE\n"
        output += "    geometry/cartesian/unit = angstrom\n"
        output += "    geometry/cartesian/input = {\n"
        output += "        " + self._get_input_geometry(self.molecule).replace('\n', '\n        ')
        output += "}\n"
        output += "\n"

        output += "    basis-set/orbital = {\n"
        basisset_str = self._get_input_basisset()
        basisset_str = basisset_str.replace('\n', '\n        ')
        output += "        " + basisset_str
        output += "}\n\n"
        
        output += "    basis-set/density-auxiliary = {\n"
        basisset_str = self._get_input_basisset_j()
        basisset_str = basisset_str.replace('\n', '\n        ')
        output += "        " + basisset_str
        output += "}\n\n"
        
        output += "    basis-set/exchange-auxiliary = {\n"
        basisset_str = self._get_input_basisset_xc()
        basisset_str = basisset_str.replace('\n', '\n        ')
        output += "        " + basisset_str
        output += "}\n\n"
        
        return output

    def _get_input_geometry(self, atom_group = None):
        output = ""
        if (atom_group == None):
            output += self._get_geometry_input(self.molecule)
        else:
            for key, group in atom_group.groups():
                output += self._get_input_geometry(group)
                output += "\n"
            for key, atom in atom_group.atoms():
                if (math.fabs(atom.charge) < 1.0E-5):
                    output += "%2s % 10.6f % 10.6f % 10.6f // %s:%s\n" %(
                        atom.symbol,
                        atom.xyz.x,
                        atom.xyz.y,
                        atom.xyz.z,
                        atom.path,
                        atom.name)
                else:
                    output += "%2s % 10.6f % 10.6f % 10.6f % 10.6f\n" %(
                        atom.symbol,
                        atom.xyz.x,
                        atom.xyz.y,
                        atom.xyz.z,
                        atom.charge)
        return output

    def _get_input_basisset(self):
        output = ""
        for atom_kind in self._pickup_atom_kinds():
            basisset = self.get_basisset(atom_kind)
            output += '%s = "%s"\n' % (atom_kind, basisset.name)
        return output

    def _get_input_basisset_j(self):
        output = ""
        for atom_kind in self._pickup_atom_kinds():
            basisset = self.get_basisset_j(atom_kind)
            output += '%s = "%s"\n' % (atom_kind, basisset.name)
        return output

    def _get_input_basisset_xc(self):
        output = ""
        for atom_kind in self._pickup_atom_kinds():
            basisset = self.get_basisset_xc(atom_kind)
            output += '%s = "%s"\n' % (atom_kind, basisset.name)
        return output
        
    def _pickup_atom_kinds(self, atom_group = None):
        answer = set()
        if (atom_group == None):
            answer = self._pickup_atom_kinds(self.molecule)
        else:
            for key, group in atom_group.groups():
                answer.update(self._pickup_atom_kinds(group))
            for key, atom in atom_group.atoms():
                atom_symbol = atom.symbol
                label = ""

                atom_kind = atom_symbol
                if (len(label) != 0):
                    atom_kind = "%s@%s" % (atom_symbol, label)
                answer.add(atom_kind)
        return answer
        
        
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
            basisset = BasisSet(**basisset_state)
            self.set_basisset(atom_label,
                              basisset)
            
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
                v = v.lower()
                if v == 'nsp':
                    v = 'rks'
                elif v == 'sp':
                    v = 'uks'
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
    

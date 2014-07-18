#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import hashlib
import pickle
import types
import math
import copy
import pprint

import bridge
import pdf

class PdfParam(object):
    """
    ProteinDF parameter file (pdfparam)を扱うクラス
    """
    def __init__(self, rhs = None):
        self._data = {}
        if isinstance(rhs, PdfParam):
            self._data = copy.deepcopy(rhs._data)
        elif isinstance(rhs, dict):
            self._set_by_dict(rhs)

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
            #runtypes = ['uks_alpha', 'uks_beta']
            runtypes = ['uks-alpha', 'uks-beta']
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
        if 'file_base_name' in self._data:
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

    def get_ovpmat_path(self):
        filename = self._get_file_base_name('Spq_matrix')
        return os.path.join(self.work_path, filename)

    def get_clomat_path(self, runtype='rks', itr=-1):
        '''
        return localized matrix (C_lo) path
        '''
        if itr == -1:
            itr = self.lo_num_of_iterations
        
        filename = self._get_file_base_name('Clo_matrix')
        filename = filename.replace('%s', '{}{}'.format(runtype, itr))
        return os.path.join(self.work_path, filename)
        
    def _get_file_base_name(self, key):
        return self._data['file_base_name'].get(key, '')
        
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
        value = str(value).lower()
        if not (value == 'rks' or value == 'uks' or value == 'roks'):
            print(value)
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

    # orbital_independence_threshold_canonical
    def _get_orbital_independence_threshold_canonical(self):
        return self._data.get('orbital_independence_threshold_canonical', 0.007)

    def _set_orbital_independence_threshold_canonical(self, value):
        value = float(value)
        self._data['orbital_independence_threshold_canonical'] = value

    orbital_independence_threshold_canonical = property(_get_orbital_independence_threshold_canonical,
                                                        _set_orbital_independence_threshold_canonical)
        
    # orbital_independence_threshold_lowdin
    def _get_orbital_independence_threshold_lowdin(self):
        return self._data.get('orbital_independence_threshold_lowdin', 0.007)

    def _set_orbital_independence_threshold_lowdin(self, value):
        value = float(value)
        self._data['orbital_independence_threshold_lowdin'] = value

    orbital_independence_threshold_lowdin = property(_get_orbital_independence_threshold_lowdin,
                                                     _set_orbital_independence_threshold_lowdin)
        
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

    # J_engine
    def _get_j_engine(self):
        return self._data.get('j_engine', 'RI_J')

    def _set_j_engine(self, value):
        value = str(value)
        self._data['j_engine'] = value

    j_engine = property(_get_j_engine, _set_j_engine)
    
    # K_engine
    def _get_k_engine(self):
        return self._data.get('k_engine', 'conventional')

    def _set_k_engine(self, value):
        value = str(value)
        self._data['k_engine'] = value

    k_engine = property(_get_k_engine, _set_k_engine)
    
    # XC_engine
    def _get_xc_engine(self):
        return self._data.get('xc_engine', 'grid')

    def _set_xc_engine(self, value):
        value = str(value)
        self._data['xc_engine'] = value

    xc_engine = property(_get_xc_engine, _set_xc_engine)
    
    # gridfree/dedicated_basis
    def _get_gridfree_dedicated_basis(self):
        return self._data.get('gridfree_dedicated_basis', 'no')

    def _set_gridfree_dedicated_basis(self, value):
        value = bool(value)
        if value == True:
            value = 'yes'
        else:
            value = 'no'
        self._data['gridfree_dedicated_basis'] = value

    gridfree_dedicated_basis = property(_get_gridfree_dedicated_basis,
                                        _set_gridfree_dedicated_basis)
        
    # gridfree/orthogonalize_method
    def _get_gridfree_orthogonalize_method(self):
        return self._data.get('gridfree_orthogonalize_method', 'lowdin')

    def _set_gridfree_orthogonalize_method(self, value):
        value = str(value)
        self._data['gridfree_orthogonalize_method'] = value

    gridfree_orthogonalize_method = property(_get_gridfree_orthogonalize_method,
                                             _set_gridfree_orthogonalize_method)
    
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
            value += self.molecule.charge
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
        return self._data.get('num_of_iterations', None)

    def _set_iterations(self, value):
        if value == None:
            value = 0
        self._data['num_of_iterations'] = int(value)

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

    # convergence_threshold_energy
    def _get_convergence_threshold_energy(self):
        return self._data.get('convergence_threshold_energy', 1.0E-4)

    def _set_convergence_threshold_energy(self, value):
        self._data['convergence_threshold_energy'] = float(value)

    convergence_threshold_energy = property(_get_convergence_threshold_energy,
                                            _set_convergence_threshold_energy)

    # convergence_threshold
    def _get_convergence_threshold(self):
        return self._data.get('convergence_threshold', 1.0E-3)

    def _set_convergence_threshold(self, value):
        self._data['convergence_threshold'] = float(value)

    convergence_threshold = property(_get_convergence_threshold,
                                     _set_convergence_threshold)
   
    # convergence_type
    def _get_convergence_type(self):
        return self._data.get('convergence_type', 'density')

    def _set_convergence_type(self, value):
        self._data['convergence_type'] = str(value)

    convergence_type = property(_get_convergence_type,
                                _set_convergence_type)

    # level_shift
    def _get_level_shift(self):
        return self._data.get('level_shift', False)

    def _set_level_shift(self, value):
        self._data['level_shift'] = bool(value)

    level_shift = property(_get_level_shift, _set_level_shift)

    def _get_level_shift_start_iteration(self):
        return self._data.get('level_shift_start_iterartion', 1)

    def _set_level_shift_start_iteration(self, value):
        self._data['level_shift_start_iteration'] = int(value)

    level_shift_start_iteration = property(_get_level_shift_start_iteration,
                                           _set_level_shift_start_iteration)

    def _get_level_shift_virtual_mo(self):
        return self._data.get('level_shift_virtual_mo', 0.0)

    def _set_level_shift_virtual_mo(self, value):
        self._data['level_shift_virtual_mo'] = float(value)

    level_shift_virtual_mo = property(_get_level_shift_virtual_mo,
                                      _set_level_shift_virtual_mo)
    
    
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

    # basisset ---------------------------------------------------------
    def get_basisset_atomlabels(self):
        """
        原子(ラベル)名のリストを返す
        """
        self._data.setdefault('basisset', {})
        for i in self._data['basisset'].keys():
            yield i

    def get_basisset_name(self, atomlabel):
        atomlabel = bridge.Utils.byte2str(atomlabel)
        self._data.setdefault('basisset_name', {})
        return self._data['basisset_name'].get(atomlabel, '')

    def set_basisset_name(self, atomlabel, value):
        atomlabel = bridge.Utils.byte2str(atomlabel)
        self._data.setdefault('basisset_name', {})
        self._data['basisset_name'][atomlabel] = bridge.Utils.byte2str(value)


    def get_basisset_j_name(self, atomlabel):
        atomlabel = bridge.Utils.byte2str(atomlabel)
        self._data.setdefault('basisset_j_name', {})
        return self._data['basisset_j_name'].get(atomlabel, '')

    def set_basisset_j_name(self, atomlabel, value):
        atomlabel = bridge.Utils.byte2str(atomlabel)
        self._data.setdefault('basisset_j_name', {})
        self._data['basisset_j_name'][atomlabel] = bridge.Utils.byte2str(value)


    def get_basisset_xc_name(self, atomlabel):
        atomlabel = bridge.Utils.byte2str(atomlabel)
        self._data.setdefault('basisset_xc_name', {})
        return self._data['basisset_xc_name'].get(atomlabel, '')

    def set_basisset_xc_name(self, atomlabel, value):
        atomlabel = bridge.Utils.byte2str(atomlabel)
        self._data.setdefault('basisset_xc_name', {})
        self._data['basisset_xc_name'][atomlabel] = bridge.Utils.byte2str(value)


    def get_basisset_gridfree_name(self, atomlabel):
        atomlabel = bridge.Utils.byte2str(atomlabel)
        self._data.setdefault('basisset_gridfree_name', {})
        return self._data['basisset_gridfree_name'].get(atomlabel, '')

    def set_basisset_gridfree_name(self, atomlabel, value):
        atomlabel = bridge.Utils.byte2str(atomlabel)
        self._data.setdefault('basisset_gridfree_name', {})
        self._data['basisset_gridfree_name'][atomlabel] = bridge.Utils.byte2str(value)

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

    def get_basisset_gridfree(self, atom_label):
        return self._get_basisset_common(atom_label, 'basisset_gridfree')
    
    def set_basisset_gridfree(self, atom_label, basisset):
        self._set_basisset_common(atom_label, basisset, 'basisset_gridfree')

    def _get_basisset_common(self, atom_label, key):
        """
        原子(ラベル)名のBasisSetがあれば、そのBasisSetオブジェクトを返す
        """
        atom_label = bridge.Utils.byte2str(atom_label)
        answer = pdf.BasisSet()
        if key in self._data:
            answer = self._data[key].get(atom_label, pdf.BasisSet())
        return answer

    def _set_basisset_common(self, atom_label, basisset, key):
        """
        原子(ラベル)名にBasisSetオブジェクトを設定する
        """
        atom_label = bridge.Utils.byte2str(atom_label)
        self._data.setdefault(key, {})
        if isinstance(basisset, str) == True:
            basis2 = pdf.Basis2()
            if (key == 'basis_set') or (key == 'basis_set_gridfree'):
                self._data[key][atom_label] = basis2.get_basisset(basisset)
            elif key == 'basis_set_j':
                self._data[key][atom_label] = basis2.get_basisset_j(basisset)
            elif key == 'basis_set_xc':
                self._data[key][atom_label] = basis2.get_basisset_xc(basisset)
            else:
                raise                              
                
        elif isinstance(basisset, pdf.BasisSet) == True:
            self._data[key][atom_label] = basisset
        else:
            raise "type mispatch"
        
    # TEs
    def _get_TEs(self):
        return self._data.get('TEs', None)

    def _set_TEs(self, value):
        self._data['TEs'] = value

    TEs = property(_get_TEs, _set_TEs)

    # counterpoise
    def _get_counterpoise(self):
        return self._data.get('counterpoise', False)

    counterpoise = property(_get_counterpoise)

    # LO ===============================================================
    def _get_lo_satisfied(self):
        yn = False
        if self._data.get('lo/satisfied', 'no').upper() == 'YES':
            yn = True
        return yn

    lo_satisfied = property(_get_lo_satisfied)

    def _get_lo_num_of_iterations(self):
        return self._data.get('lo/num_of_iterations', None)
    
    lo_num_of_iterations = property(_get_lo_num_of_iterations)
    
    # force ============================================================
    def get_force(self, atom_index):
        atom_index = int(atom_index)
        assert(0 <= atom_index < self.num_of_atoms)

        force_mat = self._data.get('force', None)
        force = None
        if force_mat != None:
            if atom_index < len(force_mat):
                force = force_mat[atom_index]
                assert(len(force) == 3)

        return force
            

    def set_force(self, atom_index, fx, fy, fz):
        atom_index = int(atom_index)
        fx = float(fx)
        fy = float(fy)
        fz = float(fz)
        assert(0 <= atom_index < self.num_of_atoms)
        self._data.setdefault('force', [[] for x in range(self.num_of_atoms)])
        self._data['force'][atom_index] = [fx, fy, fz]

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
        output += "    orbital_independence_threshold/canonical = {0}\n".format(self.orbital_independence_threshold_canonical)
        output += "    orbital_independence_threshold/lowdin = {0}\n".format(self.orbital_independence_threshold_lowdin)
        output += "    scf_acceleration/damping/damping_type = {0}\n".format(self.convergence_target)
        output += "    scf_acceleration = {0}\n".format(self.scf_acceleration)
        output += "    scf_acceleration/damping/damping_factor = {0}\n".format(self.scf_acceleration_damping_factor)
        output += "    convergence/threshold_energy = {0}\n".format(self.convergence_threshold_energy)
        output += "    convergence/threshold = {0}\n".format(self.convergence_threshold)
        output += "    convergence/type = {0}\n".format(self.convergence_type)
        output += "    xc_functional = {0}\n".format(self.xc_functional)
        output += "    J_engine = {0}\n".format(self.j_engine)
        output += "    K_engine = {0}\n".format(self.k_engine)
        output += "    XC_engine = {0}\n".format(self.xc_engine)
        output += "    gridfree/dedicated_basis = {0}\n".format(self.gridfree_dedicated_basis)
        output += "    gridfree/orthogonalize_method = {0}\n".format(self.gridfree_orthogonalize_method)
        output += "    level_shift = {0}\n".format(self.level_shift)
        output += "    level_shift/start_iteration = {0}\n".format(self.level_shift_start_iteration)
        output += "    level_shift/virtual_mo = {0}\n".format(self.level_shift_virtual_mo)
        output += "\n"
        output += ">>>>MOLECULE\n"
        output += "    geometry/cartesian/unit = angstrom\n"
        output += "    geometry/cartesian/input = {\n"
        output += "        " + self._get_input_geometry(self.molecule).replace('\n', '\n        ').rstrip()
        output += "\n"
        output += "        }\n"
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

        output += "    basis-set/gridfree = {\n"
        basisset_str = self._get_input_basisset_gridfree()
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
                output_sub = self._get_input_geometry(group)
                output += output_sub
                if len(output_sub) > 0:
                    output += "\n"
            for key, atom in atom_group.atoms():
                label = ""
                if len(atom.label) > 0:
                    label = "@{}".format(label)
                    label = re.sub("\s+", "_", label)

                charge = ""
                if math.fabs(atom.charge) > 1.0E-5:
                    charge = str(atom.charge)

                atom_symbol = atom.symbol
                if atom_symbol != 'X':
                    output += "{atom_symbol}{label} {x: f} {y: f} {z: 6f} // {path}:{name}\n".format(
                        atom_symbol = atom_symbol,
                        label = label,
                        x = atom.xyz.x,
                        y = atom.xyz.y,
                        z = atom.xyz.z,
                        path = atom.path,
                        name = atom.name)
                else:
                    output += "{atom_symbol}{label} {x: f} {y: f} {z: 6f} {charge} // {path}:{name}\n".format(
                        atom_symbol = atom_symbol,
                        label = label,
                        x = atom.xyz.x,
                        y = atom.xyz.y,
                        z = atom.xyz.z,
                        charge = charge,
                        path = atom.path,
                        name = atom.name)
                    
        return output

    def _get_input_basisset(self):
        output = ""
        for atom_kind in self._pickup_atom_kinds():
            output += '%s = "%s"\n' % (atom_kind, self.get_basisset_name(atom_kind))
        return output

    def _get_input_basisset_j(self):
        output = ""
        for atom_kind in self._pickup_atom_kinds():
            output += '%s = "%s"\n' % (atom_kind, self.get_basisset_j_name(atom_kind))
        return output

    def _get_input_basisset_xc(self):
        output = ""
        for atom_kind in self._pickup_atom_kinds():
            output += '%s = "%s"\n' % (atom_kind, self.get_basisset_xc_name(atom_kind))
        return output
        
    def _get_input_basisset_gridfree(self):
        output = ""
        for atom_kind in self._pickup_atom_kinds():
            output += '%s = "%s"\n' % (atom_kind, self.get_basisset_gridfree_name(atom_kind))
        return output

    def _pickup_atom_kinds(self, atom_group = None):
        '''
        counterpoise =
        '''
        answer = set()
        if (atom_group == None):
            answer = self._pickup_atom_kinds(self.molecule)
        else:
            for key, group in atom_group.groups():
                answer.update(self._pickup_atom_kinds(group))
            for key, atom in atom_group.atoms():
                atom_symbol = atom.symbol
                if (atom_symbol != 'X') or (self.counterpoise == True):
                    atom_label = atom.label
                    
                    atom_kind = atom_symbol
                    if (len(atom_label) != 0):
                        atom_kind = "%s@%s" % (atom_symbol, label)
                answer.add(atom_kind)
        return answer
        
        
    # --------------------------------------------------------------------------
    def _set_by_dict(self, rhs):
        assert(isinstance(rhs, dict))
        rhs = self._alias_conversion(rhs)
    
        self.method = rhs.get('method', self.method)
        self.guess = rhs.get('guess', self.guess)
        self.xc_functional = rhs.get('xc_functional', self.xc_functional)
        self.convergence_target = rhs.get('convergence_target', self.convergence_target)
        self.num_of_atoms = rhs.get('num_of_atoms', self.num_of_atoms)
        self.num_of_AOs = rhs.get('num_of_AOs', self.num_of_AOs)
        self.num_of_MOs = rhs.get('num_of_MOs', self.num_of_MOs)

        self.max_iterations = rhs.get('max_iterations', self.max_iterations)
        self.iterations = rhs.get('num_of_iterations', self.iterations)

        # basis set
        self._basissets = {}
        for atom_label, basisset_state in rhs['basis_set'].items():
            basisset = pdf.BasisSet(**basisset_state)
            self.set_basisset(atom_label,
                              basisset)
            
        # coordinates
        index = 0
        molecule = bridge.AtomGroup()
        for atomgroup_label, atomgroup_data in rhs['coordinates'].items():
            subgroup = bridge.AtomGroup()
            for atom_data in atomgroup_data:
                symbol = atom_data['symbol']
                xyz = bridge.Position(atom_data['xyz'])
                z = atom_data['charge']
                label = atom_data['label']
                atom = bridge.Atom(symbol = symbol,
                                   position = xyz,
                                   charge = z,
                                   name = label)
                subgroup.set_atom(index, atom)
                index += 1
            molecule.set_group(atomgroup_label, subgroup)
        self.molecule = molecule

        # total energy
        self.TEs = rhs.get('TEs', {})

        # force
        force_dat = rhs.get('force', None)
        if force_dat != None:
            for atom_index in range(len(force_dat)):
                force = force_dat[atom_index]
                self.set_force(atom_index, force[0], force[1], force[2])

        # lo
        self._data['lo/satisfied'] = rhs.get('lo/satisfied', 'no')
        self._data['lo/num_of_iterations'] = rhs.get('lo/num_of_iterations', None)
                
        # control
        self._data['file_base_name'] = {}
        control = rhs.get('control', None)
        if isinstance(control, dict):
            file_base_name = control.get('file_base_name', None)
            if file_base_name is not None:
                for key, value in file_base_name.items():
                    self._data['file_base_name'][key] = value
            self._data['scf_converged'] = control.get('scf_converged', False)
        
    # def __getstate__(self):
    #     return self._data

    def _alias_conversion(self, rhs):
        """
        キーおよび値に別名(alias)が使用されていた場合、処理すべきキーに変換する
        """
        answer = {}
        for k, v in rhs.items():
            if k == 'method':
                v = v.lower()
                if v == 'nsp':
                    v = 'rks'
                elif v == 'sp':
                    v = 'uks'
                answer[k] = v
            elif k == 'scf-start-guess':
                answer['guess'] = v
            elif k == 'xc-potential':
                answer['xc_functional'] = v
            elif k == 'max-iteration':
                answer['max_iteration'] = v
            elif k == 'TE':
                answer['TEs'] = v
            else:
                answer[k] = v

        return answer

    
if __name__ == '__main__':
    import doctest
    doctest.testmod()
    

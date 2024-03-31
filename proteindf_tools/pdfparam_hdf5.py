#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2019 The ProteinDF development team.
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

from .vector import Vector
from .matrix import Matrix, SymmetricMatrix
from .pdfparam_object import PdfParamObject
from .pdfcommon import run_pdf
from .functions import deprecated

import h5py
import os
import numpy
import traceback

import logging
logger = logging.getLogger(__name__)


class PdfParam_H5(PdfParamObject):
    """ ProteinDF parameter object
    
    >>> # usage
    >>> pdfparam = PdfParam_H5()
    >>> pdfparam.filepath = "test.h5"
    
    >>> pdfparam = PdfParam_H5("test.h5")
    """

    def __init__(self, rhs=None):
        super().__init__(rhs)
        
        self._state = {}

        if (isinstance(rhs, str)):
            self.open(rhs)
            self.load_model()
            self.get_TEs()


    def __del__(self):
        self._close()

    # --------------------------------------------------------------------------
    # property
    # --------------------------------------------------------------------------
    def _get_filepath(self):
        filepath = self._state.get("h5filepath", "")
        return filepath

    filepath = property(_get_filepath)


    # --------------------------------------------------------------------------
    # _h5file
    def _get_h5file(self):
        return self._state.get("h5file", None)

    def _set_h5file(self, h5file):
        self._state["h5file"] = h5file

    _h5file = property(_get_h5file, _set_h5file)


    # --------------------------------------------------------------------------
    # _h5root
    def _get_h5root(self):
        return self._state.get("h5root", None)

    def _set_h5root(self, h5root):
        self._state["h5root"] = h5root

    _h5root = property(_get_h5root, _set_h5root)


    # --------------------------------------------------------------------------
    # public
    # --------------------------------------------------------------------------
    def save_basic(self, run_type, h5_path):
        with h5py.File(h5_path, 'w') as h5:
            self._save_model(h5)
            # self._write_sp_occ(h5, run_type)
            # self._write_sp_TEs(h5)

            # self._write_sp_pop_mulliken_atom(h5, run_type, self.iterations)


    def save_standard(self, run_type, h5_path):
        self.save_basic(run_type, h5_path)

        with h5py.File(h5_path, 'w') as h5:
            self._save_s_matrix(h5)
            self._save_h_matrix(h5)
            self._save_h2_matrix(h5)
            for iteration in range(1, self.iterations + 1):
                self._save_c_matrix(h5, run_type, iteration)
                self._save_density_matrix(h5, run_type, iteration)
                self._save_energy_level(h5, run_type, iteration)

            # if self.scf_converged:
                # self._write_sp_pop_mulliken_atom(h5, run_type, self.iterations)
            # else:
            #     iteration = self.iterations -1
            #     if iteration > 0:
            #         self._write_sp_c_matrix(h5, run_type, iteration)
            #         self._write_sp_density_matrix(h5, run_type, iteration)
            #         # self._write_sp_pop_mulliken_atom(h5, run_type, iteration)


    def save_full(self, run_type, h5_path):
        self.save_standard(run_type, h5_path)


    def save_debug(self, run_type, h5_path):
        self.save_full(run_type, h5_path)
        

    # --------------------------------------------------------------------------
    # I/O
    # --------------------------------------------------------------------------
    def open(self, filepath):
        assert(isinstance(filepath, str))

        if len(filepath) > 0:
            if self.filepath != filepath:
                self._close()
            self._state["h5filepath"] = filepath
            
            is_loadable = False
            if os.path.exists(self.filepath):
                is_loadable = True
            assert(self._h5file is None)
            assert(self._h5root is None)
            self._h5file = h5py.File(self.filepath, 'a')
            self._h5root = self._h5file["/"]
            
            if is_loadable:
                self.load_model()
                TEs = self.get_TEs()
                for itr in range(self.iterations):
                    self.set_total_energy(itr +1, TEs[itr])

        else:
            raise ValueError("cannot open file: path is empty.")


    def _close(self):
        if self._h5file is not None:
            self._h5file.close()
            self._h5file = None
            self._h5root = None


    # --------------------------------------------------------------------------
    # deprecated
    # --------------------------------------------------------------------------
    # def _load_sp(self, h5grp):
    #     logger.debug("deprecated function: _load_sp")
    #     self.load_model(h5grp)
    #     self._load_TEs(h5grp)

        # run_type = "rks"
        # self.get_occ_vector(h5grp, run_type)
        # self.get_s_matrix(h5grp)
        # self.get_c_matrix(h5grp, run_type, self.iterations)

    # --------------------------------------------------------------------------
    # model
    # --------------------------------------------------------------------------
    def load_model(self):
        self._load_model()
    
    def _load_model(self, h5grp=None):
        if h5grp is None:
            h5grp = self._h5root

        self.step_control = h5grp.attrs.get("step_control")
        self.comment = h5grp.attrs.get("comment")

        self.num_of_atoms = h5grp.attrs.get("num_of_atoms")
        self.num_of_AOs = h5grp.attrs.get("num_of_AOs")
        self.num_of_MOs = h5grp.attrs.get("num_of_MOs")

        self.method = h5grp.attrs.get("method")
        self.cut_value = h5grp.attrs.get("cut_value")
        self.CDAM_tau = h5grp.attrs.get("CDAM_tau")
        self.CD_epsilon = h5grp.attrs.get("CD_epsilon")

        self.guess = h5grp.attrs.get("guess")
        self.xc_functional = h5grp.attrs.get("xc_functional")

        self.iterations = h5grp.attrs.get("iterations")
        self.max_iterations = h5grp.attrs.get("max_iterations")
        self._set_scf_converged(h5grp.attrs.get("scf_converged"))

        self.orbital_independence_threshold = h5grp.attrs.get("orbital_independence_threshold")
        self.orbital_independence_threshold_canonical = h5grp.attrs.get("orbital_independence_threshold_canonical")
        self.orbital_independence_threshold_lowdin = h5grp.attrs.get("orbital_independence_threshold_lowdin")
        self.scf_acceleration = h5grp.attrs.get("scf_acceleration")
        self.scf_acceleration_damping_factor = h5grp.attrs.get("scf_acceleration_damping_factor")

        # self._load_sp_model_molecule(h5grp) # TODO
        # self._load_sp_model_basisset(h5grp) # TODO
        self._load_model_control(h5grp)


    def _load_model_control(self, h5grp):
        file_base_name = {}
        for key, value in h5grp["control/file_base_name"].attrs.items():
            file_base_name[key] = value

        self._data['control']['file_base_name'] = file_base_name


    def save_model(self):
        self._save_model()

    def _save_model(self, h5grp=None):
        if h5grp is None:
            h5grp = self._h5root
        
        h5grp.attrs["step_control"] = self.step_control
        h5grp.attrs["comment"] = self.comment

        h5grp.attrs["num_of_atoms"] = self.num_of_atoms
        h5grp.attrs["num_of_AOs"] = self.num_of_AOs
        h5grp.attrs["num_of_MOs"] = self.num_of_MOs

        h5grp.attrs["method"] = self.method
        h5grp.attrs["cut_value"] = self.cut_value
        h5grp.attrs["CDAM_tau"] = self.CDAM_tau
        h5grp.attrs["CD_epsilon"] = self.CD_epsilon
        if self.guess:
            h5grp.attrs["guess"] = self.guess
        h5grp.attrs["xc_functional"] = self.xc_functional

        h5grp.attrs["iterations"] = self.iterations
        h5grp.attrs["max_iterations"] = self.max_iterations
        h5grp.attrs["scf_converged"] = self.scf_converged

        if self.orbital_independence_threshold != None:
            h5grp.attrs["orbital_independence_threshold"] = self.orbital_independence_threshold
        if self.orbital_independence_threshold_canonical != None:
            h5grp.attrs["orbital_independence_threshold_canonical"] = self.orbital_independence_threshold_canonical
        if self.orbital_independence_threshold_lowdin != None:
            h5grp.attrs["orbital_independence_threshold_lowdin"] = self.orbital_independence_threshold_lowdin
        h5grp.attrs["scf_acceleration"] = self.scf_acceleration
        h5grp.attrs["scf_acceleration_damping_damping_factor"] = self.scf_acceleration_damping_damping_factor

        self._save_model_molecule(h5grp)
        self._save_model_basisset(h5grp)
        self._save_model_control(h5grp)


    def _save_model_molecule(self, h5grp):
        symbols = []
        xyzc = []
        for atom_id, atom in self.molecule.atoms():
            symbol = atom.symbol.encode("utf-8")
            xyz = atom.xyz
            charge = atom.charge
            label = atom.name.encode("utf-8")
            if label == None:
                label = ""
            symbols.append(symbol)
            xyzc.append([symbol, label, xyz.x, xyz.y, xyz.z, charge])
        h5ds_molecule = h5grp.create_dataset("molecule/xyzc", data=xyzc)
        h5ds_molecule.attrs["header"] = [
            "symbol", "label", "x", "y", "z", "charge"]


    def _save_model_basisset(self, h5grp):
        h5_atomlabel = []
        h5_cgto = []
        h5_pgto = []
        for atom_label in self.get_basisset_atomlabels():
            basisset = self.get_basisset(atom_label)
            name = basisset.name.encode("utf-8")
            h5_atomlabel.append([name, atom_label.encode("utf-8")])

            for cgto_id, cgto in enumerate(basisset):
                shell_type = cgto.shell_type.encode("utf-8")
                scale_factor = cgto.scale_factor
                h5_cgto.append([name, int(cgto_id), shell_type, scale_factor])

                for pgto_id, pgto in enumerate(cgto):
                    coef = pgto.coef
                    exp = pgto.exp
                    h5_pgto.append(
                        [name, int(cgto_id), int(pgto_id), coef, exp])

        # print(h5_atomlabel)
        h5ds_name = h5grp.create_dataset("basisset/name", data=h5_atomlabel)
        h5ds_name.attrs["header"] = ["name", "atom_label"]

        h5ds_cgto = h5grp.create_dataset("basisset/cgto", data=h5_cgto)
        h5ds_cgto.attrs["header"] = [
            "name", "cgto_id", "shell_type", "scale_factor"]

        h5ds_pgto = h5grp.create_dataset("basisset/pgto", data=h5_pgto)
        h5ds_pgto.attrs["header"] = [
            "name", "cgto_id", "pgto_id", "coef", "exp"]


    def _save_model_control(self, h5grp):
        h5grp_control = h5grp.create_group("control")
        h5grp_control_file_base_name = h5grp_control.create_group(
            "file_base_name")

        control_file_base_name = self._data["control"].get(
            "file_base_name", {})
        for key, value in control_file_base_name.items():
            h5grp_control_file_base_name.attrs[key] = value


    # --------------------------------------------------------------------------
    # occ. vector
    # --------------------------------------------------------------------------
    def save_occ(self, run_type, occ):
        self._save_occ(run_type, occ)

    def _save_occ(self, run_type, occ, h5grp=None):
        assert(isinstance(occ, Vector))
        if h5grp is None:
            h5grp = self._h5root

        h5path = self._get_h5path_occ(run_type)
        self._set_vector(h5grp, h5path, occ)


    def get_occ_vector(self, run_type):
        h5path = self._get_h5path_occ(run_type)
        v = self._get_vector(self._h5root, h5path)
        return v


    def get_HOMO_level(self, run_type):
        """
        HOMOのレベルを返す
        0から始まることに注意。
        """
        answer = None
        occ = self.get_occ_vector(run_type)
        for i in range(len(occ) - 1, 0, -1):
            if occ[i] > 0.0:
                answer = i
                break
        return answer


    # --------------------------------------------------------------------------
    # energy level
    # --------------------------------------------------------------------------
    def save_energy_level(self, run_type, iteration, energy_level):
        self._save_energy_level(run_type, iteration, energy_level)

    def _save_energy_level(self, run_type, iteration, energy_level, h5grp=None):
        assert(isinstance(energy_level, Vector))
        if h5grp is None:
            h5grp = self._h5root

        h5path = self._get_h5path_energy_level(run_type, iteration)
        self._set_vector(h5grp, h5path, energy_level)            

    def get_energy_level(self, run_type, iteration):
        h5path = self._get_h5path_energy_level(run_type, iteration)
        energy_level = self._get_vector(self._h5root, h5path)
        return energy_level

    # --------------------------------------------------------------------------
    # S matrix
    # --------------------------------------------------------------------------
    def get_s_matrix(self, h5grp=None):
        if h5grp is None:
            h5grp = self._h5root
        h5path = self._get_h5path_s()

        S = self._get_matrix(h5grp, h5path)
        return S

    def save_s_matrix(self, S, h5grp=None):
        assert(isinstance(S, SymmetricMatrix))
        if h5grp is None:
            h5grp = self._h5root

        h5path = self._get_h5path_s()
        self._set_matrix(h5grp, h5path, S)


    # --------------------------------------------------------------------------
    # h matrix
    # --------------------------------------------------------------------------
    def get_h_matrix(self, h5grp=None):
        if h5grp is None:
            h5grp = self._h5root
        h5path = self._get_h5path_h()

        h = self._get_matrix(h5grp, h5path)
        return h

    def save_h_matrix(self, h, h5grp=None):
        assert(isinstance(h, SymmetricMatrix))
        if h5grp is None:
            h5grp = self._h5root

        h5path = self._get_h5path_h()
        self._set_matrix(h5grp, h5path, h)


    # --------------------------------------------------------------------------
    # h2 (h for dummy charge) matrix
    # --------------------------------------------------------------------------
    def get_h2_matrix(self, h5grp=None):
        if h5grp is None:
            h5grp = self._h5root
        h5path = self._get_h5path_h2()

        h2 = self._get_matrix(h5grp, h5path)
        return h2

    def save_h2_matrix(self, h2, h5grp=None):
        assert(isinstance(h2, SymmetricMatrix))
        if h5grp is None:
            h5grp = self._h5root
        
        h5path = self._get_h5path_h2()
        self._set_matrix(h5grp, h5path, h2)


    # --------------------------------------------------------------------------
    # C matrix
    # --------------------------------------------------------------------------
    def get_c_matrix(self, run_type, iteration, h5grp=None):
        if h5grp is None:
            h5grp = self._h5root
        h5path = self._get_h5path_c_matrix(run_type, iteration)

        C = self._get_matrix(h5grp, h5path)
        return C


    def save_c_matrix(self, run_type, iteration, C, h5grp=None):
        assert(isinstance(C, Matrix))
        if h5grp is None:
            h5grp = self._h5root

        h5path = self._get_h5path_c_matrix(run_type, iteration)
        self._set_matrix(h5grp, h5path, C, iteration=iteration)


    # --------------------------------------------------------------------------
    # density matrix
    # --------------------------------------------------------------------------
    def get_density_matrix(self, run_type, iteration, h5grp=None):
        if h5grp is None:
            h5grp = self._h5root
        h5path = self._get_h5path_density_matrix(run_type, iteration)

        P = self._get_matrix(h5grp, h5path)
        return P


    def save_density_matrix(self, run_type, iteration, P, h5grp=None):
        assert(isinstance(P, SymmetricMatrix))
        if h5grp is None:
            h5grp = self._h5root

        h5path = self._get_h5path_density_matrix(run_type, iteration)
        self._set_matrix(h5grp, h5path, P, iteration=iteration)


    # --------------------------------------------------------------------------
    # Kohn-SHam matrix
    # --------------------------------------------------------------------------
    def get_f_matrix(self, run_type, iteration, h5grp=None):
        if h5grp is None:
            h5grp = self._h5root
        h5path = self._get_h5path_f_matrix(run_type, iteration)

        F = self._get_matrix(h5grp, h5path)
        return F


    def save_f_matrix(self, run_type, iteration, F, h5grp=None):
        assert(isinstance(F, SymmetricMatrix))
        if h5grp is None:
            h5grp = self._h5root

        h5path = self._get_h5path_f_matrix(run_type, iteration)
        self._set_matrix(h5grp, h5path, F, iteration=iteration)


    # --------------------------------------------------------------------------
    # TotalEnergy
    # --------------------------------------------------------------------------
    def get_TEs(self, h5grp=None):
        if h5grp is None:
            h5grp = self._h5root
        h5path = self._get_h5path_TEs()

        TEs = self._get_vector(h5grp, h5path)
        return TEs

    # def _load_TEs(self, h5grp=None):
    #     if h5grp is None:
    #         h5grp = self._h5root
    #     h5path = self._get_h5path_TEs()
    #     v = self._get_vector(h5grp, h5path)
    #     for i in range(self.iterations):
    #         self.set_total_energy(i + 1, v[i])


    def save_TEs(self, TEs, h5grp=None):
        assert(isinstance(TEs, Vector))
        if h5grp is None:
            h5grp = self._h5root

        h5path = self._get_h5path_TEs()
        self._set_vector(h5grp, h5path, TEs)

    # def save_TEs(self, h5grp=None):
    #     TEs = Vector(self.iterations)
    #     for itr in range(1, self.iterations + 2):
    #         total_energy = self.get_total_energy(itr)
    #         if total_energy != None:
    #             TEs[itr -1] = total_energy

    #     if h5grp is None:
    #         h5grp = self._h5root
    #     h5path = self._get_h5path_TEs()
    #     self._set_vector(h5grp, h5path, TEs)

    # --------------------------------------------------------------------------
    # population
    # --------------------------------------------------------------------------
    def get_pop_mulliken_atom(self, run_type, iteration):
        h5path = self._get_h5path_pop_mulliken_atom(run_type, iteration)
        return self._load_vector(self._h5root, h5path)

    def _write_sp_pop_mulliken_atom(self, h5grp, run_type, iteration, force=False):
        file_path = self.get_pop_mulliken_path(run_type, iteration)

        if iteration > 0:
            if (os.path.exists(file_path) != True) or (force == True):
                c_mat_path = self.get_c_mat_path(run_type, iteration)
                if os.path.exists(c_mat_path):
                    logger.info("calculate pop(Mulliken) data ...")
                    run_pdf("pop-mulliken -i {}".format(iteration))
                else:
                    logger.info(
                        "not found C ({}), give up calculate pop data.".format(c_mat_path))

            vtr = Vector()
            if vtr.load(file_path):
                h5path = self._get_h5path_pop_mulliken_atom(
                    run_type, iteration)
                self._write_vector(h5grp, h5path, vtr)
        else:
            logger.debug(
                "Since the SCF iteration == 0, the charge calculation is not performed.")


    # --------------------------------------------------------------------------
    # Matrix I/O
    # --------------------------------------------------------------------------
    def _get_matrix(self, h5grp, dataset_name):
        h5ds = h5grp.get(dataset_name)

        matrix = None
        #row = h5ds.attrs["row"]
        #col = h5ds.attrs["col"]
        mat_type = h5ds.attrs["type"]
        if mat_type == "GE":
            array = numpy.array(h5ds[:])
            matrix = Matrix(array)
        elif mat_type == "SY":
            array = numpy.array(h5ds[:])
            matrix = SymmetricMatrix(array)

        return matrix

    def _set_matrix(self, h5grp, dataset_name, matrix, **extra_kwds):
        """
        """
        assert(isinstance(dataset_name, str))
        assert(isinstance(matrix, (Matrix, SymmetricMatrix)))

        h5mat = h5grp.create_dataset(dataset_name, data=matrix.data)
        h5mat.attrs["row"] = matrix.rows
        h5mat.attrs["col"] = matrix.cols
        h5mat.attrs["type"] = matrix.type
        for key, value in extra_kwds.items():
            h5mat.attrs[key] = value


    # --------------------------------------------------------------------------
    # Vector I/O
    # --------------------------------------------------------------------------
    def _get_vector(self, h5grp, dataset_name):
        h5ds = h5grp.get(dataset_name)

        vector = None
        #dim = h5ds.attrs["dim"]
        #name = h5ds.attrs["name"]
        if h5ds != None:
            array = numpy.array(h5ds[:])
            vector = Vector(array)

        return vector

    def _set_vector(self, h5grp, dataset_name, vector, **extra_kwds):
        """
        """
        assert(isinstance(dataset_name, str))
        assert(isinstance(vector, Vector))

        h5mat = h5grp.create_dataset(dataset_name, data=vector.data)
        h5mat.attrs["dim"] = len(vector)
        for key, value in extra_kwds.items():
            h5mat.attrs[key] = value


    # --------------------------------------------------------------------------
    # private functions
    # --------------------------------------------------------------------------
    def _get_h5path_occ(self, run_type):
        h5path = "occ/{run_type}".format(run_type=run_type)
        return h5path

    def _get_h5path_s(self):
        h5path = "s"
        return h5path

    def _get_h5path_h(self):
        h5path = "h"
        return h5path

    def _get_h5path_h2(self):
        h5path = "h2"
        return h5path

    def _get_h5path_c_matrix(self, run_type, iteration):
        h5path = "C/{run_type}_{iteration}".format(run_type=run_type,
                                                   iteration=iteration)
        return h5path

    def _get_h5path_density_matrix(self, run_type, iteration):
        h5path = "P/{run_type}_{iteration}".format(run_type=run_type,
                                                   iteration=iteration)
        return h5path

    def _get_h5path_f_matrix(self, run_type, iteration):
        h5path = "F/{run_type}_{iteration}".format(run_type=run_type,
                                                   iteration=iteration)
        return h5path

    def _get_h5path_TEs(self):
        h5path = "TEs"
        return h5path

    def _get_h5path_energy_level(self, run_type, iteration):
        h5path = "energy_level/{run_type}_{iteration}".format(run_type=run_type,
                                                              iteration=iteration)
        return h5path

    def _get_h5path_pop_mulliken_atom(self, run_type, iteration):
        h5path = "pop/mulliken/{run_type}_{iteration}".format(run_type=run_type,
                                                              iteration=iteration)
        return h5path

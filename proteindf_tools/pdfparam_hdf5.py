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
import h5py
import os
import numpy
import traceback

import logging
logger = logging.getLogger(__name__)


class PdfParam_H5(PdfParamObject):
    """ ProteinDF parameter object
    """

    def __init__(self, rhs=None):
        super().__init__(rhs)
        self._h5file = None
        self._h5root = None

        if (isinstance(rhs, str)):
            db_path = rhs
            self.open(db_path)

    def __del__(self):
        if self._h5file != None:
            self._h5file.close()

    # --------------------------------------------------------------------------
    # public
    # --------------------------------------------------------------------------
    def open(self, filepath):
        self._h5file = h5py.File(filepath, 'a')
        self._h5root = self._h5file["/"]
        self._load_sp(self._h5root)

    def save_basic(self, h5_path):
        with h5py.File(h5_path, 'w') as h5:
            self._save_sp_condition(h5)
            self._write_sp_TEs(h5)
            self._write_sp_pop_mulliken_atom(h5, run_type, self.iterations)

    def save_standard(self, h5_path):
        with h5py.File(h5_path, 'w') as h5:
            self._save_sp_condition(h5)
            self._write_sp_TEs(h5)

            run_type = "rks"
            self._write_sp_occ(h5, run_type)
            self._write_sp_s_matrix(h5)
            self._write_sp_h_matrix(h5)
            for i in range(1, self.iterations + 1):
                self._write_sp_energy_level(h5, run_type, i)
            if self.scf_converged:
                self._write_sp_c_matrix(h5, run_type, self.iterations)
                self._write_sp_density_matrix(h5, run_type, self.iterations)
                self._write_sp_pop_mulliken_atom(h5, run_type, self.iterations)

    def save_full(self, h5_path):
        with h5py.File(h5_path, 'w') as h5:
            self._save_sp_condition(h5)
            self._write_sp_TEs(h5)

            run_type = "rks"
            self._write_sp_occ(h5, run_type)
            self._write_sp_s_matrix(h5)
            self._write_sp_h_matrix(h5)
            for i in range(1, self.iterations + 1):
                self._write_sp_energy_level(h5, run_type, i)
                self._write_sp_c_matrix(h5, run_type, i)
                self._write_sp_density_matrix(h5, run_type, i)

    def save_debug(self, h5_path):
        with h5py.File(h5_path, 'w') as h5:
            self._save_sp_condition(h5)
            self._save_sp_TEs(h5)

            run_type = "rks"
            self._write_sp_occ(h5, run_type)
            self._write_sp_s_matrix(h5)
            self._write_sp_h_matrix(h5)
            for i in range(1, self.iterations + 1):
                self._write_sp_energy_level(h5, run_type, i)
                self._write_sp_c_matrix(h5, run_type, i)
                self._write_sp_density_matrix(h5, run_type, i)

    # --------------------------------------------------------------------------
    #
    # --------------------------------------------------------------------------

    def _save_sp_condition(self, h5grp):
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

        self._write_sp_molecule(h5grp)
        self._write_sp_basisset(h5grp)
        self._write_sp_condition_control(h5grp)

    def _write_sp_molecule(self, h5grp):
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

    def _write_sp_basisset(self, h5grp):
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

    def _write_sp_condition_control(self, h5grp):
        h5grp_control = h5grp.create_group("control")
        h5grp_control_file_base_name = h5grp_control.create_group(
            "file_base_name")

        control_file_base_name = self._data["control"].get(
            "file_base_name", {})
        for key, value in control_file_base_name.items():
            h5grp_control_file_base_name.attrs[key] = value

    def _write_sp_s_matrix(self, h5grp):
        file_path = self.get_s_mat_path()
        S = SymmetricMatrix()
        if S.load(file_path):
            h5path = self._get_h5path_s()
            self._write_matrix(h5grp, h5path, S, name="S")

    def _write_sp_h_matrix(self, h5grp):
        file_path = self.get_h_mat_path()
        H = SymmetricMatrix()
        if H.load(file_path):
            h5path = self._get_h5path_h()
            self._write_matrix(h5grp, h5path, H, name="H")

    def _write_sp_occ(self, h5grp, run_type):
        file_path = self.get_occ_path(run_type)
        occ = Vector()
        if occ.load(file_path):
            h5path = self._get_h5path_occ(run_type)
            self._write_vector(h5grp, h5path, occ, name="occ")

    def _write_sp_c_matrix(self, h5grp, run_type, iteration):
        file_path = self.get_c_mat_path(run_type, iteration)
        C = Matrix()
        if C.load(file_path):
            h5path = "C/{}_{}".format(run_type, iteration)
            self._write_matrix(h5grp, h5path, C, name="C", iteration=iteration)

    def _write_sp_energy_level(self, h5grp, run_type, iteration):
        file_path = self.get_energy_level_path(run_type, iteration)
        el = Vector()
        if el.load(file_path):
            h5path = self._get_h5path_energy_level(run_type, iteration)
            self._write_vector(h5grp, h5path, el, name="energy_level")

    def _write_sp_density_matrix(self, h5grp, run_type, iteration):
        file_path = self.get_density_matrix_path(run_type, iteration)
        P = SymmetricMatrix()
        if P.load(file_path):
            h5path = self._get_h5path_density_matrix(run_type, iteration)
            self._write_matrix(h5grp, h5path, P, name="P", iteration=iteration)

    def _write_sp_TEs(self, h5grp):
        vtr_TEs = Vector(self.iterations)
        for itr, TE in self.TEs.items():
            vtr_TEs[itr - 1] = TE

        h5path = self._get_h5path_TEs()
        self._write_vector(h5grp, h5path, vtr_TEs)

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

    def _write_matrix(self, h5grp, dataset_name, matrix, **extra_kwds):
        """
        """
        isinstance(dataset_name, str)

        h5mat = h5grp.create_dataset(dataset_name, data=matrix.data)
        h5mat.attrs["row"] = matrix.rows
        h5mat.attrs["col"] = matrix.cols
        h5mat.attrs["type"] = matrix.type
        for key, value in extra_kwds.items():
            h5mat.attrs[key] = value

    def _write_vector(self, h5grp, dataset_name, vector, **extra_kwds):
        """
        """
        isinstance(dataset_name, str)

        h5mat = h5grp.create_dataset(dataset_name, data=vector.data)
        h5mat.attrs["dim"] = len(vector)
        for key, value in extra_kwds.items():
            h5mat.attrs[key] = value

    def _load_sp(self, h5grp):
        self._load_sp_condition(h5grp)
        self._load_TEs(h5grp)

        run_type = "rks"
        #self.get_occ_vector(h5grp, run_type)
        # self.get_s_matrix(h5grp)
        #self.get_c_matrix(h5grp, run_type, self.iterations)

    def _load_sp_condition(self, h5grp):
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

        self.orbital_independence_threshold = h5grp.attrs.get(
            "orbital_independence_threshold")
        self.orbital_independence_threshold_canonical = h5grp.attrs.get(
            "orbital_independence_threshold_canonical")
        self.orbital_independence_threshold_lowdin = h5grp.attrs.get(
            "orbital_independence_threshold_lowdin")
        self.scf_acceleration = h5grp.attrs.get("scf_acceleration")
        self.scf_acceleration_damping_factor = h5grp.attrs.get(
            "scf_acceleration_damping_factor")

        self._load_sp_condition_control(h5grp)

    def _load_sp_condition_control(self, h5grp):
        file_base_name = {}
        for key, value in h5grp["control/file_base_name"].attrs.items():
            file_base_name[key] = value

        self._data['control']['file_base_name'] = file_base_name

    def _load_TEs(self, h5grp):
        h5path = self._get_h5path_TEs()
        v = self._load_vector(h5grp, h5path)

        TEs = {}
        for i in range(self.iterations):
            TEs[i + 1] = v[i]

        self._set_TEs(TEs)

    def get_occ_vector(self, run_type):
        h5path = self._get_h5path_occ(run_type)
        v = self._load_vector(self._h5root, h5path)
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

    def get_energy_level(self, run_type, iteration):
        h5path = self._get_h5path_energy_level(run_type, iteration)
        energy_level = self._load_vector(self._h5root, h5path)
        return energy_level

    def get_s_matrix(self, h5grp):
        h5path = self._get_h5path_s()
        self._load_matrix(h5grp, h5path)

    def get_c_matrix(self, h5grp, run_type, iteration):
        h5path = self._get_h5path_c_matrix(run_type, iteration)
        self._load_matrix(h5grp, h5path)

    def get_pop_mulliken_atom(self, run_type, iteration):
        h5path = self._get_h5path_pop_mulliken_atom(run_type, iteration)
        return self._load_vector(self._h5root, h5path)

    def _load_matrix(self, h5grp, dataset_name):
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

    def _load_vector(self, h5grp, dataset_name):
        h5ds = h5grp.get(dataset_name)

        vector = None
        #dim = h5ds.attrs["dim"]
        #name = h5ds.attrs["name"]
        if h5ds != None:
            array = numpy.array(h5ds[:])
            vector = Vector(array)

        return vector

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

    def _get_h5path_c_matrix(self, run_type, iteration):
        h5path = "C/{run_type}_{iteration}".format(run_type=run_type,
                                                   iteration=iteration)
        return h5path

    def _get_h5path_energy_level(self, run_type, iteration):
        h5path = "energy_level/{run_type}_{iteration}".format(run_type=run_type,
                                                              iteration=iteration)
        return h5path

    def _get_h5path_density_matrix(self, run_type, iteration):
        h5path = "P/{run_type}_{iteration}".format(run_type=run_type,
                                                   iteration=iteration)
        return h5path

    def _get_h5path_TEs(self):
        h5path = "TEs"
        return h5path

    def _get_h5path_pop_mulliken_atom(self, run_type, iteration):
        h5path = "pop/mulliken/{run_type}_{iteration}".format(run_type=run_type,
                                                              iteration=iteration)
        return h5path

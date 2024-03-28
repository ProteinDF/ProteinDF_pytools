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

from .pdfparam_hdf5 import PdfParam_H5
from .pdfparam import PdfParam
from .pdfcommon import get_default_pdfparam, run_pdf
# from .pdfarchive import PdfArchive
import proteindf_bridge as bridge
import os
import sys
import copy
import math
import logging
logger = logging.getLogger(__name__)


class PdfSim(object):
    """
    """
    def __init__(self, *args, **kwargs):
        self._db_path = kwargs.get('db_path', 'pdfresults.db')

    def setup(self, pdfparam =None, workdir ="."):
        """
        setup to run ProteinDF
        """
        # make fl_Userinput
        if pdfparam is None:
            pdfparam = get_default_pdfparam()
        assert(isinstance(pdfparam, PdfParam))

        if not os.path.exists(workdir):
            os.mkdir(workdir)
        else:
            logger.debug('already exist: {}'.format(workdir))

        input_path = os.path.join(workdir, "fl_Userinput")
        f = open(input_path, 'w')
        f.write(pdfparam.get_inputfile_contents())
        f.close()

        self._make_workdir(workdir)

    def _make_workdir(self, workdir='.'):
        # make sub-directories
        dirs = ['fl_Work']
        for d in dirs:
            path = os.path.join(workdir, d)
            if not os.path.exists(path):
                os.mkdir(path)

    # ------------------------------------------------------------------
    def sp(self, pdfparam, *args, **kwargs):
        """
        calc single-point

        Keyword arguments:
        pdfparam --- PdfParam object represented the calculation condition.
        workdir --- working directory.
        db_path --- ProteinDF ArchiveDB path (default: pdfresults.db)
        dry_run --- if True, the calculation is NOT carried out. Default is False.

        Returns:
        tuple of the number of iterations and the total energy.

        """
        workdir = kwargs.get('workdir', '')
        db_path = kwargs.get('db_path', self._db_path)
        dry_run = kwargs.get('dry_run', False)
        cmd_pdf = kwargs.get('cmd_pdf', 'serial')
        cmd_archive = kwargs.get('cmd_archive', 'archive')

        current_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or ".")

        pdf_workdir = os.path.join('.', workdir)
        if not os.path.exists(pdf_workdir):
            os.mkdir(pdf_workdir)
        else:
            logger.debug('already exist {}'.format(pdf_workdir))

        self.setup(pdfparam, pdf_workdir)
        os.chdir(pdf_workdir)

        itr = None
        total_energy = None
        if not dry_run:
            # exec pdf-serial(default)
            # run_pdf(['serial', '-o', 'pdf.log'])
            run_pdf([cmd_pdf, '-o', 'pdf.log'])

            # exec pdf-archive
            run_pdf(cmd_archive)

            entry = None
            # if cmd_archive == 'archive':
            #     logger.info('read archived file(DB)')
            #     entry = PdfArchive(db_path)
            if cmd_archive == 'archive-h5':
                logger.info('read archived file(HDF5)')
                entry = PdfParam_H5()
                entry.open('pdfresults.h5')
            itr = entry.iterations
            if itr is not None:
                total_energy = entry.get_total_energy(itr)

        os.chdir(current_dir)

        return (itr, total_energy)

    # ------------------------------------------------------------------
    def pop(self, iteration=-1, *args, **kwargs):
        """
        calc population
        """
        workdir = kwargs.get('workdir', '')
        dry_run = kwargs.get('dry_run', False)

        current_dir = os.path.abspath(".")

        pdf_workdir = os.path.join('.', workdir)
        if not os.path.exists(pdf_workdir):
            os.mkdir(pdf_workdir)
        else:
            logger.debug('already exist {}'.format(pdf_workdir))

        os.chdir(pdf_workdir)

        if not dry_run:
            if iteration != -1:
                run_pdf(['pop-mulliken', '-i', iteration])
            else:
                run_pdf(['pop-mulliken'])

        os.chdir(current_dir)

    # ------------------------------------------------------------------

    def opt(self, pdfparam, workdir=".", max_cycle=100):
        current_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or ".")

        alpha = 0.5
        accuracy = 1.0E-3

        opt_cycle = 1
        while True:
            opt_workdir = os.path.join(
                current_dir, "opt_cycle_{}".format(opt_cycle))
            if not os.path.exists(opt_workdir):
                os.mkdir(opt_workdir)
            else:
                self._logger.info('already exist: {}'.format(opt_workdir))
            (grad_mat, rms) = self.numerical_grad(
                pdfparam, opt_workdir, accuracy)

            max_force = 0.0
            for i in range(len(grad_mat)):
                for j in range(3):
                    max_force = max(max_force, math.fabs(grad_mat[i][j]))

            if (max_force < 0.001 and rms < 0.001) or (opt_cycle > max_cycle):
                break

            # update
            i = 0
            for atom_id, atom in pdfparam.molecule.atoms():
                atom.shift_by(alpha * grad_mat[i][0],
                              alpha * grad_mat[i][1],
                              alpha * grad_mat[i][2])

            opt_cycle += 1

        logger.info("opt done")
        logger.info(str(pdfparam.molecule))

    def numerical_grad(self, pdfparam, workdir=".", accuracy=1.0E-3, delta=0.001):
        direction_str = ["x", "y", "z"]

        logger.debug(">>>> start numerical grad")
        molecule = pdfparam.molecule
        num_of_atoms = molecule.get_number_of_all_atoms()

        if pdfparam.convergence_threshold_energy < accuracy * 0.5:
            logger.warning('convergence_thresold_energy={} > 0.5*accuracy(={})'.format(
                pdfparam.convergence_threshold_energy,
                accuracy*0.5))

        grad_mat = [[0.0, 0.0, 0.0] for x in range(num_of_atoms)]
        index = 0
        for atom_id, atom in molecule.atoms():
            for direction in range(3):
                logger.debug("start numerical grad for {}".format(
                    direction_str[direction]))
                h = delta

                delta_TE = 0.0
                while True:
                    atom1 = self._move(atom, direction, -h)
                    atom2 = self._move(atom, direction, +h)

                    pdfparam1 = pdfparam
                    pdfparam1.molecule[atom_id] = atom1
                    workdir1 = os.path.join(
                        workdir, "{}_d{}_h{}_1".format(atom_id, direction, h))
                    (itr, total_energy1) = self.sp(pdfparam1, workdir1)

                    pdfparam2 = pdfparam
                    pdfparam2.molecule[atom_id] = atom2
                    workdir2 = os.path.join(
                        workdir, "{}_d{}_h{}_2".format(atom_id, direction, h))
                    (itr, total_energy2) = self.sp(pdfparam2, workdir2)

                    delta_TE = total_energy2 - total_energy1
                    print("h={: e}, delta_TE={: e}, v={: e}".format(
                        h, delta_TE, delta_TE / (2.0 * h)))
                    if (math.fabs(delta_TE) < accuracy) or (h < 1.0E-2):
                        logger.debug("numerical grad condition satisfied. value={} < {}".format(
                            delta_TE, accuracy))
                        break

                    h *= 0.5
                    sys.stderr.flush()
                    sys.stdout.flush()

                # x,y,zのどれかが終了
                value = delta_TE / (2.0 * h)
                grad_mat[index][direction] = float(value)
                logger.debug("gradient value [{}][{}] = {}".format(atom_id,
                                                                   direction_str[direction],
                                                                   value))
                # self._show_grad_mat(grad_mat)

                sys.stderr.flush()
                sys.stdout.flush()
            index += 1

        logger.info("=== grad (accuracy={}) ===".format(accuracy))
        rms = 0.0
        index = 0
        for atom_id, atom in molecule.atoms():
            logger.info("[{:>8s}] {: 8.5f} {: 8.5f} {: 8.5f}".format(
                atom_id, grad_mat[index][0], grad_mat[index][1], grad_mat[index][2]))
            for i in range(3):
                rms += grad_mat[index][i] * grad_mat[index][i]
            index += 1
        rms /= num_of_atoms * 3
        rms = math.sqrt(rms)
        logger.info("RMS = {}".format(rms))
        logger.info("============\n")

        return (grad_mat, rms)


    def total_energy(self, pdfparam, workdir="."):
        (itr, total_energy) = self.calc_pdf(pdfparam, workdir)

        return (itr, total_energy)

    def _show_grad_mat(self, grad_mat):
        for i in range(len(grad_mat)):
            print("[{}] ({}, {}, {})".format(
                i, grad_mat[i][0], grad_mat[i][1], grad_mat[i][2]))

    def _move(self, atom, direction, delta):
        assert(0 <= direction < 3)
        #answer = copy.deepcopy(atom)
        answer = bridge.Atom(atom)

        if direction == 0:
            answer.xyz.x = float(answer.xyz.x) + delta
        elif direction == 1:
            answer.xyz.y = float(answer.xyz.y) + delta
        elif direction == 2:
            answer.xyz.z = float(answer.xyz.z) + delta
        else:
            logger.critical("program error.")
            exit(1)

        return answer

    # def calc_pdf(self, pdfparam, workdir):
    #     current_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or ".")

    #     pdf_workdir = os.path.join('.', workdir)
    #     if not os.path.exists(pdf_workdir):
    #         os.mkdir(pdf_workdir)
    #     else:
    #         logger.info('already exist {}'.format(pdf_workdir))

    #     self.setup(pdfparam, pdf_workdir)
    #     os.chdir(pdf_workdir)
    #     run_pdf(['-o', 'pdf.log', 'serial'])
    #     run_pdf('archive')

    #     entry = pdf.PdfArchive(self._db_path)
    #     itr = entry.iterations
    #     total_energy = entry.get_total_energy(itr)

    #     os.chdir(current_dir)

    #     return (itr, total_energy)

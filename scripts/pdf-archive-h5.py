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

"""
archive ProteinDF results.
"""

import os
import sys
import argparse
import subprocess
import tempfile
import traceback
import rich.logging
import rich.progress
#import pprint

import proteindf_bridge as bridge
import proteindf_tools as pdf
from proteindf_tools.poputils import compute_mulliken_pop, parse_mulliken_pop

import logging
logger = logging.getLogger(__name__)


def setup_logger(is_debug=False, output_file='pdf-archive-h5.log'):
    handlers = []
    level = logging.INFO
    if is_debug:
        level = logging.DEBUG

    rich_handler = rich.logging.RichHandler(rich_tracebacks=True)
    rich_handler.setLevel(level)
    rich_handler.setFormatter(logging.Formatter("%(message)s"))
    handlers.append(rich_handler)

    if len(output_file) > 0:
        file_handler = logging.FileHandler(output_file)
        file_handler.setLevel(level)
        file_handler.setFormatter(logging.Formatter("%(asctime)s@%(name)s[%(levelname)s] %(message)s"))
        handlers.append(file_handler)

    logging.basicConfig(level=logging.NOTSET,
                        handlers=handlers)



LEVELS = ('basic', 'standard', 'full', 'debug')
DEFAULT_LEVEL = 'standard'


def _last_iterations(pdfparam):
    """Iteration list holding only the converged (last) iteration."""
    if pdfparam.iterations <= 0:
        return []
    return [pdfparam.iterations]


def archive(pdfparam, entry, calc_pop=False, level=DEFAULT_LEVEL):
    if level not in LEVELS:
        raise ValueError("unknown archive level: %s" % level)
    # None means all iterations
    itrs = None if level in ('full', 'debug') else _last_iterations(pdfparam)

    logger.info("save model")
    entry.save_model()

    logger.info("save TEs")
    archive_TEs(pdfparam, entry)

    logger.info("save S")
    archive_s(pdfparam, entry)

    logger.info("save h")
    archive_h(pdfparam, entry)
    archive_h2(pdfparam, entry)

    if level != 'basic':
        logger.info("save gradient")
        archive_gradient(pdfparam, entry)

    for run_type in pdfparam.run_types():
        logger.info('run_type: %s' % run_type)

        logger.info("save occ")
        archive_occ(pdfparam, run_type, entry)

        logger.info("save P")
        archive_p(run_type, pdfparam, entry, iterations=itrs)

        if level != 'basic':
            logger.info("save input P")
            archive_input_p(run_type, pdfparam, entry)

        if level == 'debug':
            logger.info("save pop (Mulliken) for all iterations")
            for itr in range(1, pdfparam.iterations + 1):
                archive_pop_mulliken(pdfparam, run_type, entry, iteration=itr)
        elif calc_pop:
            logger.info("save pop (Mulliken)")
            archive_pop_mulliken(pdfparam, run_type, entry)

    for run_type in pdfparam.orbital_run_types():
        logger.info('orbital run_type: %s' % run_type)

        logger.info("save energy level")
        archive_energy_level(pdfparam, run_type, entry, iterations=itrs)

        logger.info("save C")
        archive_c(run_type, pdfparam, entry, iterations=itrs)

        if level in ('full', 'debug'):
            logger.info("save F")
            archive_f(run_type, pdfparam, entry)


def archive_pop_mulliken(pdfparam, run_type, entry, iteration=None):
    if pdfparam.iterations <= 0:
        logger.debug("Since the SCF iteration == 0, the charge calculation is not performed.")
        return

    if iteration is None:
        last_itr = pdfparam.iterations if pdfparam.scf_converged else pdfparam.iterations - 1
    else:
        last_itr = iteration
    if last_itr <= 0:
        return

    S = entry.get_s_matrix()
    if S is None:
        logger.warning("failed to load S from archive for Mulliken population")
        return

    # P of this iteration may not be archived at the selected level
    P = pdf.SymmetricMatrix()
    if not P.load(pdfparam.get_density_matrix_path(run_type, last_itr)):
        logger.warning(
            "failed to load P (%s, itr %d) for Mulliken population"
            % (run_type, last_itr)
        )
        return

    logger.info("calculate pop(Mulliken) data for iteration %d ..." % last_itr)
    try:
        pop = compute_mulliken_pop(entry, P, S)
    except Exception as e:
        logger.warning("failed to run pop-mulliken for iteration %d: %s" % (last_itr, e))
        return

    entry.save_pop_mulliken_atom(run_type, last_itr, pop)


def archive_occ(pdfparam, run_type, entry):
    # occupancy
    file_path = pdfparam.get_occ_path(run_type)
    occ = pdf.Vector()
    if occ.load(file_path):
        entry.save_occ(run_type, occ)
    else:
        logger.warning('failed to load occ (%s): %s' % (run_type, file_path))


def archive_TEs(pdfparam, entry):
    iterations = pdfparam.iterations
    TEs = pdf.Vector(iterations)
    for itr in range(1, iterations+1):
        itr = int(itr)
        TEs[itr - 1] = pdfparam.get_total_energy(itr)
    
    entry.save_TEs(TEs)


def archive_gradient(pdfparam, entry):
    logger.debug("save gradient")
    grad_mat = None

    num_atoms = pdfparam.num_of_atoms
    if num_atoms and num_atoms > 0:
        rows = []
        for atom_index in range(num_atoms):
            try:
                v = pdfparam.get_gradient(atom_index)
            except (AssertionError, TypeError):
                v = None
            if v is None or len(v) != 3:
                rows = []
                break
            rows.append(v)
        if len(rows) == num_atoms:
            grad_mat = pdf.Matrix(num_atoms, 3)
            for r, row in enumerate(rows):
                for c in range(3):
                    grad_mat.set(r, c, float(row[c]))

    if grad_mat is None:
        candidate_paths = ["force.mtx"]
        if getattr(pdfparam, "work_path", None):
            work_force = os.path.join(pdfparam.work_path, "force.mtx")
            if work_force not in candidate_paths:
                candidate_paths.append(work_force)

        for filepath in candidate_paths:
            if os.path.isfile(filepath):
                force_mat = pdf.Matrix()
                if force_mat.load(filepath):
                    # force.mtx holds ProteinDF's force; gradient = -force
                    grad_mat = pdf.Matrix(force_mat.rows, force_mat.cols)
                    for r in range(force_mat.rows):
                        for c in range(force_mat.cols):
                            grad_mat.set(r, c, -force_mat.get(r, c))
                else:
                    logger.warning("failed to load force matrix: %s" % filepath)
                break

    if grad_mat is None:
        logger.warning("force not found; gradient is not saved")
        return

    entry.save_gradient(grad_mat)


def _iteration_range(pdfparam, iterations):
    if iterations is None:
        return range(1, pdfparam.iterations + 1)
    return iterations


def archive_energy_level(pdfparam, run_type, entry, iterations=None):
    for itr in _iteration_range(pdfparam, iterations):
        e_level_path = pdfparam.get_energy_level_path(run_type, itr)
        e_level = pdf.Vector()
        if e_level.load(e_level_path):
            entry.save_energy_level(run_type, itr, e_level)
        else:
            logger.warning('failed to load energy level (%s): %s' % (run_type, e_level_path))


def archive_s(pdfparam, entry):
    logger.debug("save S")
    filepath = pdfparam.get_s_mat_path()
    S = pdf.SymmetricMatrix()
    S.load(filepath)
    entry.save_s_matrix(S)


def archive_h(pdfparam, entry):
    logger.debug("save h")
    filepath = pdfparam.get_h_mat_path()
    h = pdf.SymmetricMatrix()
    h.load(filepath)
    entry.save_h_matrix(h)


def archive_h2(pdfparam, entry):
    logger.debug("save h2")
    filepath = pdfparam.get_h2_mat_path()
    h2 = pdf.SymmetricMatrix()
    h2.load(filepath)
    entry.save_h2_matrix(h2)


def archive_c(run_type, pdfparam, entry, iterations=None):
    for itr in rich.progress.track(_iteration_range(pdfparam, iterations)):
        logger.debug("save c({})".format(itr))
        filepath = pdfparam.get_c_mat_path(run_type, itr)
        c = pdf.Matrix()
        if c.load(filepath):
            entry.save_c_matrix(run_type, itr, c)
        else:
            logger.warning('failed to load C (%s): %s' % (run_type, filepath))


def archive_p(run_type, pdfparam, entry, iterations=None):
    for itr in rich.progress.track(_iteration_range(pdfparam, iterations)):
        logger.debug("save p({})".format(itr))
        filepath = pdfparam.get_density_matrix_path(run_type, itr)
        p = pdf.SymmetricMatrix()
        if p.load(filepath):
            entry.save_density_matrix(run_type, itr, p)
        else:
            logger.warning('failed to load P (%s): %s' % (run_type, filepath))


def archive_input_p(run_type, pdfparam, entry):
    if pdfparam.iterations <= 0:
        return
    itr = 1
    logger.debug("save input p({})".format(itr))
    filepath = pdfparam.get_input_density_matrix_path(run_type, itr)
    if not os.path.exists(filepath):
        logger.warning('input P not found (%s): %s' % (run_type, filepath))
        return
    p = pdf.SymmetricMatrix()
    if p.load(filepath):
        entry.save_input_density_matrix(run_type, itr, p)
    else:
        logger.warning('failed to load input P (%s): %s' % (run_type, filepath))


def archive_f(run_type, pdfparam, entry, iterations=None):
    for itr in rich.progress.track(_iteration_range(pdfparam, iterations)):
        logger.debug("save F({})".format(itr))
        filepath = pdfparam.get_f_mat_path(run_type, itr)
        f = pdf.SymmetricMatrix()
        if f.load(filepath):
            entry.save_f_matrix(run_type, itr, f)
        else:
            logger.warning('failed to load F (%s): %s' % (run_type, filepath))


def main():
    parser = argparse.ArgumentParser(description='archive ProteinDF results')
    parser.add_argument('-p', '--pdfparam',
                        dest='pdfparam_path',
                        default='pdfparam.mpac',
                        help='ProteinDF parameter file')
    parser.add_argument('-o', '--output',
                        dest='output',
                        nargs='?',
                        default='pdfresults.h5',
                        help='ProteinDF results file')

    parser.add_argument('-f', '--force',
                        action='store_true',
                        default=False,
                        help='overwrite output file')
    parser.add_argument('--pop',
                        action='store_true',
                        default=False,
                        help='calculate Mulliken population for converged solution')
    parser.add_argument('-l', '--level',
                        choices=LEVELS,
                        default=DEFAULT_LEVEL,
                        help='archive level (default: %(default)s)')
    # parser.add_argument('-v', '--verbose',
    #                     action='store_true',
    #                     default=False)
    parser.add_argument('--debug',
                        action='store_true',
                        default=False,
                        help='debug output')
    args = parser.parse_args()

    # setting from command-line
    pdfparam_path = args.pdfparam_path
    output = args.output
    force_overwrite = args.force
    calc_pop = args.pop
    level = args.level
    # verbose = args.verbose
    is_debug = args.debug
    #is_little_endian = True

    # setup logger
    setup_logger(is_debug)

    # check overwrite
    if os.path.exists(output):
        if force_overwrite:
            os.remove(output)
        else:
            logger.error("output file already exists: {}".format(output))
            return
            
    # do archive
    try:
        # read ProteinDF parameters
        pdfparam = pdf.load_pdfparam(pdfparam_path)

        # setup DB
        logger.debug("setup DB")
        entry = pdf.PdfParam_H5(pdfparam)
        entry.open(output)

        logger.debug("call archive")
        archive(pdfparam, entry, calc_pop=calc_pop, level=level)

    except:
        logger.exception("failed to archive")


if __name__ == '__main__':
    main()

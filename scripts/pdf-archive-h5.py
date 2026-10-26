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
import contextlib
import io
import tempfile
import traceback
import rich.logging
import rich.progress
#import pprint

import proteindf_bridge as bridge
import proteindf_tools as pdf

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



def archive(pdfparam, entry, calc_pop=False):
    logger.info("save model")
    entry.save_model()

    logger.info("save TEs")
    archive_TEs(pdfparam, entry)

    logger.info("save S")
    archive_s(pdfparam, entry)

    logger.info("save h")
    archive_h(pdfparam, entry)
    archive_h2(pdfparam, entry)

    for run_type in pdfparam.run_types():
        logger.info('run_type: %s' % run_type)

        logger.info("save occ")
        archive_occ(pdfparam, run_type, entry)

        logger.info("save P")
        archive_p(run_type, pdfparam, entry)

        logger.info("save input P")
        archive_input_p(run_type, pdfparam, entry)

        if calc_pop:
            logger.info("save pop (Mulliken)")
            archive_pop_mulliken(pdfparam, run_type, entry)

    for run_type in pdfparam.orbital_run_types():
        logger.info('orbital run_type: %s' % run_type)

        logger.info("save energy level")
        archive_energy_level(pdfparam, run_type, entry)

        logger.info("save C")
        archive_c(run_type, pdfparam, entry)

        logger.info("save F")
        archive_f(run_type, pdfparam, entry)


def parse_mulliken_pop(output):
    """Parse Gross atom population table from pdf-pop-mulliken output.

    Parameters
    ----------
    output : str
        Standard output of pdf-pop-mulliken.

    Returns
    -------
    pdf.Vector or None
        Vector of atomic charges (net charges).
    """
    lines = output.splitlines()
    in_table = False
    charges = []
    for line in lines:
        stripped = line.strip()
        if "Gross atom population" in stripped:
            in_table = True
            continue
        if in_table:
            if not stripped:
                if charges:
                    break
                continue
            if stripped.startswith("index") or stripped.startswith("Mulliken"):
                if stripped.startswith("Mulliken") and charges:
                    break
                continue
            tokens = stripped.split()
            if len(tokens) >= 4:
                try:
                    int(tokens[0])
                    charges.append(float(tokens[3]))
                except ValueError:
                    if charges:
                        break
            elif len(tokens) == 3:
                try:
                    int(tokens[0])
                    charges.append(float(tokens[2]))
                except ValueError:
                    if charges:
                        break
    if not charges:
        return None
    return pdf.Vector(charges)


def archive_pop_mulliken(pdfparam, run_type, entry):
    if pdfparam.iterations <= 0:
        logger.debug("Since the SCF iteration == 0, the charge calculation is not performed.")
        return

    last_itr = pdfparam.iterations if pdfparam.scf_converged else pdfparam.iterations - 1
    if last_itr <= 0:
        return

    S = entry.get_s_matrix()
    if S is None:
        logger.warning("failed to load S from archive for Mulliken population")
        return

    P = entry.get_density_matrix(run_type, last_itr)
    if P is None:
        logger.warning(
            "failed to load P (%s, itr %d) from archive for Mulliken population"
            % (run_type, last_itr)
        )
        return

    with tempfile.TemporaryDirectory() as tmpdir:
        mpac_path = os.path.join(tmpdir, "pdfparam.mpac")
        s_path = os.path.join(tmpdir, "S.mat")
        p_path = os.path.join(tmpdir, "P.mat")

        entry.export_pdfparam_mpac(mpac_path)
        S.save(s_path)
        P.save(p_path)

        logger.info("calculate pop(Mulliken) data for iteration %d ..." % last_itr)
        buf = io.StringIO()
        try:
            with contextlib.redirect_stdout(buf):
                pdf.run_pdf("pop-mulliken -p {} -P {} -S {}".format(mpac_path, p_path, s_path))
        except Exception as e:
            logger.warning("failed to run pop-mulliken for iteration %d: %s" % (last_itr, e))
            return

        output = buf.getvalue()
        pop = parse_mulliken_pop(output)
        if pop is not None:
            entry.save_pop_mulliken_atom(run_type, last_itr, pop)
        else:
            logger.warning(
                "failed to parse pop mulliken from output (%s, itr %d)" % (run_type, last_itr)
            )


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


def archive_energy_level(pdfparam, run_type, entry):
    iterations = pdfparam.iterations
    energy_levels = pdf.Vector(iterations)
    for itr in range(1, iterations+1):
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


def archive_c(run_type, pdfparam, entry):
    iterations = pdfparam.iterations
    for itr in rich.progress.track(range(1, iterations+1)):
        logger.debug("save c({})".format(itr))
        filepath = pdfparam.get_c_mat_path(run_type, itr)
        c = pdf.Matrix()
        if c.load(filepath):
            entry.save_c_matrix(run_type, itr, c)
        else:
            logger.warning('failed to load C (%s): %s' % (run_type, filepath))


def archive_p(run_type, pdfparam, entry):
    iterations = pdfparam.iterations
    for itr in rich.progress.track(range(1, iterations+1)):
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


def archive_f(run_type, pdfparam, entry):
    iterations = pdfparam.iterations
    for itr in rich.progress.track(range(1, iterations+1)):
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
        archive(pdfparam, entry, calc_pop=calc_pop)

    except:
        logger.exception("failed to archive")


if __name__ == '__main__':
    main()

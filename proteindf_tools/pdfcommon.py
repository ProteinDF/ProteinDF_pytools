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

from .process import Process
from .pdfparam import PdfParam
from .basis2 import Basis2

import proteindf_bridge as bridge

import traceback
import io
import sys
import os
import shlex
import tempfile

import logging

logger = logging.getLogger(__name__)

epsilon = 1.0e-10  # machine epsilon
error = 1.0e-4  # allowable error


def pdf_home():
    """
    return PDF_HOME environment parameter value
    """
    answer = os.environ.get("PDF_HOME", "")
    return answer


def get_default_pdfparam(verbose=False):
    """
    Return the default pdfparam.
    """
    # make temp dir & filepath
    tempfile_fd, tempfile_path = tempfile.mkstemp()
    os.close(tempfile_fd)

    # read the initialization info from the temp file
    args = ["init-param"]
    if verbose:
        args.append("-v")
    args.append("-o")
    args.append(tempfile_path)
    run_pdf(args)
    tempdata = bridge.load_msgpack(tempfile_path)

    # remove temp
    os.remove(tempfile_path)

    pdfparam = PdfParam(tempdata)
    pdfparam.step_control = "create integral guess scf"

    # pdfparam.guess = 'harris'
    # pdfparam.orbital_independence_threshold = 0.007
    # pdfparam.orbital_independence_threshold_canonical = 0.007
    # pdfparam.orbital_independence_threshold_lowdin = 0.007
    # pdfparam.scf_acceleration = 'damping'
    # pdfparam.scf_acceleration_damping_factor = 0.85
    # pdfparam.convergence_threshold_energy = 1.0E-4
    # pdfparam.convergence_threshold = 1.0E-3
    # pdfparam.scf_acceleration_damping_damping_type = 'density_matrix'
    # pdfparam.xc_functional = "b3lyp"
    # pdfparam.j_engine = "CD"
    # pdfparam.k_engine = "CD"
    # pdfparam.xc_engine = "grid"
    # pdfparam.gridfree_orthogonalize_method = "canonical"

    return pdfparam


def set_basisset(
    pdfparam,
    basisset_name_ao="DZVP2",
    basisset_name_rij="DZVP2",
    basisset_name_rixc="DZVP2",
    basisset_name_gridfree="cc-pVDZ-SP",
):
    """
    Set the basisset on pdfparam.
    """
    assert isinstance(pdfparam, PdfParam)

    basis2 = Basis2()
    atoms = ["C", "H", "N", "O", "S"]
    basisset = {}

    if basisset_name_ao == basisset_name_gridfree:
        pdfparam.gridfree_dedicated_basis = False
    else:
        pdfparam.gridfree_dedicated_basis = True

    for atom in atoms:
        basisset_ao = basis2.get_basisset("O-{}.{}".format(basisset_name_ao, atom))
        basisset_j = basis2.get_basisset_j("A-{}.{}".format(basisset_name_rij, atom))
        basisset_xc = basis2.get_basisset_xc("A-{}.{}".format(basisset_name_rixc, atom))
        basisset_gf = basis2.get_basisset("O-{}.{}".format(basisset_name_gridfree, atom))

        pdfparam.set_basisset(atom, basisset_ao)
        pdfparam.set_basisset_j(atom, basisset_j)
        pdfparam.set_basisset_xc(atom, basisset_xc)
        pdfparam.set_basisset_gridfree(atom, basisset_gf)

    return pdfparam


def run_pdf(subcmd):
    """
    run ProteinDF command
    """
    logger.debug("run_pdf({})".format(subcmd))

    try:
        if isinstance(subcmd, list):
            subcmd_tmp = [str(x) for x in subcmd]
            subcmd = " ".join(subcmd_tmp)
    except:
        print(subcmd)
        raise

    # cmd = os.path.join(pdf_home(), "bin", "pdf") + " " + subcmd
    cmd = "pdf" + " " + subcmd
    logger.debug("run: {0}".format(cmd))

    p = Process()
    return_code = p.cmd(cmd).commit()

    logger.debug("return code={}".format(return_code))
    if return_code != 0:
        sys.stderr.write("Failed to execute command: %s" % cmd)
        sys.stderr.write("return code = {}".format(return_code))
        logger.critical("Failed to execute command: %s" % cmd)
        logger.critical("return code = {}".format(return_code))
        raise RuntimeError(
            "Failed to execute command: {} (return code={})".format(cmd, return_code)
        )


def mpac2py(path):
    """
    load message pack binary file to python dictionary data
    """
    assert isinstance(path, str) == True

    data = None
    with open(path, "rb") as f:
        contents = f.read()
        unpacked_data = msgpack.unpackb(contents)
        data = bridge.StrUtils.to_unicode_dict(unpacked_data)

    return data


def load_pdfparam(pdfparam_path="pdfparam.mpac"):
    data = bridge.load_msgpack(pdfparam_path)
    param = PdfParam(data)

    return param


def save_pdfparam(pdfparam_data, pdfparam_path):
    assert isinstance(pdfparam_path, str)
    raw_data = pdfparam_data.get_raw_data()
    bridge.save_msgpack(raw_data, pdfparam_path)


def open_archive(path):
    """Open an archive file (.h5/.hdf5 as PdfParam_H5, otherwise error)."""
    path_str = os.fspath(path)
    ext = os.path.splitext(path_str)[1].lower()
    if ext in (".h5", ".hdf5"):
        from .pdfparam_hdf5 import PdfParam_H5

        entry = PdfParam_H5()
        entry.open(path_str)
        return entry
    else:
        raise ValueError(
            f"SQLite 形式は廃止されました。pdf-archive-h5.py で pdfresults.h5 を作り直してください: {path_str}"
        )

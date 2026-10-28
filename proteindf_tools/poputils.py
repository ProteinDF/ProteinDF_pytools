#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import os
import subprocess
import tempfile

from .pdfcommon import load_pdfparam, open_archive, save_pdfparam
from .vector import Vector
import proteindf_bridge as bridge

import logging
logger = logging.getLogger(__name__)


def parse_mulliken_pop(output):
    """Parse Gross atom population table from pdf-pop-mulliken output.

    Parameters
    ----------
    output : str
        Standard output of pdf-pop-mulliken.

    Returns
    -------
    Vector or None
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
    return Vector(charges)


def compute_mulliken_pop(pdfparam, P, S):
    """Compute Mulliken population using pdf-pop-mulliken.

    Parameters
    ----------
    pdfparam : PdfParamObject or PdfParam_H5
        ProteinDF parameter object.
    P : Matrix or SymmetricMatrix
        Density matrix.
    S : SymmetricMatrix
        Overlap matrix S.

    Returns
    -------
    Vector
        Vector of atomic Mulliken charges.
    """
    with tempfile.TemporaryDirectory() as tmpdir:
        mpac_path = os.path.join(tmpdir, "pdfparam.mpac")
        s_path = os.path.join(tmpdir, "S.mat")
        p_path = os.path.join(tmpdir, "P.mat")

        save_pdfparam(pdfparam, mpac_path)
        S.save(s_path)
        P.save(p_path)

        cmd = ["pdf", "pop-mulliken", "-p", mpac_path, "-P", p_path, "-S", s_path]
        try:
            res = subprocess.run(cmd, capture_output=True, text=True)
            if res.returncode != 0:
                raise RuntimeError(
                    "pop-mulliken failed with exit code {}: {}".format(
                        res.returncode, res.stderr
                    )
                )
        except Exception as e:
            raise RuntimeError("failed to run pop-mulliken: {}".format(e))

        pop = parse_mulliken_pop(res.stdout)
        if pop is None:
            raise RuntimeError("failed to parse Mulliken population from pop-mulliken output")
        return pop


class PopUtils(object):
    '''Population Utilities
    '''
    @staticmethod
    def get_atomlist_by_db(db_path):
        entry = open_archive(db_path)
        molecule = entry.get_molecule()
        return PopUtils._get_atomlist(molecule)

    @staticmethod
    def get_atomlist_by_param(param_path):
        pdfparam = load_pdfparam(param_path)
        molecule = pdfparam.molecule

        return PopUtils._get_atomlist(molecule)

    @staticmethod
    def _get_atomlist(molecule):
        nuc_chargeX = molecule.nuclei_charge
        nuc_charge = molecule.real_nuclei_charge
        charge = nuc_charge - nuc_chargeX

        atoms = molecule.get_atom_list()
        atomlist = []
        for atom in atoms:
            if atom.symbol != 'X':
                atomlist.append(atom)

        return (atomlist, charge)

    @staticmethod
    def set_atomlist(atom_charges, atomlist):
        '''substitute atom_charges for atomlist
        '''
        atom_charges = list(atom_charges)
        atomlist = list(atomlist)

        q_total = 0.0
        for i in range(len(atom_charges)):
            q = atom_charges[i]
            atomlist[i].charge = q
            q_total += q
            print(atomlist[i])
        if abs(q_total) < 1e-12:
            q_total = 0.0
        print('total charges: {: .3f}'.format(q_total))


class PopEspUtils(PopUtils):
    '''Population Utilities for ESP
    '''

    def get_RRMS(self, mpac_path, atoms):
        grids, ESPs = self._load_ESPs(mpac_path)
        rrms = self._calcRRMS(grids, ESPs, atoms)

        return rrms

    def _load_ESPs(self, mpac_path):
        logger.debug('load: {}'.format(mpac_path))
        data = bridge.load_msgpack(mpac_path)
        if data != None:
            grids = []
            for p in data['grids']:
                pos = bridge.Position(p[0], p[1], p[2])
                pos *= (1.0 / 0.5291772108)  # Angstroam to a.u.
                grids.append(pos)

            ESPs = []
            for v in data['ESP']:
                ESPs.append(v)

            assert(len(grids) == len(ESPs))
        return (grids, ESPs)

    def _calcRRMS(self, grids, ESPs, atoms):
        sum_delta2 = 0.0
        sum_v2 = 0.0
        num_of_grids = len(grids)
        logger.debug('# of grids: {}'.format(num_of_grids))
        assert(len(ESPs) == num_of_grids)
        for i in range(num_of_grids):
            estimate_esp = self._calc_esp(grids[i], atoms)
            exact_esp = ESPs[i]
            delta = estimate_esp - exact_esp
            delta2 = delta * delta
            sum_delta2 += delta2
            sum_v2 += exact_esp * exact_esp
        rrms2 = sum_delta2 / sum_v2
        rrms = math.sqrt(rrms2)
        logger.debug('delta2={: .3f} sum_v2={: .3f} RRMS2={: .3f} RRMS={: .3f}'.format(
            sum_delta2, sum_v2, rrms2, rrms))

        return rrms

    def _calc_esp(self, pos, atoms):
        esp = 0.0
        num_of_atoms = len(atoms)
        for i in range(num_of_atoms):
            d = pos.distance_from(atoms[i].xyz)
            esp += atoms[i].charge / d

        return esp

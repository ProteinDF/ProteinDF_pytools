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
make msgpack files from the results of ProteinDF calculation for Harris functional DB.
"""

import os
import sys
import argparse
import array

import proteindf_bridge as bridge
import proteindf_tools as pdf


def main():
    # initialize
    parser = argparse.ArgumentParser(description='register Harris DB')
    parser.add_argument('density_matrix',
                        nargs=1,
                        help='density matrix')
    parser.add_argument('atomlabel',
                        nargs=1,
                        help='target atom label')
    parser.add_argument('-o', '--output',
                        nargs=1,
                        dest='db_path',
                        default='harris.mpac',
                        help='update database file path')
    parser.add_argument('-p', '--param',
                        dest='param',
                        default='pdfparam.mpac',
                        help='ProteinDF parameter path')
    args = parser.parse_args()

    # set up parameters
    atomlabel = args.atomlabel[0]
    harris_db_path = args.db_path[0]
    density_matrix_path = args.density_matrix[0]

    atom = ''
    label = ''
    atomlabel_parts = atomlabel.split('@')
    atom = atomlabel_parts[0]
    if len(atomlabel_parts) > 1:
        label = atomlabel_parts[1]
    if len(label) > 0:
        print('target atom %s with label %s' % (atom, label))
    else:
        print('target atom %s' % (atom))

    # read pdfparam
    pdfparam = pdf.load_pdfparam(args.param)
    orb_info = pdf.OrbInfo(pdfparam)

    # prepare harris DB
    harris_db = {}
    if (os.path.exists(harris_db_path)):
        print("load database: %s" % (harris_db_path))
        harris_db = bridge.load_msgpack(harris_db_path)

    # set basis_set
    basisset = pdfparam.get_basisset(atomlabel)
    harris_db.setdefault('basis_sets', {})
    harris_db['basis_sets'][atom] = basisset.get_raw_data()

    # coord
    #mol = pdfparam.molecule
    # for obj in mol:
    #    obj_atomlabel = '%s' % (obj.symbol)
    #    if len(obj.name) > 0:
    #        obj_atomlabel += '@%s' % (obj.name)
    #
    #    if obj_atomlabel == atomlabel:
    #        harris_db['harris'][atom] = obj
    #        break

    # density matrix
    print("density matrix: %s" % (density_matrix_path))
    mat = pdf.SymmetricMatrix()
    mat.load(density_matrix_path)

    # 密度行列の切り出し処理
    target_orbs = []
    num_of_orbitals = mat.rows
    assert(num_of_orbitals == mat.cols)
    assert(num_of_orbitals == orb_info.get_num_of_orbitals())
    for i in range(num_of_orbitals):
        obj = orb_info.get_atom(i)
        obj_atomlabel = '%s' % (obj.symbol)
        if len(obj.name) > 0:
            obj_atomlabel += '@%s' % (obj.name)

        if obj_atomlabel == atomlabel:
            target_orbs.append(i)

    num_of_target_orbs = len(target_orbs)
    atom_dens_mat = pdf.SymmetricMatrix(num_of_target_orbs)
    for i in range(num_of_target_orbs):
        for j in range(i + 1):
            atom_dens_mat.set(i, j, mat.get(target_orbs[i],
                                            target_orbs[j]))

    # 登録
    harris_db.setdefault('density_matrix', {})
    harris_db['density_matrix'][atom] = atom_dens_mat.get_raw_data()

    # output
    bridge.save_msgpack(harris_db, harris_db_path)


if __name__ == '__main__':
    main()

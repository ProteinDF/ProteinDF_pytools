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

import array
import logging
import math
import sqlite3

import proteindf_bridge as bridge
from .pdfparam import PdfParam
from . import pdfcommon as pdf


class PdfArchive_Sqlite3(object):
    """
    ProteinDFの結果をDBに格納するクラス
    """
    _format_version = 20121230  # DB format version

    def __init__(self,
                 db_path='pdfentry.db',
                 pdfparam=None):
        """
        オブジェクトを作成する
        """
        self._logger = logging.getLogger(__name__)

        self._pdf_id = "0"
        self._db = bridge.DbManager(db_path)
        user_version = self._db.get_user_version()
        if user_version > PdfArchive._format_version:
            self._logger.warning('DB format version is later: code=%d < db=%d' % (
                PdfArchive._format_version, user_version))
        self._create_tables()
        self._read_db()

        if isinstance(pdfparam, PdfParam):
            pdf_id = pdfparam.digest()
            self._set_pdfparam_conditions(pdfparam)
            self._set_pdfparam_coordinates(pdfparam)
            self._set_pdfparam_basisset(pdfparam)
            self._set_pdfparam_total_energies(pdfparam)
            self._set_pdfparam_gradient(pdfparam)

    # private ------------------------------------------------------------------
    def _create_tables(self):
        self._create_table_conditions()
        self._create_table_vectors()
        self._create_table_matrices()
        self._create_table_gradient()

    def _create_table_conditions(self):
        """
        pdf_id は多重登録を防ぐため
        """
        table_name = 'conditions'
        if not self._db.has_table(table_name):
            self._db.create_table(table_name,
                                  ['pdf_id',
                                   'comment',
                                   'num_of_atoms',
                                   'num_of_AOs',
                                   'num_of_MOs',
                                   'method',
                                   'guess',
                                   'xc_functional',
                                   'iterations',
                                   'max_iterations',
                                   'scf_converged'],
                                  ['pdf_id'])

    def _create_table_vectors(self):
        """
        ベクトル情報のDBを作成

        dtype(data type): occ, elevelなど
        """
        table_name = 'vectors'
        if not self._db.has_table(table_name):
            self._db.create_table(table_name,
                                  ['dtype',
                                   'runtype',
                                   'iteration',
                                   'size',
                                   'data'],
                                  ['dtype',
                                   'runtype',
                                   'iteration'])

    def _create_table_matrices(self):
        """
        行列情報のDBを作成

        dtype(data type): C, Pなど
        """
        table_name = 'matrices'
        if not self._db.has_table(table_name):
            self._db.create_table(table_name,
                                  ['dtype',
                                   'runtype',
                                   'iteration',
                                   'rows',
                                   'cols',
                                   'format',
                                   'data'],
                                  ['dtype',
                                   'runtype',
                                   'iteration'])

    def _create_table_gradient(self):
        """
        gradient情報のDBを作成
        """
        table_name = 'gradient'
        if not self._db.has_table(table_name):
            self._db.create_table(table_name,
                                  ['atom_id',
                                   'x',
                                   'y',
                                   'z'],
                                  ['atom_id'])

    def _read_db(self):
        # set pdf_id
        table_name = 'conditions'
        if self._db.has_table(table_name):
            lines = self._db.select(table_name, fields=['pdf_id'])
            if len(lines) > 0:
                self._pdf_id = lines[0]['pdf_id']

    def _set_pdfparam_conditions(self, pdfparam):
        assert(isinstance(pdfparam, PdfParam))
        table_name = 'conditions'

        # self._pdf_id = pdfparam.digest()
        comment = pdfparam.comment
        num_of_atoms = pdfparam.num_of_atoms
        num_of_AOs = pdfparam.num_of_AOs
        num_of_MOs = pdfparam.num_of_MOs
        method = pdfparam.method
        guess = pdfparam.guess
        xc_functional = pdfparam.xc_functional
        iterations = pdfparam.iterations
        max_iterations = pdfparam.max_iterations
        scf_converged = pdfparam.scf_converged

        where_str = 'pdf_id = "%s"' % (self._pdf_id)
        check_record = self._db.select(table_name, where=where_str)
        if check_record:
            self._db.update(table_name,
                            contents={'comment': comment,
                                      'num_of_atoms': num_of_atoms,
                                      'num_of_AOs': num_of_AOs,
                                      'num_of_MOs': num_of_MOs,
                                      'method': method,
                                      'guess': guess,
                                      'xc_functional': xc_functional,
                                      'iterations': iterations,
                                      'max_iterations': max_iterations,
                                      'scf_converged': scf_converged},
                            where=where_str)
        else:
            self._db.insert(table_name,
                            {'pdf_id': self._pdf_id,
                             'comment': comment,
                             'num_of_atoms': num_of_atoms,
                             'num_of_AOs': num_of_AOs,
                             'num_of_MOs': num_of_MOs,
                             'method': method,
                             'guess': guess,
                             'xc_functional': xc_functional,
                             'iterations': iterations,
                             'max_iterations': max_iterations,
                             'scf_converged': scf_converged})

    def _set_pdfparam_coordinates(self, pdfparam):
        assert(isinstance(pdfparam, PdfParam))
        table_name = 'coordinates'
        if not self._db.has_table(table_name):
            self._db.create_table(table_name,
                                  ['atom_id',
                                   'symbol',
                                   'x', 'y', 'z',
                                   'charge',
                                   'label'],
                                  ['atom_id'])
        # pdfparam.moleculeはbridge.AtomGroupを返す
        for atom_id, atom in pdfparam.molecule.atoms():
            symbol = atom.symbol
            xyz = atom.xyz
            charge = atom.charge
            label = atom.name

            where_str = 'atom_id = "%s"' % (atom_id)
            check_record = self._db.select(table_name, where=where_str)
            if check_record:
                self._db.update(table_name,
                                contents={'symbol': symbol,
                                          'x': xyz.x,
                                          'y': xyz.y,
                                          'z': xyz.z,
                                          'charge': charge,
                                          'label': label},
                                where=where_str)
            else:
                self._db.insert(table_name,
                                {'atom_id': atom_id,
                                 'symbol': symbol,
                                 'x': xyz.x,
                                 'y': xyz.y,
                                 'z': xyz.z,
                                 'charge': charge,
                                 'label': label})

    def _set_pdfparam_basisset(self, pdfparam):
        """
        BasisSet情報を格納する
        """
        assert(isinstance(pdfparam, PdfParam))
        table_name = 'basis'
        if not self._db.has_table(table_name):
            self._db.create_table(table_name,
                                  ['atom_label',
                                   'name'],
                                  ['atom_label'])
        table_name = 'basisset_CGTO'
        if not self._db.has_table(table_name):
            self._db.create_table(table_name,
                                  ['name',
                                   'CGTO_ID',
                                   'shell_type',
                                   'scale_factor'],
                                  ['name', 'CGTO_ID'])
        table_name = 'basisset_PGTO'
        if not self._db.has_table(table_name):
            self._db.create_table(table_name,
                                  ['name',
                                   'CGTO_ID',
                                   'PGTO_ID',
                                   'coef',
                                   'exp'],
                                  ['name', 'CGTO_ID', 'PGTO_ID'])

        for atom_label in pdfparam.get_basisset_atomlabels():
            basisset = pdfparam.get_basisset(atom_label)
            name = basisset.name

            has_atom_label = self._db.select('basis',
                                             where='atom_label="{atom_label}"'.format(
                                                 atom_label=atom_label))
            if has_atom_label:
                self._db.update('basis',
                                contents={'name': name},
                                where='atom_label="{}"'.format(atom_label))
            else:
                self._db.insert('basis',
                                {'atom_label': atom_label,
                                 'name': name})

            for cgto_id, cgto in enumerate(basisset):
                shell_type = cgto.shell_type
                scale_factor = cgto.scale_factor

                has_CGTO = self._db.select('basisset_CGTO',
                                           where='name="{name}" and CGTO_ID={CGTO_ID}'.format(
                                               name=name,
                                               CGTO_ID=cgto_id))
                if has_CGTO:
                    self._db.update('basisset_CGTO',
                                    {'shell_type': shell_type,
                                     'scale_factor': scale_factor
                                     },
                                    where='name="{name}" and CGTO_ID={CGTO_ID}'.format(
                                        name=name,
                                        CGTO_ID=cgto_id))
                else:
                    self._db.insert('basisset_CGTO',
                                    {'name': name,
                                     'CGTO_ID': cgto_id,
                                     'shell_type': shell_type,
                                     'scale_factor': scale_factor
                                     })

                for pgto_id, pgto in enumerate(cgto):
                    coef = pgto.coef
                    exp = pgto.exp

                    has_PGTO = self._db.select('basisset_PGTO',
                                               where='name="{name}" and CGTO_ID={CGTO_ID} and PGTO_ID={PGTO_ID}'.format(
                                                   name=name,
                                                   CGTO_ID=cgto_id,
                                                   PGTO_ID=pgto_id))

                    if has_PGTO:
                        self._db.update('basisset_PGTO',
                                        {'coef': coef,
                                         'exp': exp
                                         },
                                        where='name="{name}" and CGTO_ID="{CGTO_ID}" and PGTO_ID="{PGTO_ID}"'.format(
                                            name=name,
                                            CGTO_ID=cgto_id,
                                            PGTO_ID=pgto_id))
                    else:
                        self._db.insert('basisset_PGTO',
                                        {'name': name,
                                         'CGTO_ID': cgto_id,
                                         'PGTO_ID': pgto_id,
                                         'coef': coef,
                                         'exp': exp
                                         })

    def _set_pdfparam_total_energies(self, pdfparam):
        assert(isinstance(pdfparam, PdfParam))
        table_name = 'total_energies'
        if not self._db.has_table(table_name):
            self._db.create_table(table_name,
                                  ['iteration',
                                   'energy'],
                                  ['iteration'])

        TEs = pdfparam.TEs
        if TEs is not None:
            for iteration, energy in TEs.items():
                iteration = int(iteration)
                where_str = 'iteration = {}'.format(iteration)
                check_record = self._db.select(table_name, where=where_str)
                if check_record:
                    self._db.update(table_name,
                                    contents={'energy': energy},
                                    where=where_str)
                else:
                    self._db.insert(table_name,
                                    {'iteration': iteration,
                                     'energy': energy})

    def _set_pdfparam_gradient(self, pdfparam):
        assert(isinstance(pdfparam, PdfParam))
        table_name = 'gradient'
        if not self._db.has_table(table_name):
            self._create_table_gradient()

        num_of_atoms = pdfparam.num_of_atoms
        for atom_index in range(num_of_atoms):
            v = pdfparam.get_gradient(atom_index)
            if v != None:
                data = {'atom_id': atom_index,
                        'x': v[0],
                        'y': v[1],
                        'z': v[2]}
                where_str = 'atom_id = %d' % (atom_index)
                check_record = self._db.select(table_name, where=where_str)
                if check_record:
                    self._db.update(table_name,
                                    contents=data,
                                    where=where_str)
                else:
                    self._db.insert(table_name, data)
        #

    # ------------------------------------------------------------------
    # pdf_id
    def _get_pdf_id(self):
        return self._pdf_id

    pdf_id = property(_get_pdf_id)

    # ------------------------------------------------------------------
    # num_of_atoms
    def _get_num_of_atoms(self):
        value = None
        table_name = 'conditions'
        if self._db.has_table(table_name):
            values = self._db.select(table_name,
                                     fields=['num_of_atoms'])
            if len(values) > 0:
                value = int(values[0]['num_of_atoms'])
        return value

    num_of_atoms = property(_get_num_of_atoms)

    # ------------------------------------------------------------------
    # num_of_AOs
    def _get_num_of_AOs(self):
        value = None
        table_name = 'conditions'
        if self._db.has_table(table_name):
            values = self._db.select(table_name,
                                     fields=['num_of_AOs'])
            if len(values) > 0:
                value = int(values[0]['num_of_AOs'])
        return value

    num_of_AOs = property(_get_num_of_AOs)

    # ------------------------------------------------------------------
    # num_of_MOs
    def _get_num_of_MOs(self):
        value = None
        table_name = 'conditions'
        if self._db.has_table(table_name):
            values = self._db.select(table_name,
                                     fields=['num_of_MOs'])
            if len(values) > 0:
                value = int(values[0]['num_of_MOs'])
        return value

    num_of_MOs = property(_get_num_of_MOs)

    # ------------------------------------------------------------------
    # method
    def _get_method(self):
        value = None
        table_name = 'conditions'
        if self._db.has_table(table_name):
            values = self._db.select(table_name,
                                     fields=['method'])
            if len(values) > 0:
                value = values[0]['method']
        return value

    method = property(_get_method)

    # ------------------------------------------------------------------
    # iterations
    def _get_iterations(self):
        value = None
        table_name = 'conditions'
        if self._db.has_table(table_name):
            values = self._db.select(table_name,
                                     fields=['iterations'])
            if len(values) > 0:
                value = values[0]['iterations']
        return value

    iterations = property(_get_iterations)

    # ------------------------------------------------------------------
    # scf_converged
    def _get_scf_converged(self):
        value = None
        table_name = 'conditions'
        if self._db.has_table(table_name):
            values = self._db.select(table_name,
                                     fields=['scf_converged'])
            if len(values) > 0:
                value = values[0]['scf_converged']
        return value

    scf_converged = property(_get_scf_converged)

    # 出力 ----------
    def get_molecule(self):
        """
        AtomGroupを返す
        """
        mol = bridge.AtomGroup()
        table_name = 'coordinates'
        if self._db.has_table(table_name):
            atoms = self._db.select(table_name,
                                    fields=['atom_id',
                                            'symbol',
                                            'x',
                                            'y',
                                            'z',
                                            'charge',
                                            'label'])
            num_of_atoms = len(atoms)
            for atom in atoms:
                atom_id = atom['atom_id']
                symbol = atom['symbol']
                x = atom['x']
                y = atom['y']
                z = atom['z']
                charge = atom['charge']
                label = atom['label']
                a = bridge.Atom(symbol=symbol,
                                position=bridge.Position([x, y, z]),
                                charge=charge,
                                name=label)
                mol.set_atom(atom_id, a)
        return mol

    def get_basisset_name(self, atom_label):
        """
        原子名(atom_label)に対応する名前を返す
        """
        atom_label = str(atom_label)
        name = ""
        table_name = 'basis'
        if self._db.has_table(table_name):
            results = self._db.select(table_name,
                                      fields=['name'],
                                      where={'atom_label': atom_label})
            if len(results) > 0:
                name = results[0]['name']
        return name

    def get_basisset(self, name):
        """
        与えられた名前の基底関数情報(Basisオブジェクト)を返す
        """
        name = str(name)
        basisset = BasisSet()
        cgto_table_name = 'basisset_CGTO'
        pgto_table_name = 'basisset_PGTO'
        if self._db.has_table(cgto_table_name):
            cgto_entries = self._db.select(cgto_table_name,
                                           fields=['CGTO_ID',
                                                   'shell_type',
                                                   'scale_factor'],
                                           where={'name': name})
            num_of_CGTOs = len(cgto_entries)
            #print("num_of_CGTOs = %d" % (num_of_CGTOs))
            basisset = BasisSet(name, num_of_CGTOs)
            for cgto_entry in cgto_entries:
                CGTO_ID = cgto_entry['CGTO_ID']
                assert(CGTO_ID < num_of_CGTOs)
                shell_type = cgto_entry['shell_type']
                scale_factor = cgto_entry['scale_factor']
                pgto_entries = self._db.select(pgto_table_name,
                                               fields=['PGTO_ID',
                                                       'coef',
                                                       'exp'],
                                               where={'name': name,
                                                      'CGTO_ID': CGTO_ID})
                num_of_PGTOs = len(pgto_entries)
                basisset[CGTO_ID] = ContractedGTO(shell_type, num_of_PGTOs)
                for pgto_entry in pgto_entries:
                    PGTO_ID = int(pgto_entry['PGTO_ID'])
                    assert(PGTO_ID < num_of_PGTOs)
                    coef = pgto_entry['coef']
                    exp = pgto_entry['exp']
                    basisset[CGTO_ID][PGTO_ID] = PrimitiveGTO(exp, coef)
        return basisset

    def get_HOMO_level(self, runtype):
        """
        HOMOのレベルを返す
        0から始まることに注意。
        """
        answer = None
        v = self.get_occupations(runtype)
        for i in range(len(v) - 1, 0, -1):
            if v[i] > 0.0:
                answer = i
                break
        return answer

    def get_total_energy(self, iteration):
        iteration = int(iteration)
        value = None
        table_name = 'total_energies'
        if self._db.has_table(table_name):
            values = self._db.select(table_name,
                                     fields=['energy'],
                                     where={'iteration': iteration})
            if len(values) > 0:
                value = values[0]['energy']
        return value

    # 入力 ---------------

    # vectors ------------------------------------------------------------------
    def _set_vector(self, dtype, runtype, iteration, v):
        """
        ベクトルを保存する
        """
        dtype = str(dtype)
        dtype = dtype.upper()
        runtype = str(runtype)
        runtype = runtype.upper()
        iteration = int(iteration)

        table_name = 'vectors'
        where_str = 'dtype = "%s" and runtype = "%s" and iteration = %d' % (
            dtype, runtype, iteration)
        check_record = self._db.select(table_name, where=where_str)
        if check_record:
            self._db.update(table_name,
                            contents={'size': len(v),
                                      'data': sqlite3.Binary(v.get_buffer())},
                            where=where_str)
        else:
            self._db.insert(table_name,
                            {'dtype': dtype,
                             'runtype': runtype,
                             'iteration': iteration,
                             'size': len(v),
                             'data': sqlite3.Binary(v.get_buffer())})

    def _get_vector(self, dtype, runtype, iteration):
        """
        ベクトルを取得する
        """
        dtype = str(dtype)
        dtype = dtype.upper()
        runtype = str(runtype)
        runtype = runtype.upper()
        iteration = int(iteration)

        answer = None
        table_name = 'vectors'
        where_str = 'dtype = "%s" and runtype = "%s" and iteration = %d' % (
            dtype, runtype, iteration)
        check_record = self._db.select(table_name,
                                       ['size', 'data'],
                                       where=where_str)
        if check_record:
            size = check_record[0]['size']
            buf = check_record[0]['data']
            answer = bridge.Vector()
            answer.set_buffer(buf)
            assert(size == len(answer))

        return answer

    # matrix -------------------------------------------------------------------
    def _set_matrix(self, dtype, runtype, iteration, m):
        """
        行列を設定する
        """
        dtype = str(dtype)
        dtype = dtype.upper()
        runtype = str(runtype)
        runtype = runtype.upper()
        iteration = int(iteration)
        assert(isinstance(m, bridge.Matrix))

        if m.rows <= 0:
            self._logger.warn(
                'The row size of the input matrix is too small.: {}'.format(m.rows))
            self._logger.warn('cannnot register the matrix.')
            return
        if m.cols <= 0:
            self._logger.warn(
                'The col size of the input matrix is too small.: {}'.format(m.cols))
            self._logger.warn('cannnot register the matrix.')
            return

        table_name = 'matrices'
        where_str = 'dtype = "%s" and runtype = "%s" and iteration = %d' % (
            dtype, runtype, iteration)
        check_record = self._db.select(table_name, where=where_str)
        if check_record:
            self._db.update(table_name,
                            contents={'rows': m.rows,
                                      'cols': m.cols,
                                      'format': m.type,
                                      'data': sqlite3.Binary(m.get_buffer())
                                      },
                            where=where_str)
        else:
            self._db.insert(table_name,
                            {'dtype': dtype,
                             'runtype': runtype,
                             'iteration': iteration,
                             'rows': m.rows,
                             'cols': m.cols,
                             'format': m.type,
                             'data': sqlite3.Binary(m.get_buffer())})

    def _get_matrix(self, dtype, runtype, iteration):
        """
        行列を取得する
        """
        dtype = str(dtype)
        dtype = dtype.upper()
        runtype = str(runtype)
        runtype = runtype.upper()
        iteration = int(iteration)

        answer = None
        table_name = 'matrices'
        where_str = 'dtype = "%s" and runtype = "%s" and iteration = %d' % (
            dtype, runtype, iteration)
        check_record = self._db.select(table_name,
                                       ['rows',
                                        'cols',
                                        'format',
                                        'data'],
                                       where=where_str)
        if check_record:
            rows = check_record[0]['rows']
            cols = check_record[0]['cols']
            mat_format = check_record[0]['format']
            buf = check_record[0]['data']
            if mat_format == 'SY':
                assert(rows == cols)
                answer = bridge.SymmetricMatrix(rows)
            else:
                assert(mat_format == 'GE')
                answer = bridge.Matrix(rows, cols)
            answer.set_buffer(buf)

        return answer

    # occ ----------------------------------------------------------------------
    def set_occupations(self, runtype, occ_level):
        """
        占有軌道情報を設定する
        """
        self._set_vector('occ', runtype, 0, occ_level)

    def get_occupations(self, runtype):
        """
        占有軌道情報を取得する
        """
        return self._get_vector('occ', runtype, 0)

    # energy level -------------------------------------------------------------
    def set_energylevel(self, runtype, iteration, energy_levels):
        """
        エネルギー準位を設定する
        """
        self._set_vector('energy_level', runtype, iteration, energy_levels)

    def get_energylevel(self, runtype, iteration):
        """
        エネルギー準位を取得する
        """
        return self._get_vector('energy_level', runtype, iteration)

    # LCAO ---------------------------------------------------------------------
    def set_lcao(self, runtype, iteration, lcao):
        """
        LCAOを設定する
        """
        self._set_matrix('C', runtype, iteration, lcao)

    # population ------------------------------------------------------------------
    def set_population(self, pop_type, iteration, pop_mat):
        pop_type = pop_type.upper()
        assert(pop_type in ['MULLIKEN'])
        self._set_matrix(dtype=pop_type,
                         runtype='',
                         iteration=iteration,
                         m=pop_mat)

    def get_population(self, pop_type, iteration):
        pop_type = pop_type.upper()
        assert(pop_type in ['MULLIKEN'])
        mat = self._get_matrix(dtype=pop_type,
                               runtype='',
                               iteration=iteration)

        return mat

    # gradient --------------------------------------------------------------------
    def get_gradient(self, atom_id):
        atom_id = int(atom_id)
        gradient = None
        table_name = 'gradient'
        if self._db.has_table(table_name):
            results = self._db.select(
                table_name, fields=['x', 'y', 'z'], where={'atom_id': atom_id})
            if len(results) > 0:
                gradient = [results[0]['x'], results[0]['y'], results[0]['z']]
        return gradient

    def get_gradient_rms(self):
        """
        gradientのRMSを返す
        """
        rms = 0.0
        for atom_index in range(self.num_of_atoms):
            v = self.get_gradient(atom_index)
            if v != None:
                rms += v[0]*v[0] + v[1]*v[1] + v[2]*v[2]
            else:
                return None
        return math.sqrt(rms / (3 * self.num_of_atoms))

    # compare ------------------------------------------------------------------

    def __eq__(self, other):
        answer = True
        answer = answer & self.compare_info(other)
        answer = answer & self.compare_energy(other)

        # check Total Energy
        iterations = max(self.iterations, other.iterations)
        for itr in range(1, iterations + 1):
            TE1 = self.get_total_energy(itr)
            TE2 = other.get_total_energy(itr)
            answer = answer & self._check(TE1, TE2, '#%d TE' % (itr))

        return answer

    def _check(self, val1, val2, msg, allowable_error=pdf.error):
        answer = False
        if isinstance(val1, float) and isinstance(val2, float):
            v = math.fabs(val1 - val2)
            answer = (v < allowable_error)
        elif isinstance(val1, list) and isinstance(val2, list):
            if len(val1) == len(val2):
                for i in range(len(val1)):
                    answer = answer ^ self._check(val1[i], val2[i], '{0}[{1}]'.format(msg, i), allowable_error)
        else:
            answer = (val1 == val2)

        if answer:
            self._logger.debug('test: {msg}; {val1} == {val2} (threshold={threshold})'.format(
                msg=str(msg), val1=val1, val2=val2, threshold=allowable_error))
        else:
            self._logger.error('test: {msg}; {val1} != {val2} (threshold={threshold})'.format(
                msg=str(msg), val1=val1, val2=val2, threshold=allowable_error))

        return answer

    def compare_info(self, rhs):
        answer = True

        answer = answer & self._check(self.num_of_AOs, rhs.num_of_AOs,
                                      'num_of_AOs')
        answer = answer & self._check(self.num_of_MOs, rhs.num_of_MOs,
                                      'num_of_MOs')
        answer = answer & self._check(self.scf_converged, rhs.scf_converged,
                                      'scf_converged')

        return answer

    def compare_energy(self, rhs, threshold=pdf.error):
        answer = True

        TE_self = self.get_total_energy(self.iterations)
        TE_other = rhs.get_total_energy(rhs.iterations)
        answer = answer & self._check(TE_self, TE_other,
                                      'TE', threshold)

        return answer

    def compare_pop(self, rhs):
        answer = True

        answer = answer & self._check(self.num_of_atoms,
                                      rhs.num_of_atoms,
                                      'num of atoms')
        num_of_atoms = self.num_of_atoms

        itr1 = self.iterations
        itr2 = rhs.iterations
        pop1 = self.get_population('mulliken', itr1)
        pop2 = rhs.get_population('mulliken', itr2)

        if isinstance(pop1, bridge.Matrix):
            if isinstance(pop2, bridge.Matrix):
                answer = answer & self._check(pop1.rows,
                                              pop2.rows,
                                              'pop rows')
                answer = answer & self._check(pop1.cols,
                                              pop2.cols,
                                              'pop cols')
                answer = answer & self._check(pop1, pop2,
                                              'pop matrix')
            else:
                self._logger.warning('type mismatch(rhs): POP not found')
                answer = False
        else:
            self._logger.warning('type mismatch(lhs): POP not found')
            answer = False

        return answer

    def compare_grad(self, rhs):
        answer = True

        # check gradient elements
        max1 = 0.0
        max2 = 0.0
        for atom_index in range(self.num_of_atoms):
            xyz1 = self.get_gradient(atom_index)
            xyz2 = rhs.get_gradient(atom_index)
            answer = answer & self._check(xyz1, xyz2, 'grad')
            max1 = max(max1, max(max(xyz1), -min(xyz1)))
            max2 = max(max2, max(max(xyz2), -min(xyz2)))
        answer = answer & self._check(max1, max2, 'grad max')

        # check gradient RMS
        answer = answer & self._check(self.get_gradient_rms(),
                                      rhs.get_gradient_rms(),
                                      'grad RMS')
        return answer


if __name__ == '__main__':
    pass

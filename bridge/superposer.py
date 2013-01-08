#!/usr/bin/env python
# -*- coding: utf-8 -*-

import math
import copy

from position import Position

class Superposer(object):
    """
    """
    
    def __init__(self, atom_group1, atom_group2):
        self._translation_vct2 = None
        self._rotation_mat = None
        self._rmsd = None

        self._calc(copy.deepcopy(atom_group1),
                   copy.deepcopy(atom_group2))
        
    @property
    def rotation_mat(self):
        return self._rotation_mat

    @property
    def rmsd(self):
        return self._rmsd

    def _calc(self, atom_group1, atom_group2):
        (positions1, positions2) = self._match_positions(atom_group1, atom_group2)
        num_of_positions = len(positions1)
        assert(num_of_positions == len(positions2))
        
        (translation_vct1, self._translation_vct2) = self._fix_positions(positions1, positions2)
        self._rotation_mat = self._get_rotation_matrix(num_of_positions,
                                                       positions1, positions2)
        
        self._update_positions(positions1, positions2,
                               self._rotation_mat,
                               self._translation_vct2)

        self._rmsd = self._calc_rmsd(positions1, positions2)
        
    def _match_positions(self, atom_group1, atom_group2):
        """
        compare two atom_group objects, and select common points.
        return the list of list corresponding two points.
        """
        positions1 = []
        positions2 = []
        for key, ag1 in atom_group1.groups():
            if atom_group2.has_group(key):
                (p1, p2) = self._match_positions(ag1,
                                                 atom_group2.get_group(key))
                positions1 += p1
                positions2 += p2

        for key, atom1 in atom_group1.atoms():
            if atom_group2.has_atom(key):
                positions1 += [atom1.position]
                positions2 += [atom_group2.get_atom(key).position]

        return (positions1, positions2)

    def _fix_positions(self, positions1, positions2):
        center1 = self._calc_center(positions1)
        center2 = self._calc_center(positions2)

        translation_vct1 = - center1
        translation_vct2 =   center2

        for i in range(len(positions1)):
            positions1[i] -= center1
        for i in range(len(positions2)):
            positions2[i] -= center2

        return (translation_vct1, translation_vct2)


    def _calc_center(self, positions):
        """
        return the center position of input positions
        """
        c = Position()

        num_of_positions = len(positions)
        for i in range(num_of_positions):
            c += positions[i]
        c /= num_of_positions
        
        return c
            
    def _get_rotation_matrix(self, num_of_points,
                             positions1, positions2):
        r = Matrix(3, 3)

        # r_ij = Sum_over_k{p2(k, i) * p1(k, j)}
        for k in range(num_of_points):
            x1 = positions1[k].x
            y1 = positions1[k].y
            z1 = positions1[k].z
            x2 = positions2[k].x
            y2 = positions2[k].y
            z2 = positions2[k].z

            r.add(0, 0, x2 * x1)
            r.add(0, 1, x2 * y1)
            r.add(0, 2, x2 * z1)
            r.add(1, 0, y2 * x1)
            r.add(1, 1, y2 * y1)
            r.add(1, 2, y2 * z1)
            r.add(2, 0, z2 * x1)
            r.add(2, 1, z2 * y1)
            r.add(2, 2, z2 * z1)

        tr = r.copy()
        tr.transpose()
        trr = tr * r
        trr = trr.get_symmetric_matrix()
        
        eigval, eigvec = trr.eig()
        a = self._make_right_handed(eigvec)
        
        b = Matrix(3, 3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    v = r.get(j , k) * a.get(i, k)
                    b.add(i, j, v) 

            # normalize b[i]
            w = 0.0
            for j in range(3):
                w += b.get(i, j) * b.get(i, j)

            t = math.sqrt(1.0 / w)
            for j in range(3):
                v = b.get(i, j)
                b.set(i, j, v * t)
                
        # b[2] = b[0] x b[1]
        tmp_vct = self._calc_vector_product(b.get_row_vector(0),
                                            b.get_row_vector(1))
        for i in range(3):
            b.set(2, i, tmp_vct[i])

        # rotation matrix r_ij = b_ki * a_kj
        mat = self._set_rotation(a, b)
        return mat
        
    def _make_right_handed(self, mat):
        assert(mat.get_num_of_rows() == 3)
        assert(mat.get_num_of_cols() == 3)

        v1 = Vector(3)
        v2 = Vector(3)
        for i in range(3):
            v1[i] = mat.get(0, i)
            v2[i] = mat.get(1, i)

        v3 = self._calc_vector_product(v1, v2)
        
        answer = Matrix(3, 3)
        for i in range(3):
            answer.set(0, i, v1[i])
            answer.set(1, i, v2[i])
            answer.set(2, i, v3[i])

        return answer
            
    def _calc_vector_product(self, v1, v2):
        assert(len(v1) == 3)
        assert(len(v2) == 3)

        v3 = Vector(3)
        v3[0] = v1[1] * v2[2] - v1[2] * v2[1]
        v3[1] = v1[2] * v2[0] - v1[0] * v2[2]
        v3[2] = v1[0] * v2[1] - v1[1] * v2[0]

        return v3
    
    def _set_rotation(self, a, b):
        assert(a.get_num_of_rows() == 3)
        assert(a.get_num_of_cols() == 3)
        assert(b.get_num_of_rows() == 3)
        assert(b.get_num_of_cols() == 3)
    
        r = Matrix(3, 3)
        for i in range(3):
            for j in range(3):
                for k in range(3):
                    v = b.get(k, i) * a.get(k, j)
                    r.add(i, j, v)
        return r

    def _update_positions(self, positions1, positions2,
                          rotation_mat, translation_vct2):
        for i in range(len(positions1)):
            positions1[i].rotate(rotation_mat)
            
        for i in range(len(positions1)):
            positions1[i] += translation_vct2
        for i in range(len(positions2)):
            positions2[i] += translation_vct2

    def _calc_rmsd(self, positions1, positions2):
        """
        calc rmsd.
        store the value to self._rmsd
        """
        num_of_positions = min(len(positions1),
                               len(positions2))
        msd = 0.0
        for i in range(num_of_positions):
            msd += positions1[i].square_distance_from(positions2[i])
        msd /= float(num_of_positions)
        rmsd = math.sqrt(msd)

        return rmsd
    

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

import os
import sys
import argparse
try:
    import msgpack
except:
    import msgpack_pure as msgpack
    
import pdfpytools as pdf

def main():
    A = get_matrix_A()
    B = get_matrix_B()

    A.save("A.mat")
    B.save("B.mat")
    
def get_matrix_A():
    num_of_rows = 3
    num_of_cols = 3
    A = pdf.Matrix(num_of_rows, num_of_cols)
    count = 0.0
    for r in range(num_of_rows):
        for c in range(num_of_cols):
            A.set(r, c, count)
            count += 1.0
    return A

def get_matrix_B():
    num_of_rows = 3
    num_of_cols = 3
    B = pdf.Matrix(num_of_rows, num_of_cols)
    count = 0.0
    for c in range(num_of_cols):
        for r in range(num_of_rows):
            B.set(r, c, count)
            count += 1.0
    return B

if __name__ == '__main__':
    main()

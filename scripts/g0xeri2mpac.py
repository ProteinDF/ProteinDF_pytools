#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2002-2014 The ProteinDF project
# see also AUTHORS and README.
# 
# This file is part of ProteinDF.
# 
# ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
# 
# ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

"""
Gaussianのログファイルから2電子積分部分を抽出し、
message pack形式で出力する。

Gaussian入力ファイルの例:
#p hf/gen sp 5D 
Pop=Full GFPrint GFInput 
Integral(NoJEngine)
ExtraLinks=L316
NoRaff NoSymm
SCF=Conventional
IOP(3/33=3)

"""

import sys
import re
import argparse
try:
    import msgpack
except:
    import msgpack_pure as msgpack



def parse_gau_output(path):
    """
    read Gaussian output file
    """
    re_begin_dump_2e = re.compile('^.*Dumping Two-Electron integrals.*$')
    re_end_dump_2e = re.compile('^.*Leave Link  316.*$')
    re_dump_2e_int = re.compile('^.*I=\s*(\d+)\s+J=\s*(\d+)\s+K=\s*(\d+)\s+L=\s*(\d+)\s+Int=\s*(\S+)')

    integrals = []
    
    is_2e_block = False
    fin = open(path, 'r')
    for line in fin:
        line = line.rstrip()
        
        MatchObj = re_begin_dump_2e.match(line)
        if MatchObj:
            is_2e_block = True

        MatchObj = re_end_dump_2e.match(line)
        if MatchObj:
            is_2e_block = False

        if is_2e_block:
            MatchObj = re_dump_2e_int.match(line)
            if MatchObj:
                I = int(MatchObj.group(1)) -1
                J = int(MatchObj.group(2)) -1 
                K = int(MatchObj.group(3)) -1
                L = int(MatchObj.group(4)) -1
                Int = float(MatchObj.group(5).replace('D', 'E'))

                #sys.stdout.write('I=%d J=%d K=%d L=%d: % e\n' % (I, J, K, L, Int))
                integral = [I, J, K, L, Int]
                integrals.append(integral)
        
    return integrals


def main():
    parser = argparse.ArgumentParser(description='Gaussian output to mspack')
    parser.add_argument('GAU_FILE',
                        nargs=1,
                        help='Gaussian output file')
    parser.add_argument('MPAC_FILE',
                        nargs=1,
                        help='msgpack file')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    verbose = args.verbose
    gau_file = args.GAU_FILE[0]
    mpac_file = args.MPAC_FILE[0]

    integrals = parse_gau_output(gau_file)

    mpac_data = msgpack.packb(integrals)
    
    fout = open(mpac_file, 'wb')
    fout.write(mpac_data)
    fout.close()


if __name__ == '__main__':
    main()


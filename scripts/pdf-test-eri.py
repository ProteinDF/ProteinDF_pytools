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
ERIの正しい解とProteinDFのERI値を比較する

ERIの正解はGaussianなどによって求める
"""

import os
import sys
import re
import math
import subprocess
import argparse

import proteindf_bridge as bridge
import proteindf_tools as pdf

# parameters
PDF_INDEX_FILE = 'pdfindex.mpac'


def make_orbinfo_table(pdfparam_path):
    """
    """
    # get orbital information
    orb_info_path = 'orbinfo.mpac'
    PDF_HOME = os.environ['PDF_HOME']
    cmd = '%s/bin/pdf-info-orb' % (PDF_HOME)
    args = [cmd,
            '-p', pdfparam_path,
            '-w', orb_info_path
            ]
    subproc_args = {'stdin': subprocess.PIPE,
                    'stdout': subprocess.PIPE,
                    'stderr': subprocess.STDOUT,
                    'close_fds': True,
                    }
    try:
        proc = subprocess.Popen(args, **subproc_args)
    except OSError:
        sys.exit('Failed to execute command: %s\n' % args[0])

    (stdouterr, stdin) = (proc.stdout, proc.stdin)
    while True:
        line = stdouterr.readline()
        if not line:
            break
        print line.rstrip()
    ret = proc.wait()
    #sys.stdout.write('%s retrun code: %d\n' % (args[0], ret))

    orbinfo = load_mpac(orb_info_path)
    return orbinfo


def translate_5d_index(orbinfo, p, q, r, s):
    """
    Gaussianの5dが(z2, xz, yz, x2-y2, xy)の順なのに対し、
    ProteinDFでは(xy, xz, yz, x2-y2, z2)の順なので、
    入力されたインデックスがGaussian 5d軌道の場合、
    対応するProteinDF 5d軌道のインデックスを返す。
    """
    p = gau2pdf_corresponding_index(p, orbinfo[p]['basis_type_id'])
    q = gau2pdf_corresponding_index(q, orbinfo[q]['basis_type_id'])
    r = gau2pdf_corresponding_index(r, orbinfo[r]['basis_type_id'])
    s = gau2pdf_corresponding_index(s, orbinfo[s]['basis_type_id'])

    return (p, q, r, s)


def gau2pdf_corresponding_index(orb, basis_type_id):
    new_orb = orb
    if basis_type_id == 4:
        new_orb = orb + 4
    elif basis_type_id == 5:
        pass
    elif basis_type_id == 6:
        pass
    elif basis_type_id == 7:
        pass
    elif basis_type_id == 8:
        new_orb = orb - 4

    return new_orb


def get_pdf_eri_list(pdfparam_path, pdf_eri_indeces_path):
    re_eri_val = re.compile('^.*=\s*(\S+)')

    PDF_HOME = os.environ['PDF_HOME']
    cmd = '%s/bin/pdf-eri' % (PDF_HOME)
    args = [cmd,
            '-p', pdfparam_path,
            str(pdf_eri_indeces_path)]
    subproc_args = {'stdin': subprocess.PIPE,
                    'stdout': subprocess.PIPE,
                    'stderr': subprocess.STDOUT,
                    'close_fds': True,
                    }
    try:
        proc = subprocess.Popen(args, **subproc_args)
    except OSError:
        sys.exit('Failed to execute command: %s' % args[0])

    eri_values = []
    (stdouterr, stdin) = (proc.stdout, proc.stdin)
    while True:
        line = stdouterr.readline()

        matchObj = re_eri_val.match(line)
        if matchObj:
            value = float(matchObj.group(1))
            eri_values.append(value)
        if not line:
            break
    ret = proc.wait()
    return eri_values


def check_eri(value1, value2, threshold, p, q, r, s,
              verbose=False):
    answer = True
    err = math.fabs(value1 - value2)
    if err > threshold:
        sys.stderr.write(
            'NG (%2d %2d|%2d %2d) % 8.5f != % 8.5f (% 8.5f)\n' % (
                p, q, r, s, value1, value2, err))
        answer = False
    else:
        if verbose:
            sys.stderr.write(
                'OK (%2d %2d|%2d %2d) % 8.5f == % 8.5f\n' % (
                    p, q, r, s, value1, value2))

    return answer


def main():
    parser = argparse.ArgumentParser(description='test ERI value')
    parser.add_argument('ERI_MPAC_FILE',
                        nargs=1,
                        help='ERI msgpack file')
    parser.add_argument('PDF_PARAM_FILE',
                        nargs=1,
                        help='ProteinDF parameter file')
    parser.add_argument('-t', '--threshold',
                        nargs=1,
                        action='store',
                        default=1.0E-5)
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    verbose = args.verbose
    eri_mpac_file = args.ERI_MPAC_FILE[0]
    pdfparam_file = args.PDF_PARAM_FILE[0]
    threshold = args.threshold

    # ERI check
    error_count = 0
    total_count = 0
    eri_data = load_mpac(eri_mpac_file)

    num_of_tests = len(eri_data)

    # make index list
    orbinfo = make_orbinfo_table(pdfparam_file)
    pdf_eri_indeces = [[0, 0, 0, 0] for x in range(num_of_tests)]
    ref_eri_values = [0.0 for x in range(num_of_tests)]
    for i in range(num_of_tests):
        # 入力インデックスはGaussianのインデックス
        index4 = translate_5d_index(orbinfo,
                                    eri_data[i][0],
                                    eri_data[i][1],
                                    eri_data[i][2],
                                    eri_data[i][3])
        pdf_eri_indeces[i] = index4
        ref_eri_values[i] = eri_data[i][4]

    # file out
    bridge.save_msgpack(pdf_eri_indeces, PDF_INDEX_FILE)

    # calc eri
    pdf_eri_values = get_pdf_eri_list(pdfparam_file, PDF_INDEX_FILE)
    assert(len(pdf_eri_values) == num_of_tests)

    # check
    for i in range(num_of_tests):
        ref_eri_value = ref_eri_values[i]
        pdf_eri_value = pdf_eri_values[i]

        if not check_eri(pdf_eri_value, ref_eri_value,
                         threshold,
                         pdf_eri_indeces[i][0],
                         pdf_eri_indeces[i][1],
                         pdf_eri_indeces[i][2],
                         pdf_eri_indeces[i][3],
                         verbose):
            error_count += 1
        total_count += 1

    if error_count:
        sys.exit('%d/%d ERROR FOUND' % (error_count, total_count))
    else:
        sys.stderr.write('pass %d tests.\n' % (total_count))


if __name__ == '__main__':
    main()

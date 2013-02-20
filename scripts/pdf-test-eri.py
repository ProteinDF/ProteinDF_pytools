#!/usr/bin/env python
# -*- coding: utf-8 -*-

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
try:
    import msgpack
except:
    import msgpack_pure as msgpack

import bridge
import pdf


def translate_5d_index(pdfparam_path, p, q, r, s):
    """
    Gaussianの5dが(z2, xz, yz, x2-y2, xy)の順なのに対し、
    ProteinDFでは(xy, yz, xz, x2-y2, z2)の順なので、
    入力されたインデックスがGaussian 5d軌道の場合、
    対応するProteinDF 5d軌道のインデックスを返す。
    """

    # get orbital information
    orb_info_path = 'orbinfo.mpac'
    PDF_HOME = os.environ['PDF_HOME']
    cmd = '%s/bin/pdf-info-orb' % (PDF_HOME)
    args = [cmd,
            '-p', pdfparam_path,
            '-w', orb_info_path,
            str(p), str(q), str(r), str(s)]
    subproc_args = { 'stdin': subprocess.PIPE,
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

    orb_info = bridge.mpac2py(orb_info_path);
    #print(orb_info)
    #print("p=", p)

    p = gau2pdf_corresponding_index(p, orb_info[p]['basis_type_id'])
    q = gau2pdf_corresponding_index(q, orb_info[q]['basis_type_id'])
    r = gau2pdf_corresponding_index(r, orb_info[r]['basis_type_id'])
    s = gau2pdf_corresponding_index(s, orb_info[s]['basis_type_id'])

    return (p, q, r, s)


def gau2pdf_corresponding_index(orb, basis_type_id):
    new_orb = orb
    if basis_type_id == 4:
        # dxy to z2
        new_orb = orb + 4
    elif basis_type_id == 5:
        # dyz to xz
        new_orb = orb +1 
    elif basis_type_id == 6:
        # dxz to yz
        new_orb = orb - 1
    elif basis_type_id == 7:
        # dx2-y2 to dx2-y2
        pass
    elif basis_type_id == 8:
        # dz2 to xy
        new_orb = orb - 4

    return new_orb


def get_pdf_eri(pdfparam_path,
                p, q, r, s,
                verbose = False):
    re_eri_val = re.compile('^.*=\s*(\S+)')

    PDF_HOME = os.environ['PDF_HOME']
    cmd = '%s/bin/pdf-eri' % (PDF_HOME)
    args = [cmd,
            '-p', pdfparam_path,
            str(p), str(q), str(r), str(s)]
    subproc_args = { 'stdin': subprocess.PIPE,
                     'stdout': subprocess.PIPE,
                     'stderr': subprocess.STDOUT,
                     'close_fds': True,
                     }
    try:
        proc = subprocess.Popen(args, **subproc_args)
    except OSError:
        sys.exit('Failed to execute command: %s' % args[0])

    pdf_integral = 0.0
    (stdouterr, stdin) = (proc.stdout, proc.stdin)
    while True:
        line = stdouterr.readline()

        matchObj = re_eri_val.match(line)
        if matchObj:
            pdf_integral = float(matchObj.group(1))
        if not line:
            break
    ret = proc.wait()
    #sys.stdout.write('%s retrun code: %d' % (args[0], ret))
    return pdf_integral


def check_eri(value1, value2, threshold, p, q, r, s,
              verbose = False):
    answer = True
    err = math.fabs(value1 - value2)
    if err > threshold:
        sys.stderr.write(
            'NG (%d %d|%d %d) % f != %f (% f)\n' % (
                p, q, r, s, value1, value2, err))
        answer = False
    else:
        if verbose:
            sys.stderr.write(
                'OK (%d %d|%d %d) % f == %f\n' % (
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
    eri_data = bridge.mpac2py(eri_mpac_file)
    for eri in eri_data:
        gau_eri_value = eri[4]
        (p, q, r, s) = translate_5d_index(pdfparam_file,
                                          eri[0], eri[1], eri[2], eri[3])
        pdf_eri_value = get_pdf_eri(pdfparam_file,
                                    p, q, r, s,
                                    verbose)
        if not check_eri(pdf_eri_value, gau_eri_value, threshold, p, q, r, s,
                         verbose):
            error_count += 1
        total_count += 1

    if error_count:
        sys.exit('%d/%d ERROR FOUND' % (error_count, total_count))
    else:
        sys.stderr.write('pass %d tests.\n' % (total_count))


if __name__ == '__main__':
    main()

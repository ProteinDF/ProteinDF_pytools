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
ProteinDF コマンドを実行する

* このコマンドが実行されたディレクトリにある(サブ)コマンドを実行する
"""

import sys
import os.path
import argparse
import copy
import pprint

import proteindf_bridge as bridge
import proteindf_tools as pdf


def main():
    # parse args
    parser = argparse.ArgumentParser(description='ProteinDF command helper')
    parser.add_argument('cmd',
                        nargs='+',
                        help='command')
    parser.add_argument("--DEBUG",
                        action="store_true",
                        default=False)
    parser.add_argument("--OUTPUT",
                        nargs=1,
                        action="store",
                        default="")
    args, unknown = parser.parse_known_args()

    debug = args.DEBUG
    output = ''
    if len(args.OUTPUT) > 0:
        output = args.OUTPUT[0]
    arg_array = args.cmd + unknown

    # PDF_HOMEを設定する
    #pdfcmd_path = os.path.abspath(os.path.dirname(sys.argv[0]))
    #pdfcmd_dirname = os.path.dirname(pdfcmd_path)
    #pdf_home = pdfcmd_dirname
    # if pdfcmd_dirname[-4:] == '/bin':
    #    pdf_home = pdfcmd_dirname[0:-5]
    # if debug:
    #    sys.stderr.write('PDF_HOME environment variable is set to \'{}\'\n'.format(pdf_home))
    #os.environ['PDF_HOME'] = pdf_home

    (yn, pdfcmd, pdfargs, ext) = get_cmd(arg_array)
    if yn == False:
        sys.stderr.write("cannot command: {}".format(
            ' '.join(args.cmd + unknown)))
    subproc_cmd = " ".join([pdfcmd] + pdfargs)

    if debug:
        sys.stderr.write('cmd={}\n'.format(subproc_cmd))

    p = pdf.Process(logfile_path=output)
    return_code = p.cmd(subproc_cmd).commit()

    if return_code != 0:
        sys.stderr.write(
            'Failed to execute pdf command: {}\n'.format(subproc_cmd))
        sys.stderr.write('return code = {}\n'.format(return_code))

    return return_code


def get_cmd(in_args):
    paths = []

    # add the directory where this command is
    current_path = os.path.abspath(os.path.dirname(sys.argv[0]))
    paths.append(current_path)

    # add PDF_HOME/bin
    if 'PDF_HOME' in os.environ:
        PDF_HOME_bin = os.path.join(os.environ["PDF_HOME"], "bin")
        paths.append(PDF_HOME_bin)

    # add PATH
    if 'PATH' in os.environ:
        env_PATHs = os.environ['PATH'].split(os.pathsep)
        paths.extend(env_PATHs)

    # search
    yn = False
    cmd = None
    args = None
    ext = None
    for i in paths:
        args = copy.copy(in_args)
        (yn, cmd, args, ext) = get_cmd_in_dir(args, i)
        if yn == True:
            break

    return (yn, cmd, args, ext)


def get_cmd_in_dir(arg_array, path):
    """
    引数で指定された文字列配列から'pdf-'で始まるコマンドをpathディレクトリから検索する。
    成否, コマンドとその引数(リスト)、拡張子をタプルにして返す
    """
    assert(len(arg_array) > 0)

    subcmd = arg_array[0]
    arg_array.pop(0)
    if subcmd == 'help':
        subcmd = arg_array[0]
        arg_array.pop(0)
        arg_array.insert(0, '-h')
        if subcmd[0] == '-':
            print('illegal command format: %s' % (subcmd))
            return (None, None)

    cmd = '%s/pdf-%s' % (path, subcmd)

    found = False
    ext_list = ['', '.sh', '.py', '.exe']
    for ext in ext_list:
        if os.path.isfile(cmd+ext):
            cmd += ext
            found = True
            break

    return (found, cmd, arg_array, ext)


if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)

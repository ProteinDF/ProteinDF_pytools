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
ProteinDF コマンドを実行する

* このコマンドが実行されたディレクトリにある(サブ)コマンドを実行する
"""

import sys
import os.path
import argparse
import subprocess
import types

def main():
    # parse args
    parser = argparse.ArgumentParser(description='ProteinDF command helper')
    parser.add_argument('cmd',
                        nargs='+',
                        help='command')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    args, unknown = parser.parse_known_args()
    
    verbose = args.verbose
    arg_array = args.cmd + unknown

    # PDF_HOMEを設定する
    pdfcmd_path = os.path.abspath(os.path.dirname(sys.argv[0]))
    pdfcmd_dirname = os.path.dirname(pdfcmd_path)
    pdf_home = pdfcmd_dirname
    if pdfcmd_dirname[-4:] == '/bin':
        pdf_home = pdfcmd_dirname[0:-5]
    if verbose:
        print('PDF_HOME environment variable is set to \'%s\'' % (pdf_home))
    os.environ['PDF_HOME'] = pdf_home
    
    (pdfcmd, pdfargs) = get_cmd(arg_array)
    if verbose:
        print('exec cmd: \'%s\'' % (pdfcmd))
        print('args    : \'%s\'' % (pdfargs))
    
    args = [pdfcmd] + pdfargs
    #subproc_args = {'stdin': None,
    #                'stdout': subprocess.PIPE,
    #                'stderr': subprocess.PIPE,
    #                }
    subproc_args = {'stdin': None,
                    'stderr': subprocess.STDOUT,
                    }
    try:
        proc = subprocess.Popen(args, **subproc_args)
    except OSError:
        print('Failed to execute command: %s' % args[0])
        sys.exit(1)

    return_code = proc.wait()
    (stdouterr, stdin) = (proc.stdout, proc.stdin)

    return return_code

def get_cmd(arg_array):
    """
    引数で指定された文字列配列から'pdf-'で始まるコマンドを検索、
    コマンドとその引数(リスト)をタプルにして返す
    """
    assert(isinstance(arg_array, types.ListType))
    assert(len(arg_array) > 0)

    current_path = os.path.abspath(os.path.dirname(sys.argv[0]))

    subcmd = arg_array[0]
    arg_array.pop(0)
    if subcmd == 'help':
        subcmd = arg_array[0]
        arg_array.pop(0)
        arg_array.insert(0, '-h')
        if subcmd[0] == '-':
            print('illegal command format: %s' % (subcmd))
            return (None, None)
        
    cmd = '%s/pdf-%s' % (current_path, subcmd)
    
    ext_list = ['.sh', '.py', '.exe']
    for ext in ext_list:
        if os.path.isfile(cmd + ext):
            cmd += ext
            break
    
    return (cmd, arg_array)
            
if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)

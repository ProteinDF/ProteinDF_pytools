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
import shlex
import subprocess
import types
import errno
import pprint

def main():
    # parse args
    parser = argparse.ArgumentParser(description='ProteinDF command helper')
    parser.add_argument('cmd',
                        nargs='+',
                        help='command')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    parser.add_argument("-o", "--output",
                        nargs=1,
                        action="store",
                        default="")
    args, unknown = parser.parse_known_args()
    
    verbose = args.verbose
    output = None
    if len(args.output) > 0:
        output = args.output[0]
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
    
    (pdfcmd, pdfargs, ext) = get_cmd(arg_array)
    subproc_cmd = " ".join([pdfcmd] + pdfargs)
    subproc_cmd = shlex.split(subproc_cmd)
    if verbose:
        print('cmd="{}"'.format(subproc_cmd))
    
    use_shell = False
    if ext == '.sh':
        use_shell = True
    subproc_args = {'stdin': None,
                    'stdout': None,
                    'stderr': subprocess.STDOUT,
                    'shell': use_shell
                    }
    output_file = None
    if output != None:
        output_file = open(output, "a")
        subproc_args['stdout'] = output_file
        
    try:
        proc = subprocess.Popen(args = subproc_cmd, **subproc_args)
    except OSError as e:
        print('Failed to execute command: %s' % subproc_cmd)
        print(errno.errorcode[e.errno])
        print(os.strerror(e.errno))
        raise e
        #sys.exit(1)

    return_code = proc.wait()
    (stdouterr, stdin) = (proc.stdout, proc.stdin)
    if verbose:
        print('return code={}'.format(return_code))
    
    if output_file != None:
        output_file.close()

    return return_code

def get_cmd(arg_array):
    """
    引数で指定された文字列配列から'pdf-'で始まるコマンドを検索、
    コマンドとその引数(リスト)、拡張子をタプルにして返す
    """
    #assert(isinstance(arg_array, types.ListType))
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
    
    return (cmd, arg_array, ext)
            
if __name__ == '__main__':
    return_code = main()
    sys.exit(return_code)

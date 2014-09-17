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
    pdfcmd_path = os.path.abspath(os.path.dirname(sys.argv[0]))
    pdfcmd_dirname = os.path.dirname(pdfcmd_path)
    pdf_home = pdfcmd_dirname
    if pdfcmd_dirname[-4:] == '/bin':
        pdf_home = pdfcmd_dirname[0:-5]
    if debug:
        print('PDF_HOME environment variable is set to \'%s\'' % (pdf_home))
    os.environ['PDF_HOME'] = pdf_home

    (pdfcmd, pdfargs, ext) = get_cmd(arg_array)
    subproc_cmd = " ".join([pdfcmd] + pdfargs)
    subproc_cmd = shlex.split(subproc_cmd)
    if debug:
        print('cmd={}'.format(subproc_cmd))
        
    use_shell = False
    if ext == '.sh':
        use_shell = True
    subproc_args = {'stdin': subprocess.PIPE,
                    'stdout': subprocess.PIPE,
                    'stderr': subprocess.PIPE,
                    'shell': use_shell,
                    'universal_newlines': True
                    }
        
    try:
        proc = subprocess.Popen(args = subproc_cmd, **subproc_args)
    except OSError as e:
        print('Failed to execute command: %s' % subproc_cmd)
        print(errno.errorcode[e.errno])
        print(os.strerror(e.errno))
        raise e

    return_code = proc.wait()
    stdout_lines = proc.stdout.readlines() 
    stderr_lines = proc.stderr.readlines() 
    for line in stdout_lines:
        sys.stdout.write(line)
    for line in stderr_lines:
        sys.stderr.write(line)
    if len(output) > 0:
        output_file = open(output, "a")
        for line in stdout_lines:
            output_file.write(line + '\n')
        for line in stderr_lines:
            output_file.write(line + '\n')
        output_file.close()
        
    if return_code != 0:
        sys.stderr.write('Failed to execute pdf command: {}\n'.format(subproc_cmd))
        sys.stderr.write('return code = {}\n'.format(return_code))
        
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

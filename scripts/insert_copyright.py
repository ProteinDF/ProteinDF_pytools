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
insert authorize info to source files
"""

import os
import sys
import shutil
import argparse
import logging
import re

def main():
    # parse args
    parser = argparse.ArgumentParser(description='insert copyright info into the top of file')
    parser.add_argument('source_files',
                        nargs='+',
                        help='source files')
    parser.add_argument('-i', '--insert',
                        nargs=1,
                        default='copyright.txt',
                        help='copyright filepath to insert')
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
        
    # setting
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    source_files = args.source_files
    copyright_filepath = args.insert[0]

    for source_file in source_files:
        logging.debug('processing...: {}'.format(source_file))
        (basename, ext) = os.path.splitext(source_file)

        data = None
        if ext in ['.h', '.cpp', 'cxx']:
            data = prepare_cpp_source(source_file, copyright_filepath)
        elif ext == '.py':
            data = prepare_py_source(source_file, copyright_filepath)

        if data:
            shutil.copyfile(source_file, source_file + '.orig')
            OUT_FILE = open(source_file, 'w')
            for l in data:
                l = l.rstrip('\n')
                OUT_FILE.write(l + '\n')
            OUT_FILE.close()

def get_copyright_lines(copyright_filepath, comment_out):
    logging.debug('open copyright file: {}'.format(copyright_filepath))
    INS_FILE = open(copyright_filepath, 'r')
    lines = INS_FILE.readlines()
    INS_FILE.close()

    for i in range(len(lines)):
        lines[i] = comment_out + ' ' + lines[i].rstrip('\n')
    
    return lines
    
def prepare_cpp_source(in_path, copyright_filepath):
    copyright_lines = get_copyright_lines(copyright_filepath, '//')
    re_begin_copyright_block = re.compile('^// Copyright')

    IN_FILE = open(in_path, 'r')
    
    # check if the copyright block has been inserted
    found_copyright_block = False
    lines = IN_FILE.readlines()
    copyright_block_start_line = -1
    copyright_block_end_line = -1
    for i, line in enumerate(lines):
        line = line.rstrip('\n')
        #print(line)
        
        begin_copyright_block_result= re_begin_copyright_block.match(line)
        #print(begin_copyright_block_result)
        if begin_copyright_block_result:
            assert(found_copyright_block == False)
            found_copyright_block = True
            copyright_block_start_line = i
            continue

        if found_copyright_block:
            if line[0:2] != '//':
                # end of the block
                copyright_block_end_line = i
                break
                
    # remove previous block
    new_lines = []
    for i, line in enumerate(lines):
        line = line.rstrip('\n')
        
        if copyright_block_start_line <= i < copyright_block_end_line:
            continue
        else:
            new_lines.append(line)

    lines = new_lines
            
    # insert block
    answer = []
    answer.extend(copyright_lines)
    answer.extend(lines)

    IN_FILE.close()
    return answer
            
def prepare_py_source(in_path, copyright_filepath):
    copyright_lines = get_copyright_lines(copyright_filepath, '#')
    re_begin_copyright_block = re.compile('^# Copyright')

    IN_FILE = open(in_path, 'r')
    
    # check if the copyright block has been inserted
    found_copyright_block = False
    lines = IN_FILE.readlines()
    copyright_block_start_line = -1
    copyright_block_end_line = -1
    for i, line in enumerate(lines):
        line = line.rstrip('\n')
        #print(line)
        
        begin_copyright_block_result= re_begin_copyright_block.match(line)
        #print(begin_copyright_block_result)
        if begin_copyright_block_result:
            assert(found_copyright_block == False)
            found_copyright_block = True
            copyright_block_start_line = i
            continue

        if found_copyright_block:
            if line[0:2] != '#':
                # end of the block
                copyright_block_end_line = i
                break
                
    # remove previous block
    new_lines = []
    for i, line in enumerate(lines):
        line = line.rstrip('\n')
        
        if copyright_block_start_line <= i < copyright_block_end_line:
            continue
        else:
            new_lines.append(line)

    lines = new_lines
            
    # insert block
    answer = []
    #  processing shebang
    if lines[0][0:2] == '#!':
        answer.append(lines[0])
        lines.pop(0)
        if lines[0][0:6] == '# -*- ':
            answer.append(lines[0])
            lines.pop(0)
        answer.append('\n')
    
    answer.extend(copyright_lines)
    answer.extend(lines)

    IN_FILE.close()
    return answer
            
if __name__ == '__main__':
    main()
    

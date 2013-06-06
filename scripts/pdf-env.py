#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
report machine environments
"""

import sys
import os.path
import argparse
import logging

def main():
    parser = argparse.ArgumentParser(description='report environments')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    verbose = args.verbose
    logging.basicConfig()
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    else:
        logging.basicConfig(level=logging.WARNING)

        
    # 環境変数の表示
    for k, v in os.environ.items():
        print("{key} : {value}".format(key=k, value=v))
        

if __name__ == '__main__':
    main()

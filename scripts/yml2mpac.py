#!/usr/bin/env python
# -*- coding: utf-8 -*-
"""
convert YAML file to MsgPack.
"""

import sys
import argparse
import yaml
try:
    import msgpack
except:
    import msgpack_pure as msgpack

    
def main():
    # initialize
    parser = argparse.ArgumentParser()
    parser.add_argument('YAML_FILE',
                        nargs=1,
                        help='YAML file')
    parser.add_argument('MPAC_FILE',
                        nargs=1,
                        help='message pack file')
    args = parser.parse_args()

    yaml_path = args.YAML_FILE[0]
    mpac_path = args.MPAC_FILE[0]

    fin = open(yaml_path, "rb")
    contents = fin.read()
    fin.close()

    contents = contents.decode('utf8') # for Japanese
    yaml_data = yaml.load(contents)

    mpac_data = msgpack.packb(yaml_data)

    fout = open(mpac_path, "wb")
    fout.write(mpac_data)
    fout.close()

if __name__ == '__main__':
    main()

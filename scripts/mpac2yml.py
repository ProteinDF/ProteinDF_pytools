#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
import yaml
try:
    import msgpack
except:
    import msgpack_pure as msgpack

    
def main():
    # initialize
    parser = argparse.ArgumentParser(description='display file formatted by MsgPack using YAML.')
    parser.add_argument('FILE',
                        nargs=1,
                        help='Amber PDB file')
    args = parser.parse_args()

    file_path = args.FILE[0]

    f = open(file_path, "rb")
    contents = f.read()
    data = msgpack.unpackb(contents)
    f.close()

    yaml_str = yaml.dump(data,
                         encoding='utf8',
                         allow_unicode=True,
                         default_flow_style=False)
    print(yaml_str)

if __name__ == '__main__':
    main()

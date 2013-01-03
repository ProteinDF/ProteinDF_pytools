#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import bridge


def main():
    # parse args
    parser = argparse.ArgumentParser(description='print DB(sqlite3) file')
    parser.add_argument('FILE',
                        nargs=1,
                        help='DB(SQLite3) file')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default = False)
    args = parser.parse_args()

    # setting
    db_file = args.FILE[0]
    verbose = args.verbose

    db = bridge.DbManager(db = db_file, sql_debugout = verbose)
    print(db)
    
    # end

if __name__ == '__main__':
    main()

    

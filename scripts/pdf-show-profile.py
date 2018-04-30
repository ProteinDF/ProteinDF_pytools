#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import pstats

def main():
    # parse args
    parser = argparse.ArgumentParser(description='[expert] show profile from a stats file')
    parser.add_argument("stats_file",
                        nargs=1,
                        help="stats file")
    args = parser.parse_args()
    
    stats_file = args.stats_file[0]
    ps = pstats.Stats(stats_file)

    ps.sort_stats('cumulative').print_stats()
    ps.sort_stats('time').print_stats()
    

if __name__ == "__main__":
    main()
    

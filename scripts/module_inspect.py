#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import pprint
import inspect
import pdf

def main():
    pp = pprint.PrettyPrinter()

    obj = inspect.getmembers(pdf, inspect.isclass)
    pp.pprint(obj)
    
if __name__ == '__main__':
    main()
    

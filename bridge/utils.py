#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re

def sort_nicely(l):
    """
    Sort the given list in the way that humans expect.
    ref: http://www.codinghorror.com/blog/2007/12/sorting-for-humans-natural-sort-order.html
    """
    convert = lambda text: int(text) if text.isdigit() else text
    alphanum_key = lambda key: [ convert(c) for c in re.split('([0-9]+)', key) ]
    l.sort( key=alphanum_key )

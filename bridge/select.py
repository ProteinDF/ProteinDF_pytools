#!/usr/bin/env python
# -*- coding: utf-8 -*-

import re

class Select(object):
    def is_match(self, obj):
        return False


class Select_Name(Select):
    def __init__(self, query):
        assert(isinstance(query, str) == True)
        self.query = query
    
    def is_match(self, obj):
        answer = False
        name = obj.name.strip().rstrip()
        if (name == self.query):
            answer = True
        return answer


class Select_Regex(Select):
    def __init__(self, query):
        assert(isinstance(query, str) == True)
        self.query = query
        self.regex = re.compile(query)

    def is_match(self, obj):
        answer = False
        path = obj.get_path()
        if (self.regex.search(path) != None):
            print("path=[%s] regex=[%s]" % (path, self.query))
            answer = True
        return answer


class Select_Atom(Select):
    def __init__(self, query):
        assert(isinstance(query, str) == True)
        self.query = query.upper()

    def is_match(self, obj):
        answer = False
        if (isinstance(obj, Atom) == True):
            symbol = obj.get_symbol().upper()
            if (symbol == self.query):
                answer = True
        return answer


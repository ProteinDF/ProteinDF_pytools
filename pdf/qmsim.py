#!/usr/bin/env python
# -*- coding: utf-8 -*-

class QmSim(object):
    """
    QmSimulatorの基本クラス
    """
    def __init__(self):
        self._data = {}

    # properies ================================================================
    # method
    def _get_method(self):
        return self._data.get('method', 'hf')
    def _set_method(self, method):
        self._data['method'] = str(method)
    method = property(_get_method, _set_method)


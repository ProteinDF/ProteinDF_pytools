#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2014 The ProteinDF development team.
# see also AUTHORS and README if provided.
#
# This file is a part of the ProteinDF software package.
#
# The ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# The ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.


import numpy
import matplotlib
matplotlib.use('Agg')

import matplotlib.pyplot as plt
import matplotlib.cm as cm

import math
import types

import logging
logger = logging.getLogger(__name__)

import proteindf_bridge as bridge


class DfGraph(object):
    def __init__(self):
        plt.clf()
        self._fig = None
        self._ax = None

        self._make_default_subplots()


    def _make_default_subplots(self):
        self._fig = plt.figure()
        self._ax = self._fig.add_subplot(111)


    # ==========================================================================
    # property
    # ==========================================================================
    def _set_title(self, title):
        self._title = str(title)
    def _get_title(self):
        if not '_title' in self.__dict__:
            self._title = ''
        return self._title
    title = property(_get_title, _set_title)

    def _set_xlabel(self, xlabel):
        self._xlabel = str(xlabel)
    def _get_xlabel(self):
        if not '_xlabel' in self.__dict__:
            self._xlabel = ''
        return self._xlabel
    xlabel = property(_get_xlabel, _set_xlabel)

    def _set_ylabel(self, ylabel):
        self._ylabel = str(ylabel)
    def _get_ylabel(self):
        if not '_ylabel' in self.__dict__:
            self._ylabel = ''
        return self._ylabel
    ylabel = property(_get_ylabel, _set_ylabel)

    def _set_is_draw_grid(self, TF):
        self._is_draw_grid = bool(TF)
    def _get_is_draw_grid(self):
        if not '_is_draw_grid' in self.__dict__:
            self._is_draw_grid = False
        return self._is_draw_grid
    is_draw_grid = property(_get_is_draw_grid, _set_is_draw_grid)

    def _set_xmin(self, xmin):
        self._xmin = xmin
    def _get_xmin(self):
        if not '_xmin' in self.__dict__:
            self._xmin = None
        return self._xmin
    xmin = property(_get_xmin, _set_xmin)

    def _set_xmax(self, xmax):
        self._xmax = xmax
    def _get_xmax(self):
        if not '_xmax' in self.__dict__:
            self._xmax = None
        return self._xmax
    xmax = property(_get_xmax, _set_xmax)

    def _set_ymin(self, ymin):
        self._ymin = ymin
    def _get_ymin(self):
        if not '_ymin' in self.__dict__:
            self._ymin = None
        return self._ymin
    ymin = property(_get_ymin, _set_ymin)

    def _set_ymax(self, ymax):
        self._ymax = ymax
    def _get_ymax(self):
        if not '_ymax' in self.__dict__:
            self._ymax = None
        return self._ymax
    ymax = property(_get_ymax, _set_ymax)

    def _set_aspect(self, value):
        self._aspect = value
    def _get_aspect(self):
        if not '_aspect' in self.__dict__:
            self._aspect = None
        return self._aspect
    aspect = property(_get_aspect, _set_aspect)

    def _set_xticks(self, ticks):
        assert(isinstance(ticks, types.ListType))
        self._xticks = ticks
    def _get_xticks(self):
        if not '_xticks' in self.__dict__:
            self._xticks = None
        return self._xticks
    xticks = property(_get_xticks, _set_xticks)

    def _set_yticks(self, ticks):
        assert(isinstance(ticks, types.ListType))
        self._yticks = ticks
    def _get_yticks(self):
        if not '_yticks' in self.__dict__:
            self._yticks = None
        return self._yticks
    yticks = property(_get_yticks, _set_yticks)

    def _set_xticklabels(self, labels):
        self._xticklabels = labels
    def _get_xticklabels(self):
        if not '_xticklabels' in self.__dict__:
            self._xticklabels = None
        return self._xticklabels
    xticklabels = property(_get_xticklabels, _set_xticklabels)

    def _set_yticklabels(self, labels):
        self._yticklabels = labels
    def _get_yticklabels(self):
        if not '_yticklabels' in self.__dict__:
            self._yticklabels = None
        return self._yticklabels
    yticklabels = property(_get_yticklabels, _set_yticklabels)

    def _set_legend(self, yn):
        self._legend = bool(yn)
    def _get_legend(self):
        if not '_legend' in self.__dict__:
            self._legend = False
        return self._legend
    legend = property(_get_legend, _set_legend)

    # ==========================================================================
    # method
    # ==========================================================================
    def _draw(self):
        self._ax.use_sticky_edges = False
        self._draw_data()

        self._ax.set_title(self.title)
        self._ax.set_xlabel(self.xlabel)
        self._ax.set_ylabel(self.ylabel)
        self._ax.grid(self.is_draw_grid)

        self._ax.set_xlim(left = self.xmin,
                          right = self.xmax)
        self._ax.set_ylim(bottom = self.ymin,
                          top = self.ymax)
        # self._ax.relim(visible_only=True)
        self._ax.autoscale_view(tight=True)

        if self.xticks != None:
            self._ax.set_xticks(self.xticks)
        if self.yticks != None:
            self._ax.set_yticks(self.yticks)
        if self.xticklabels != None:
            self._ax.set_xticklabels(self.xticklabels)
        if self.yticklabels != None:
            self._ax.set_yticklabels(self.yticklabels)

        if self.aspect != None:
            self._ax.set_aspect(self.aspect)

        if self.legend == True:
            self._ax.legend()



    def _draw_data(self):
        pass

    def save(self, path):
        self._draw()
        self._fig.savefig(path)


class DfTotalEnergyHistGraph(DfGraph):
    def __init__(self):
        DfGraph.__init__(self)
        self.title = 'total energy history'
        self.xlabel = 'iteration'
        self.ylabel = 'energy / a.u.'
        self.is_draw_grid = True

    def load_data(self, path):
        f = open(path)
        lines = f.readlines();
        f.close()

        self._x = [];
        self._y = [];
        for line in lines:
            items = line.strip().split(',')
            self._x.append(float(items[0]))
            self._y.append(float(items[1]))

    def _draw_data(self):
        # see. http://matplotlib.sourceforge.net/api/pyplot_api.html#matplotlib.pyplot.plot
        self._ax.plot(self._x, self._y, 'bo-')

        # self._ax.legend(('total energy'), loc='best')


class DfEnergyLevelHistoryGraphH(DfGraph):
    """
    エネルギー準位を水平線で描画します
    """
    def __init__(self):
        DfGraph.__init__(self)
        self.title = 'Energy Level'
        self.xlabel = 'iteration'
        self.ylabel = 'energy / eV'
        self.is_draw_grid = True
        self.xmin = None
        self.xmax = None
        self.ymin = -20.0
        self.ymax =   5.0

        self._HOMO_level = -1
        self._LUMO_level = -1
        self._max_itr = 0
        self._select_iterations = []

    def set_HOMO_level(self, level):
        self._HOMO_level = int(level)

    def set_LUMO_level(self, level):
        self._LUMO_level = int(level)

    def load_data(self, path):
        f = open(path)
        lines = f.readlines();
        f.close()

        self._data = []
        for line in lines:
            line = line.strip()
            if len(line) == 0:
                continue
            (itr, level, value) = line.strip().split(',')
            itr = int(itr)
            self._max_itr = max(self._max_itr, itr)
            self._data.append([itr, float(level), float(value)])

    def select_iterations(self, iterations):
        self._select_iterations = set(iterations)

    def _draw_data(self):
        draw_iterations = []
        if len(self._select_iterations) == 0:
            for itr in  range(self._max_itr +1):
                draw_iterations.append(itr +1)
        else:
            for order, itr in enumerate(self._select_iterations):
                draw_iterations.append(itr)

        for d in self._data:
            itr = d[0]
            level = d[1]
            value = d[2]

            if itr in draw_iterations:
                order = draw_iterations.index(itr)
            else:
                continue

            strike = False
            color = 'k'
            if level == self._HOMO_level:
                strike = True
                color = 'r'
            elif level == self._LUMO_level:
                strike = True
                color = 'b'
            self._draw_data_line(order, value,
                                 strike=strike,
                                 color=color)

    def _draw_data_line(self, order, value, strike=False, color='k'):
        if (self.ymin < value) and (value < self.ymax):
            width = 0.8
            height = 0.01
            if strike:
                width = 1.0
                height = 0.01

            left = order - width / 2.0
            bottom = value - height / 2.0

            self._ax.bar(left, height, width, bottom,
                         edgecolor=color, color=color)


class DfEnergyLevelHistoryGraphV(DfEnergyLevelHistoryGraphH):
    """
    エネルギー準位を垂直線で描画します。
    """
    def __init__(self):
        DfEnergyLevelHistoryGraphH.__init__(self)
        self.title = 'Energy Level'
        self.xlabel = 'energy / eV'
        self.ylabel = 'iteration'
        self.is_draw_grid = True
        self.xmin = -20.0
        self.xmax =   5.0
        self.ymin = None
        self.ymax = None

        self._HOMO_level = -1
        self._LUMO_level = -1

    def _draw_data_line(self, order, value, strike=False, color='k'):
        if (self.xmin < value) and (value < self.xmax):
            width = 0.01
            height = 0.8
            if strike:
                width = 0.01
                height = 1.0

            bottom = order - height / 2.0
            left = value - width / 2.0

            self._ax.bar(left, height, width, bottom,
                         edgecolor=color, color=color)

    def _draw(self):
        self._ax.tick_params(labelleft=False)
        super(DfEnergyLevelHistoryGraphV, self)._draw()


class DfPopulationGraph(DfGraph):
    def __init__(self):
        super().__init__()

    def load_data(self, path):
        self._x = []
        self._y = []
        with open(path) as f:
            for line in f:
                items = line.strip().split(',')
                self._x.append(float(items[0]))
                self._y.append(float(items[1]))

    def _draw_data(self):
        self._ax.bar(self._x, self._y)


class DfLineChart(DfGraph):
    def __init__(self):
        DfGraph.__init__(self)
        self._data = []

    def add_data(self, X, Y):
        d = (X, Y)
        self._data.append(d)

    def _draw_data(self):
        # print('data size=', len(self._data))
        # print(self._data)
        for i, (X, Y) in enumerate(self._data):
            self._ax.plot(X, Y, linewidth=1.0, linestyle="-", label='{}'.format(i))


class DfMatrixGraph(DfGraph):
    def __init__(self, matrix, is_diverging = True):
        assert(isinstance(matrix, bridge.Matrix))
        DfGraph.__init__(self)
        self._matrix = matrix
        self._is_diverging = is_diverging

        self.title = 'Matrix Value'
        self.xlabel = ''
        self.ylabel = ''
        self.is_draw_grid = True

        if self._is_diverging:
            self._cmap = cm.coolwarm
        else:
            self._cmap = cm.bone_r

    def _draw_data(self):
        data = self._matrix.data
        if not self._is_diverging:
            data = numpy.absolute(data)

        cax = self._ax.imshow(data,
                              cmap=self._cmap,
                              interpolation='none',
                              norm=None,
                              #vmax=1.0,
                              #vmin=0.0,
                              origin='upper')

        # self._ax.tick_params(axis='x', labeltop=True, labelbottom=False)
        self._ax.tick_params(axis='x', labeltop=False, labelbottom=True)

        # self._fig.colorbar(cax, ticks=[ 0, 1], shrink=0.92)
        self._fig.colorbar(cax, shrink=0.92)


class DfDistanceVsElementGraph(DfGraph):
    def __init__(self, log = False):
         DfGraph.__init__(self)
         self._log = log

         self.xlabel = 'distance'
         self.ylabel = 'magnitude'

    def _make_default_subplots(self):
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h = left_h = left+width+0.02

        rect_scatter = [left, bottom, width, height]
        rect_hist_x = [left, bottom_h, width, 0.2]
        rect_hist_y = [left_h, bottom, 0.2, height]

        self._fig = plt.figure(1)

        self._ax_scatter = plt.axes(rect_scatter)
        self._ax_hist_x = plt.axes(rect_hist_x)
        self._ax_hist_y = plt.axes(rect_hist_y)

        nullfmt = matplotlib.ticker.NullFormatter()
        #self._ax_hist_x.xaxis.set_major_formatter(nullfmt)
        #self._ax_hist_y.yaxis.set_major_formatter(nullfmt)


    def load(self, path):
        self._data_x = []
        self._data_y = []
        data_x = []
        data_y = []

        with open(path) as f:
            for line in f:
                line = line.strip()
                d = line.split(',')

                x = float(d[0])
                y = float(d[1])
                y = math.fabs(y)

                data_x.append(x)
                data_y.append(y)

        self._data_x = numpy.array(data_x, 'd')
        self._data_y = numpy.array(data_y, 'd')

    def _draw(self):
        self._draw_data()
        self._ax_scatter.set_xlabel(self.xlabel)
        self._ax_scatter.set_ylabel(self.ylabel)
        self._ax_scatter.grid(self.is_draw_grid)

    def _draw_data(self):
        self._ax_scatter.scatter(self._data_x, self._data_y)
        if self._log:
            self._ax_scatter.set_ylim(1.0E-8, 100.0)
            self._ax_scatter.set_yscale('log')

        self._draw_histx()
        self._draw_histy()


    def _draw_histx(self):
        self._ax_hist_x.hist(self._data_x)
        self._ax_hist_x.set_xlim(self._ax_scatter.get_xlim())

        (y_min, y_max) = self._ax_hist_x.get_ylim()
        range_y = y_max - y_min
        width_y = range_y / 2
        self._ax_hist_x.set_yticks(numpy.arange(y_min, y_max, width_y))


    def _draw_histy(self):
        self._ax_hist_y.hist(self._data_y, orientation='horizontal', log=self._log)
        self._ax_hist_y.set_ylim(self._ax_scatter.get_ylim())
        if self._log:
            self._ax_hist_y.set_yscale('log')

        (x_min, x_max) = self._ax_hist_y.get_xlim()
        range_x = x_max - x_min
        width_x = range_x / 2
        self._ax_hist_y.set_xticks(numpy.arange(x_min, x_max, width_x))


class DfLineChart(DfGraph):
    def __init__(self):
        DfGraph.__init__(self)
        self._matrix = matrix
        self._is_diverging = is_diverging

        self.title = 'Matrix Value'
        self.xlabel = ''
        self.ylabel = ''
        self.is_draw_grid = True

        if self._is_diverging:
            self._cmap = cm.coolwarm
        else:
            self._cmap = cm.bone_r

    def _draw_data(self):
        data = self._matrix.data
        if not self._is_diverging:
            data = numpy.absolute(data)

        cax = self._ax.imshow(data,
                              cmap=self._cmap,
                              interpolation='none',
                              norm=None,
                              #vmax=1.0,
                              #vmin=0.0,
                              origin='upper')

        # self._ax.tick_params(axis='x', labeltop=True, labelbottom=False)
        self._ax.tick_params(axis='x', labeltop=False, labelbottom=True)

        # self._fig.colorbar(cax, ticks=[ 0, 1], shrink=0.92)
        self._fig.colorbar(cax, shrink=0.92)


class DfDistanceVsElementGraph(DfGraph):
    def __init__(self, log = False):
         DfGraph.__init__(self)
         self._log = log

         self.xlabel = 'distance'
         self.ylabel = 'magnitude'

    def _make_default_subplots(self):
        left, width = 0.1, 0.65
        bottom, height = 0.1, 0.65
        bottom_h = left_h = left+width+0.02

        rect_scatter = [left, bottom, width, height]
        rect_hist_x = [left, bottom_h, width, 0.2]
        rect_hist_y = [left_h, bottom, 0.2, height]

        self._fig = plt.figure(1)

        self._ax_scatter = plt.axes(rect_scatter)
        self._ax_hist_x = plt.axes(rect_hist_x)
        self._ax_hist_y = plt.axes(rect_hist_y)

        nullfmt = matplotlib.ticker.NullFormatter()
        #self._ax_hist_x.xaxis.set_major_formatter(nullfmt)
        #self._ax_hist_y.yaxis.set_major_formatter(nullfmt)


    def load(self, path):
        self._data_x = []
        self._data_y = []
        data_x = []
        data_y = []

        with open(path) as f:
            for line in f:
                line = line.strip()
                d = line.split(',')

                x = float(d[0])
                y = float(d[1])
                y = math.fabs(y)

                data_x.append(x)
                data_y.append(y)

        self._data_x = numpy.array(data_x, 'd')
        self._data_y = numpy.array(data_y, 'd')

    def _draw(self):
        self._draw_data()
        self._ax_scatter.set_xlabel(self.xlabel)
        self._ax_scatter.set_ylabel(self.ylabel)
        self._ax_scatter.grid(self.is_draw_grid)

    def _draw_data(self):
        self._ax_scatter.scatter(self._data_x, self._data_y)
        if self._log:
            self._ax_scatter.set_ylim(1.0E-8, 100.0)
            self._ax_scatter.set_yscale('log')

        self._draw_histx()
        self._draw_histy()


    def _draw_histx(self):
        self._ax_hist_x.hist(self._data_x)
        self._ax_hist_x.set_xlim(self._ax_scatter.get_xlim())

        (y_min, y_max) = self._ax_hist_x.get_ylim()
        range_y = y_max - y_min
        width_y = range_y / 2
        self._ax_hist_x.set_yticks(numpy.arange(y_min, y_max, width_y))


    def _draw_histy(self):
        self._ax_hist_y.hist(self._data_y, orientation='horizontal', log=self._log)
        self._ax_hist_y.set_ylim(self._ax_scatter.get_ylim())
        if self._log:
            self._ax_hist_y.set_yscale('log')

        (x_min, x_max) = self._ax_hist_y.get_xlim()
        range_x = x_max - x_min
        width_x = range_x / 2
        self._ax_hist_y.set_xticks(numpy.arange(x_min, x_max, width_x))



# ******************************************************************************

class LineChart(object):
    def set_data(self, x, y):
        self.__x = x
        self.__y = y

    def savefig(self, fname):
        plt.figure()
        plt.subplot(111)

        # clear
        plt.clf()

        # line
        line = plt.plot(self.__x,
                        self.__y)
        # range
        is_tight = True
        range = plt.axis()
        if '__xmin' in self.__dict__:
            range[0] = self.__xmin
            is_tight = False
        if '__xmax' in self.__dict__:
            range[1] = self.__xmax
            is_tight = False
        if '__ymin' in self.__dict__:
            range[2] = self.__ymin
            is_tight = False
        if '__ymax' in self.__dict__:
            range[3] = self.__ymax
            is_tight = False
        if is_tight:
            plt.axis('tight')
        else:
            plt.axis(range)

        plt.savefig(fname)


class EnergyLevelChart(object):
    __values = {}

    def set_data(self, itr, values):
        self.__data[itr] = values

    def savefig(self, fname):
        plt.figure()
        plt.subplot(111)

        # clear
        plt.clf()

        # line
        for itr in self.__data.keys():
            for value in self.__data[itr]:
                bottom = value
                width = 0.8
                height = 0
                left = itr - 0.45
                line = plt.barh(bottom, width, height, left)

        # range
        is_tight = True
        range = plt.axis()
        if '__xmin' in self.__dict__:
            range[0] = self.__xmin
            is_tight = False
        if '__xmax' in self.__dict__:
            range[1] = self.__xmax
            is_tight = False
        if '__ymin' in self.__dict__:
            range[2] = self.__ymin
            is_tight = False
        if '__ymax' in self.__dict__:
            range[3] = self.__ymax
            is_tight = False
        if is_tight:
            plt.axis('tight')
        else:
            plt.axis(range)

        plt.savefig(fname)



class Graph(object):
    colors = ['blue', 'green', 'red', 'cyan', 'magenta',
              'yellow']

    def __init__(self, size = None):
        """
        sizeはinchで指定する(width x height)
        """
        # clear the figure
        #self.fig.clf()

    def file_out(self, file_path):
        self.__prepare()
        self.fig.savefig(file_path)


    def set_xrange(self, min_value, max_value):
        self.param['min_x'] = min_value
        self.param['max_x'] = max_value


    def set_yrange(self, min_value, max_value):
        self.param['min_y'] = min_value
        self.param['max_y'] = max_value


    def set_xscale(self, value):
        assert(value in ['linear', 'log', 'symlog'])
        self.param['xscale'] = value


    def set_yscale(self, value):
        assert(value in ['linear', 'log', 'symlog'])
        self.param['yscale'] = value


    def set_xlabel(self, label):
        self.param['xlabel'] = label


    def set_ylabel(self, label):
        self.param['ylabel'] = label


    def draw_grid(self, yn):
        assert((yn == True) or (yn == False))
        self.param['draw_grid'] = yn


    def draw_legend(self, yn):
        assert((yn == True) or (yn == False))
        self.param['draw_legend'] = yn


    def plot_basis(self, type = 's', coef = [], exponent = [], color = 'blue'):
        """
        print Gaussian
        """
        assert(len(coef) == len(exponent))

        max_value = 5
        x = []
        for i in range(max_value * 100):
            v = 0.01 * float(i)
            x.append(v)

        y = []
        if (type == 's'):
            for i in x:
                i = float(i)
                numOfItems = len(coef)
                v = 0.0
                for t in range(numOfItems):
                    tmp = float(exponent[t]) * (i * i)
                    if (fabs(tmp) > log(10000.0)):
                        continue
                    tmp2 = exp(- tmp)
                    v += float(coef[t]) * tmp2
                y.append(v)
        elif (type == 'p'):
            for i in x:
                i = float(i)
                numOfItems = len(coef)
                v = 0.0
                for t in range(numOfItems):
                    tmp = float(exponent[t]) * (i * i)
                    if (fabs(tmp) > log(10000.0)):
                        continue
                    tmp2 = exp(- tmp)
                    v += float(coef[t]) * tmp2 * i
                y.append(v)
        elif (type == 'd'):
            for i in x:
                i = float(i)
                numOfItems = len(coef)
                v = 0.0
                for t in range(numOfItems):
                    tmp = float(exponent[t]) * (i * i)
                    if (fabs(tmp) > log(10000.0)):
                        continue
                    tmp2 = exp(- tmp)
                    v += float(coef[t]) * tmp2 * i * i
                y.append(v)

        self.ax.plot(x, y, color)


    def __prepare(self):
        plt.figure()
        plt.subplot(111)

        # clear
        plt.clf()

        # range
        range = plt.axis()
        if self.__xmin:
            range[0] = self.__xmin
        if self.__xmax:
            range[1] = self.__xmax
        if self.__ymin:
            range[2] = self.__ymin
        if self.__ymax:
            range[3] = self.__ymax

        # locator
        if (self.param.get('xlocator', True) == False):
            nullLocator = NullLocator()
            self.ax.xaxis.set_major_locator(nullLocator)
            self.ax.xaxis.set_minor_locator(nullLocator)
        if (self.param.get('ylocator', True) == False):
            nullLocator = NullLocator()
            self.ax.yaxis.set_major_locator(nullLocator)
            self.ax.yaxis.set_minor_locator(nullLocator)

        # axis label
        if (self.param.has_key('xlabel') == True):
            self.ax.set_xlabel(self.param['xlabel'])
        if (self.param.has_key('ylabel') == True):
            self.ax.set_ylabel(self.param['ylabel'])

        # grid
        if (self.param.setdefault('draw_grid', False) == True):
            self.ax.grid(True)

        # legend
        if (self.param.setdefault('draw_legend', False) == True):
            self.ax.legend(loc='best')


class GraphEnergyLevelHistory(Graph):
    def __init__(self, size = None):
        Graph.__init__(self, size)

    def plot_energy_levels(self, iteration, levels, homo_level = -1):
        """
        levels is the list of energy level.
        """
        if (levels == None):
            return

        min_y = self.param['min_y']
        max_y = self.param['max_y']

        # this value is for negligence
        last_value = min_y
        negligence_range = (max_y - min_y) * 0.001 # 0.1%

        for level, value in enumerate(levels):
            # draw HOMO
            if (level == homo_level):
                self.ax.broken_barh([(iteration -0.50, 1.0)], (value, 0),
                                    edgecolors='red', facecolors='red')
                continue

            # negligence
            if ((value < min_y) or (max_y < value)):
                continue
            if (abs(last_value - value) < negligence_range):
                continue
            else:
                last_value = value

            # register data
            self.ax.broken_barh([(iteration -0.45, 0.9)], (value, 0),
                                edgecolors='black', facecolors='black')


class GraphEnergyLevelSingle(Graph):
    def __init__(self, size = None):
        Graph.__init__(self, size)

    def plot_energy_levels(self, levels, homo_level = -1):
        """
        levels is the list of energy level.
        """
        self.param['ylocator'] = False

        min_y = self.param['min_y']
        max_y = self.param['max_y']

        # this value is for negligence
        last_value = min_y
        negligence_range = (max_y - min_y) * 0.001 # 0.1%

        # plot
        for level, value in enumerate(levels):
            # draw HOMO
            if (level == homo_level):
                self.ax.broken_barh([(value, 0)], (1.0 -0.50, 1.0),
                                    edgecolors='red', facecolors='red'
                                    )

            # negligence
            if ((value < min_y) or (max_y < value)):
                continue
            if (abs(last_value - value) < negligence_range):
                continue
            else:
                last_value = value

            # register data
            self.ax.broken_barh([(value, 0)], (1.0 -0.50, 1.0),
                                edgecolors='black', facecolors='black'
                                )


class GraphTotalEnergyHistory(Graph):
    def __init__(self, size = None):
        Graph.__init__(self, size)

    def plot(self, values, label = "", marker = '+', linestyle = '-'):
        """
        plot total energy history data
        """
        x = []
        y = []
        for i in range(len(values)):
            if (values[i] != None):
                x.append(i)
                y.append(values[i])
        self.ax.plot(x, y, label = label, marker = marker, linestyle = linestyle)


class GraphConvergenceCheck(Graph):
    def __init__(self, size = None):
        Graph.__init__(self, size)

    def plot(self, values, label = "", marker = '+', linestyle = '-'):
        """
        plot convergence history data
        """
        x = []
        y = []
        for i in range(len(values)):
            if (values[i] != None):
                x.append(i)
                y.append(values[i])
        self.ax.plot(x, y, label = label, marker = marker, linestyle = linestyle)

    def plotLog(self, values, label = "", marker = '+', linestyle = '-'):
        """
        plot convergence history data (log)
        """
        x = []
        y = []
        for i in range(len(values)):
            if (values[i] != None):
                x.append(i)
                y.append(values[i])
        self.ax.semilogy(x, y, label = label, marker = marker, linestyle = linestyle)

class BarGraph(DfGraph):
    def __init__(self):
        DfGraph.__init__(self)
        #self.ax = pylab.subplot(121)
        #params = {
        #    'legend.fontsize' : 9
        #    }
        #pylab.rcParams.update(params)

    def _draw_data(self):
        pass

    def plot(self, data, labels = None):
        # prepare
        rows = len(data)
        cols = len(data[0]) if rows else 0
        for i in range(rows):
            assert cols == len(data[i])

        colors = get_colours(cols)
        colors.reverse()

        for row in range(rows):
            left = row
            y_offset = 0.0
            for index, value in enumerate(data[row]):
                self.ax.bar(left, value,
                            bottom=y_offset,
                            color=colors[index])
                y_offset += value

        if labels:
            #self.ax.legend(labels, bbox_to_anchor=(1.05, 1), loc=2)
            self.ax.legend(labels, loc='best')

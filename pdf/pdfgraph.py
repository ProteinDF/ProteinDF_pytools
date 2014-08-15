#!/usr/bin/env python
# -*- coding: utf-8 -*-

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import types

class DfGraph(object):
    def __init__(self):
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

    # ==========================================================================
    # method
    # ==========================================================================
    def _draw(self):
        self._draw_data()

        self._ax.set_title(self.title)
        self._ax.set_xlabel(self.xlabel)
        self._ax.set_ylabel(self.ylabel)
        self._ax.grid(self.is_draw_grid)

        self._ax.set_xlim(left = self.xmin,
                          right = self.xmax)
        self._ax.set_ylim(bottom = self.ymin,
                          top = self.ymax)
        self._ax.autoscale_view(tight=True)

        if self.xticks != None:
            self._ax.set_xticks(self.xticks)
        if self.yticks != None:
            self._ax.set_yticks(self.yticks)

        if self.aspect != None:
            self._ax.set_aspect(self.aspect)
        
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

        self._ax.legend(('total energy', loc='best'))

        
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
        self._max_itr = 0
        self._select_iterations = []
        
    def set_HOMO_level(self, level):
        self._HOMO_level = int(level)
        
    def load_data(self, path):
        f = open(path)
        lines = f.readlines();
        f.close()
        
        self._data = []
        for line in lines:
            (itr, level, value) = line.strip().split(',')
            itr = int(itr)
            self._max_itr = max(self._max_itr, itr)
            self._data.append([itr, float(level), float(value)])

    def select_iterations(self, iterations):
        self._select_iterations = set(iterations)
            
    def _draw_data(self):
        itr_vs_list = []
        if len(self._select_iterations) == 0:
            itr_vs_list = range(self._max_itr +1)
        else:
            itr_vs_list = [ 0 for x in range(self._max_itr +1) ]
            for order, itr in enumerate(self._select_iterations):
                itr_vs_list[itr] = order +1
                
        for d in self._data:
            itr = d[0]
            level = d[1]
            value = d[2]

            order = itr_vs_list[itr]
            if order == 0:
                continue
            
            self._draw_data_line(order, value, (level == self._HOMO_level))

    def _draw_data_line(self, order, value, is_HOMO):
        if (self.ymin < value) and (value < self.ymax):
            width = 0.8
            height = 0.01
            c = 'k'
            if is_HOMO:
                width = 1.0
                height = 0.01
                c = 'r'

            left = order - width / 2.0
            bottom = value - height / 2.0

            self._ax.bar(left, height, width, bottom,
                         edgecolor=c, color=c)


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

    def _draw_data_line(self, order, value, is_HOMO):
        if (self.xmin < value) and (value < self.xmax):
            width = 0.01
            height = 0.8
            c = 'k'
            if is_HOMO:
                width = 0.01
                height = 1.0
                c = 'r'

            bottom = order - height / 2.0
            left = value - width / 2.0

            self._ax.bar(left, height, width, bottom,
                         edgecolor=c, color=c)

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

class BarGraph(Graph):
    def __init__(self, size = None):
        Graph.__init__(self, size)
        #self.ax = pylab.subplot(121)
        params = {
            'legend.fontsize' : 9
            }
        pylab.rcParams.update(params)
        
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
            
        

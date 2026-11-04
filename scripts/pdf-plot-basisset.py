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

import os
import sys
from types import *
import argparse
import yaml

import math
import numpy
import matplotlib as mpl
#import matplotlib.pyplot as plt

import proteindf_bridge as bridge
import proteindf_tools as pdf

import logging


class _FlowMapping(dict):
    """Mapping serialized in YAML inline flow style."""

    pass


yaml.add_representer(
    _FlowMapping,
    lambda dumper, data: dumper.represent_mapping(
        "tag:yaml.org,2002:map", data, flow_style=True
    ),
    Dumper=yaml.SafeDumper,
)


def _prepare_plot_data_for_yaml(plot_data):
    if not isinstance(plot_data, dict) or "series" not in plot_data:
        return plot_data
    series_list = plot_data.get("series")
    if not isinstance(series_list, list):
        return plot_data
    need_copy = False
    for s in series_list:
        if isinstance(s, dict) and "points" in s:
            pts = s.get("points")
            if isinstance(pts, list) and any(not isinstance(p, _FlowMapping) for p in pts):
                need_copy = True
                break
    if not need_copy:
        return plot_data
    new_data = dict(plot_data)
    new_series = []
    for s in series_list:
        if isinstance(s, dict) and "points" in s:
            new_s = dict(s)
            new_s["points"] = [
                p if isinstance(p, _FlowMapping) else _FlowMapping(p)
                for p in s.get("points", [])
            ]
            new_series.append(new_s)
        else:
            new_series.append(s)
    new_data["series"] = new_series
    return new_data


def make_basisset_plot_data(basisset_name, bs, X, ys=None):
    if ys is None:
        ys = [func(bs[i], X) for i in range(len(bs))]
    series = []
    for cgto_id, Y in enumerate(ys):
        cgto = bs[cgto_id]
        points = [
            _FlowMapping({"x": float(x), "y": float(y)})
            for x, y in zip(X, Y)
        ]
        series.append(
            {
                "name": "{}".format(cgto_id),
                "shell_type": str(cgto.shell_type),
                "points": points,
            }
        )
    return {
        "schema": "pdf-plot/1",
        "kind": "basisset",
        "title": basisset_name,
        "source": {
            "tool": "pdf-plot-basisset.py",
            "basisset_name": basisset_name,
            "num_of_cgtos": len(bs),
        },
        "axes": {
            "x": {
                "label": "distance / a.u.",
                "unit": "a.u.",
            },
            "y": {
                "label": "intensity",
                "unit": None,
            },
        },
        "series": series,
    }


def save_basisset_yaml(plot_data, output_path):
    data = _prepare_plot_data_for_yaml(plot_data)
    content = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)


def main():
    # parse args
    parser = argparse.ArgumentParser(description='plot basis set')
    parser.add_argument('basisset_name',
                        nargs=1,
                        help='basisset name')
    parser.add_argument('--yaml',
                        nargs=1,
                        default=None,
                        help='output yaml path')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    parser.add_argument('-D', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    basisset_name = args.basisset_name[0]
    yaml_path = args.yaml[0] if args.yaml else None

    N = 100
    X = numpy.linspace(0, 10, N, endpoint=True)

    basis2 = pdf.Basis2()
    bs = basis2.get_basisset(basisset_name)
    if len(bs) == 0:
        sys.stderr.write('unknown basisset: {}\n'.format(basisset_name))
        sys.exit(1)

    num_of_CGTOs = len(bs)
    ys = []
    for cgto_id in range(num_of_CGTOs):
        CGTO = bs[cgto_id]
        Y = func(CGTO, X)
        ys.append(Y)

    if yaml_path:
        plot_data = make_basisset_plot_data(basisset_name, bs, X, ys=ys)
        save_basisset_yaml(plot_data, yaml_path)
    else:
        chart = pdf.DfLineChart()
        for Y in ys:
            chart.add_data(X, Y)

        chart.ymin = 0.0
        chart.ymax = 2.0
        # chart.draw_legend(True)

        chart.title = basisset_name
        chart.xlabel = 'distance / a.u.'
        chart.ylabel = 'intensity'

        chart.save("basis_{}.png".format(basisset_name))

    # debug
    # for cgto_id in range(num_of_CGTOs):
    #    CGTO = bs[cgto_id]
    #    CGTO_norm = CGTO.normalize()
    #    print("shell_type: {}".format(CGTO.shell_type))
    #    for pgto_id in range(len(CGTO)):
    #        PGTO = CGTO[pgto_id]
    #        PGTO_norm = PGTO.normalize(CGTO.shell_type)
    #        print("{} {} -> {} ".format(PGTO.exp, PGTO.coef,
    #                                    CGTO_norm * PGTO_norm * PGTO.coef))


def func(CGTO, X):
    shell_type = CGTO.shell_type.upper()
    pow_x = 1
    if shell_type == 'S':
        pow_x = 0
    elif shell_type == 'P':
        pow_x = 1
    elif shell_type == 'D':
        pow_x = 2
    elif shell_type == 'F':
        pow_x = 3
    elif shell_type == 'G':
        pow_x = 4

    CGTO_norm = CGTO.normalize()
    #CGTO_norm = 1.0
    Y = numpy.ndarray(shape=(len(X)), dtype=float)
    for i in range(len(X)):
        x = X[i]
        xx = x * x

        coef_x = 1.0
        for cycle in range(pow_x):
            coef_x *= x

        value = 0.0
        for j in range(len(CGTO)):
            PGTO_norm = CGTO[j].normalize(CGTO.shell_type)
            #PGTO_norm = 1.0
            pgto_y = (CGTO_norm * PGTO_norm *
                      CGTO[j].coef * coef_x *
                      math.exp(-1.0 * CGTO[j].exp * xx))

            value += xx * math.fabs(pgto_y)
        Y[i] = value

    return Y


if __name__ == '__main__':
    main()

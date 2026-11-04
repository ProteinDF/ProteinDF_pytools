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
import argparse
import yaml

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


def make_elevel_plot_data(source_path, itr, HOMO_level, eigvals, run_type=None):
    values = [float(e) * 27.2116 for e in eigvals]
    return {
        "schema": "pdf-plot/1",
        "kind": "elevel",
        "title": "Energy Level",
        "source": {
            "tool": "pdf-plot-elevel.py",
            "input": source_path,
            "run_type": run_type,
            "iteration": int(itr),
            "homo_level": int(HOMO_level),
        },
        "axes": {
            "x": {
                "label": "energy / eV",
                "unit": "eV",
            },
            "y": {
                "label": "",
                "unit": None,
            },
        },
        "series": [
            {
                "name": "energy_level",
                "values": values,
            }
        ],
    }


def save_elevel_yaml(plot_data, output_path):
    data = _prepare_plot_data_for_yaml(plot_data)
    content = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)



def main():
    # parse args
    parser = argparse.ArgumentParser(
        description='plot energy level(single; vertical)')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-e', '--elevel_vector',
                       nargs=1,
                       help='ProteinDF energy level vector file')
    group.add_argument(
        "-d",
        "--db",
        nargs="?",
        action="store",
        const="pdfresults.h5",
        help="ProteinDF results file (.h5 は HDF5、それ以外は SQLite(従来))",
    )
    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=None,
                        help='output graph path')
    parser.add_argument('--yaml',
                        nargs=1,
                        default=None,
                        help='output yaml path')
    parser.add_argument('-l', '--level',
                        nargs=1,
                        default=['0'],
                        help='HOMO level')
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
    yaml_path = args.yaml[0] if args.yaml else None
    output_path = args.output[0] if args.output else (None if yaml_path else 'elevel.png')
    HOMO_level = int(args.level[0])

    if args.elevel_vector:
        eigval_path = args.elevel_vector[0]

        itr = 1
        data_path = 'tmp.dat'

        eigvals = pdf.Vector()
        if verbose:
            print('load: %s' % (eigval_path))
            if output_path:
                print('output: %s' % (output_path))
            if yaml_path:
                print('yaml: %s' % (yaml_path))
            print('HOMO: %d' % (HOMO_level))
        eigvals.load(eigval_path)
        if yaml_path:
            plot_data = make_elevel_plot_data(eigval_path, itr, HOMO_level, eigvals)
            save_elevel_yaml(plot_data, yaml_path)
        if output_path:
            elevel2dat(itr, eigvals, data_path)
            plot_elevel_single(data_path, itr, HOMO_level, output_path)
    else:
        entry = pdf.open_archive(args.db)

        runtypes = entry.get_energy_level_runtypes()
        if len(runtypes) == 0:
            raise ValueError('ENERGY_LEVEL is not found in database')
        elif len(runtypes) > 1:
            raise ValueError(
                'unsupported runtypes %s: UKS/ROKS are not supported' % runtypes
            )
        runtype = runtypes[0]

        itr = entry.iterations
        HOMO_level = entry.get_HOMO_level(runtype)
        data_path = 'tmp.dat'

        eigvals = entry.get_energylevel(runtype, itr)
        if verbose:
            if output_path:
                print('output: %s' % (output_path))
            if yaml_path:
                print('yaml: %s' % (yaml_path))
            print('HOMO: %d' % (HOMO_level))
        if yaml_path:
            plot_data = make_elevel_plot_data(
                args.db, itr, HOMO_level, eigvals, run_type=runtype
            )
            save_elevel_yaml(plot_data, yaml_path)
        if output_path:
            elevel2dat(itr, eigvals, data_path)
            plot_elevel_single(data_path, itr, HOMO_level, output_path)


def elevel2dat(itr, eigvals, data_path):
    """
    エネルギー準位ベクトルファイルを
    グラフ描画用データに変換する

    itr: iteraton(int)
    eigvals: energy level(bridge.Vector object)
    data_path: output path of data
    """
    assert(isinstance(itr, int))
    assert(isinstance(eigvals, bridge.Vector))
    assert(isinstance(data_path, str))

    dat = open(data_path, 'w')
    for level, e in enumerate(eigvals):
        e *= 27.2116
        dat.write('%d, %d, % 16.10f\n' % (itr, level, e))
    dat.close()


def plot_elevel_single(data_path, itr, HOMO_level, output_path):
    """
    横向きのエネルギー準位グラフを描画する

    data_path: グラフ描画用データファイルのパス
    HOMO: 最高占有軌道(0からスタート)
    itr: iteration
    output_path:
    """
    graphV = pdf.DfEnergyLevelHistoryGraphV()
    graphV.set_HOMO_level(HOMO_level)  # option base 0
    graphV.load_data(data_path)
    graphV.select_iterations([itr])
    graphV.is_draw_grid = False
    graphV.yticks = []
    graphV.ylabel = ''
    graphV.ymin = itr - 0.6
    graphV.ymax = itr + 0.6
    graphV.aspect = 3.0

    graphV.save(output_path)


if __name__ == '__main__':
    main()

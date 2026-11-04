#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import csv
import logging
import math
import os
import sys
from types import *

import numpy
import proteindf_bridge as bridge
import proteindf_tools as pdf
import yaml
from matplotlib import use


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


def make_gvalues_plot_data(
    h5_path,
    method,
    iteration,
    HOMO_level,
    data_x_occ,
    data_y_occ,
    data_x_vir,
    data_y_vir,
):
    return {
        "schema": "pdf-plot/1",
        "kind": "gvalues",
        "title": "g-values",
        "source": {
            "tool": "pdf-plot-gvalues.py",
            "result": h5_path,
            "method": method,
            "iteration": int(iteration),
            "homo_level": int(HOMO_level),
        },
        "axes": {
            "x": {
                "label": "energy level / eV",
                "unit": "eV",
            },
            "y": {
                "label": "g-values",
                "unit": None,
            },
        },
        "series": [
            {
                "name": "occ",
                "points": [
                    _FlowMapping({"x": float(x), "y": float(y)})
                    for x, y in zip(data_x_occ, data_y_occ)
                ],
            },
            {
                "name": "vir",
                "points": [
                    _FlowMapping({"x": float(x), "y": float(y)})
                    for x, y in zip(data_x_vir, data_y_vir)
                ],
            },
        ],
    }


def save_gvalues_yaml(plot_data, output_path):
    data = _prepare_plot_data_for_yaml(plot_data)
    content = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)



def output_csv(csv_data, path):
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(csv_data)


def main():
    # parse args
    parser = argparse.ArgumentParser(description='plot G-values')
    parser.add_argument('-d', '--db',
                        nargs=1,
                        default=['pdfresults.h5'],
                        help='ProteinDF results file')

    parser.add_argument('-i', '--iteration',
                        nargs=1,
                        default=[0],
                        help='SCF iteration')

    parser.add_argument('--occ',
                        action="store_true",
                        default=False,
                        help='occ only')

    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=None,
                        help='output graph path')
    parser.add_argument('--csv',
                        nargs=1,
                        default=['dos.csv'],
                        help='output csv path')
    parser.add_argument('--yaml',
                        nargs=1,
                        default=None,
                        help='output yaml path')

    parser.add_argument("--frontier",
                        action="store_true",
                        default=False)

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    parser.add_argument('-D', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    # print(args)

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    pdfresults_db = args.db[0]

    iteration = int(args.iteration[0])
    yaml_path = args.yaml[0] if args.yaml else None
    output_path = args.output[0] if args.output else (None if yaml_path else 'g-values.png')
    csv_path = args.csv[0]
    is_frontier = args.frontier

    AU2EV = 27.2116

    # load DB
    pdfparam = pdf.PdfParam_H5()
    pdfparam.open(pdfresults_db)

    # eLevels
    method = pdfparam.method
    run_type = "rks"
    HOMO_level = pdfparam.get_HOMO_level(run_type)

    if iteration == 0:
        iteration = pdfparam.iterations
    print("iteration: {}".format(iteration))
    eigvals = pdfparam.get_energy_level(method, iteration)
    eigvals *= AU2EV

    # make CSC and g-values
    CSC_mat_path = "CSC.mat"
    g_vtr_path = "g-values.vtr"
    args1 = [
        "component",
        "-v",
        "-i", iteration,
        "-S", CSC_mat_path,
        "-g", g_vtr_path
    ]
    # print(args)
    pdf.run_pdf(args1)

    vtr = pdf.Vector()
    vtr.load(g_vtr_path)

    num_of_MOs = len(vtr)
    print("HOMO: {}".format(HOMO_level))

    data_x_occ = [ eigvals[i] for i in range(HOMO_level +1)]
    data_y_occ = [ vtr.get(i) for i in range(HOMO_level +1)]
    data_x_vir = [ eigvals[i] for i in range(HOMO_level +1, len(vtr))]
    data_y_vir = [ vtr.get(i) for i in range(HOMO_level +1, len(vtr))]

    if yaml_path:
        plot_data = make_gvalues_plot_data(
            pdfresults_db,
            method,
            iteration,
            HOMO_level,
            data_x_occ,
            data_y_occ,
            data_x_vir,
            data_y_vir,
        )
        save_gvalues_yaml(plot_data, yaml_path)

    # plot data
    if output_path:
        graph = pdf.DfGvalsGraph()
        if is_frontier:
            graph.xmin = -20.0
            graph.xmax = 5.0

        graph.set_data(data_x_occ, data_y_occ, data_x_vir, data_y_vir)
        print(output_path)
        graph.save(output_path)


if __name__ == '__main__':
    main()

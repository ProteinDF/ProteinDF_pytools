#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
from types import *
import argparse
import csv
import yaml

from matplotlib import use
import numpy

import proteindf_bridge as bridge
import proteindf_tools as pdf

import logging
logger = logging.getLogger(__name__)


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


def make_dos_plot_data(
    h5_path,
    method,
    iteration,
    use_func,
    resolution,
    occ_only,
    output_data,
    num_of_groups,
):
    series = []
    for i in range(num_of_groups):
        points = [
            _FlowMapping({"x": float(row[0]), "y": float(row[i + 1])})
            for row in output_data
        ]
        series.append(
            {
                "name": "y_{}".format(i),
                "points": points,
            }
        )
    return {
        "schema": "pdf-plot/1",
        "kind": "dos",
        "title": "Density of State",
        "source": {
            "tool": "pdf-plot-dos.py",
            "result": h5_path,
            "method": method,
            "iteration": int(iteration),
            "function": use_func,
            "resolution": float(resolution),
            "occ_only": bool(occ_only),
        },
        "axes": {
            "x": {
                "label": "energy level / eV",
                "unit": "eV",
            },
            "y": {
                "label": "DOS",
                "unit": None,
            },
        },
        "series": series,
    }


def save_dos_yaml(plot_data, output_path):
    data = _prepare_plot_data_for_yaml(plot_data)
    content = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)


def dos(func, eps, omega=0.005, eLevels=[], mo_groups=[]):
    answer = [0.0] * (len(mo_groups) + 1)
    for index, level in enumerate(eLevels):
        v = func(eps, omega, level)

        answer[0] += v
        for group_index in range(len(mo_groups)):
            # print("group index: {}".format(group_index))
            if mo_groups[group_index] is not None:
                if index in mo_groups[group_index]:
                    answer[group_index + 1] += v
    return answer


def func_Gaussian(eps, omega, level):
    answer = 0.0
    coef1 = omega * math.sqrt(0.5 * math.pi)
    coef2 = -2.0 / (omega * omega)

    x = eps - level
    x2 = x * x
    answer += math.exp(coef2 * x2)

    answer *= coef1
    return answer


def func_Lorentzian(eps, omega, level):
    answer = 0.0
    omega2 = omega * omega

    x = eps - level
    x2 = x * x
    answer += 1.0 / (omega2 + 4.0 * x2)

    answer *= (2.0 * omega)/(math.pi)
    return answer


def dosGaussian(eps, omega=0.005, eLevels=[]):
    answer = 0.0
    coef1 = omega * math.sqrt(0.5 * math.pi)
    coef2 = -2.0 / (omega * omega)

    for el in eLevels:
        x = eps - el
        x2 = x * x
        answer += math.exp(coef2 * x2)

    answer *= coef1
    return answer


def dosLorentzian(eps, omega=0.005, eLevels=[]):
    answer = 0.0
    omega2 = omega * omega

    for el in eLevels:
        x = eps - el
        x2 = x * x
        answer += 1.0 / (omega2 + 4.0 * x2)

    answer *= (2.0 * omega)/(math.pi)
    return answer


def output_csv(csv_data, path):
    with open(path, 'w', newline='') as f:
        writer = csv.writer(f)
        writer.writerow(csv_data)


def main():
    # parse args
    parser = argparse.ArgumentParser(description='plot density of state')
    parser.add_argument('-d', '--db',
                        nargs=1,
                        default=['pdfresults.h5'],
                        help='ProteinDF results file')

    group_func = parser.add_mutually_exclusive_group()
    group_func.add_argument('--lorentzian',
                            action="store_true",
                            default=False,
                            help='use Lorentzian')
    group_func.add_argument('--gaussian',
                            action="store_true",
                            default=False,
                            help='use Gaussian')

    parser.add_argument('-r', '--resolution',
                        nargs=1,
                        default=[0.1],
                        help='resolution (default: 0.1 eV)')
    parser.add_argument('--occ',
                        action="store_true",
                        default=False,
                        help='occ only')
    parser.add_argument('--mogroup',
                        nargs=1,
                        default=[""],
                        help='msgpack describing mo list for grouping')

    parser.add_argument('-i', '--iteration',
                        nargs=1,
                        default=[0],
                        help='SCF iteration')

    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=None,
                        help='output graph path')
    parser.add_argument('--csv',
                        nargs=1,
                        default=None,
                        help='output csv path')
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
    # print(args)

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    pdfresults_db = args.db[0]

    iteration = args.iteration[0]
    mo_group_path = args.mogroup[0]

    use_func = "gaussian"
    if args.lorentzian:
        use_func = "lorentzian"

    resolution = float(args.resolution[0])
    occ_only = args.occ
    yaml_path = args.yaml[0] if args.yaml else None
    output_path = args.output[0] if args.output else (None if yaml_path else 'dos.png')
    csv_path = args.csv[0] if args.csv else (None if (yaml_path and not output_path) else 'dos.csv')

    AU2EV = 27.2116

    #
    if verbose:
        print("resolution: {}".format(resolution))
        print("func: {}".format(use_func))

    # load DB
    if os.path.exists(pdfresults_db) != True:
        sys.exit("file not found: {}".format(pdfresults_db))
    pdfparam = pdf.PdfParam_H5()
    pdfparam.open(pdfresults_db)

    # mo_group
    mo_groups = []
    if len(mo_group_path) > 0:
        if verbose:
            print("MO group path: {}".format(mo_group_path))
        mo_groups = bridge.load_msgpack(mo_group_path)

    # eLevels
    method = pdfparam.method
    run_type = "rks"
    HOMO_level = pdfparam.get_HOMO_level(run_type)

    if iteration == 0:
        iteration = pdfparam.iterations
    if verbose:
        print("iteration: {}".format(iteration))
    eigvals = pdfparam.get_energy_level(method, iteration)

    # make data
    e_levels = []
    if occ_only:
        e_levels = [None] * HOMO_level
        for i in range(HOMO_level):
            e_levels[i] = eigvals[i] * AU2EV
    else:
        e_levels = [i * AU2EV for i in eigvals]

    # calc intensity
    output_data = []
    minLevel = -20.0
    maxLevel = 5.0
    steps = int((maxLevel - minLevel) / resolution)
    func = func_Lorentzian if use_func == "lorentzian" else func_Gaussian
    for step in range(steps + 1):
        row = []
        e = minLevel + resolution * step
        intensities = dos(func, e, resolution, e_levels, mo_groups)

        row.append(e)
        row.extend(intensities)
        output_data.append(row)

    num_of_groups = len(mo_groups) + 1  # `+1` means `all`
    if yaml_path:
        plot_data = make_dos_plot_data(
            pdfresults_db,
            method,
            iteration,
            use_func,
            resolution,
            occ_only,
            output_data,
            num_of_groups,
        )
        save_dos_yaml(plot_data, yaml_path)

    # output data as csv
    if csv_path:
        with open(csv_path, 'w') as f:
            writer = csv.writer(f)
            writer.writerows(output_data)

    # plot data
    if output_path:
        graph = pdf.DfDosGraph()
        graph.load_data(csv_path)

        print("# groups: ", num_of_groups)
        for i in range(num_of_groups):
            graph.add_draw_series(i)

        graph.title = "Density of State"
        graph.xmax = maxLevel
        graph.xmin = minLevel
        graph.xlabel = "energy level / eV"
        graph.ylabel = "DOS"
        if verbose:
            print("output graph: {}".format(output_path))
        graph.save(output_path)


if __name__ == '__main__':
    main()

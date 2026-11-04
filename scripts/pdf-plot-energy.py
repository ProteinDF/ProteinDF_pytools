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

import argparse
import logging
import os
import sys
import yaml

import proteindf_tools as pdf

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


def make_total_energy_plot_data(source_path, points):
    series_points = [
        _FlowMapping({"iteration": int(itr), "energy": float(ene)})
        for itr, ene in points
    ]
    return {
        "schema": "pdf-plot/1",
        "kind": "total_energy",
        "title": "total energy history",
        "source": source_path,
        "axes": {
            "x": {
                "label": "iteration",
                "unit": None,
            },
            "y": {
                "label": "energy / a.u.",
                "unit": "a.u.",
            },
        },
        "series": [
            {
                "name": "total_energy",
                "points": series_points,
            }
        ],
    }


def save_total_energy_yaml(plot_data, output_path):
    data = _prepare_plot_data_for_yaml(plot_data)
    content = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)


def plot_total_energy(points, output_path):
    graph = pdf.DfTotalEnergyHistGraph()
    graph._x = [float(p[0]) for p in points]
    graph._y = [float(p[1]) for p in points]
    graph.save(output_path)


def read_total_energies(pdfparam):
    iterations = pdfparam.iterations
    points = []
    for itr in range(1, iterations + 1):
        try:
            val = pdfparam.get_total_energy(itr)
        except Exception as e:
            logger.warning("iteration %d: failed to read total energy: %s, skipping", itr, e)
            continue
        if val is None:
            logger.warning("iteration %d: total energy is not available, skipping", itr)
            continue
        points.append((itr, float(val)))
    return points


def main():
    parser = argparse.ArgumentParser(description="plot total energy history")
    parser.add_argument(
        "pdfresults_db",
        nargs="?",
        default="pdfresults.h5",
        help="ProteinDF results file (default: pdfresults.h5)",
    )
    parser.add_argument(
        "-o",
        "--output",
        nargs=1,
        default=None,
        help="output graph path (default: energy.png)",
    )
    parser.add_argument(
        "--yaml",
        nargs=1,
        default=None,
        help="output yaml path",
    )
    parser.add_argument(
        "-v",
        "--verbose",
        action="store_true",
        default=False,
    )
    parser.add_argument(
        "-D",
        "--debug",
        action="store_true",
        default=False,
    )
    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    elif args.verbose:
        logging.basicConfig(level=logging.INFO)
    else:
        logging.basicConfig(level=logging.WARNING)

    db_path = args.pdfresults_db
    yaml_path = args.yaml[0] if args.yaml else None
    output_path = args.output[0] if args.output else (None if yaml_path else "energy.png")

    if args.verbose:
        print("source: %s" % db_path)
        if output_path:
            print("output: %s" % output_path)
        if yaml_path:
            print("yaml: %s" % yaml_path)

    pdfparam = pdf.PdfParam_H5()
    pdfparam.open(db_path)

    points = read_total_energies(pdfparam)

    if yaml_path:
        plot_data = make_total_energy_plot_data(db_path, points)
        save_total_energy_yaml(plot_data, yaml_path)

    if output_path:
        plot_total_energy(points, output_path)


if __name__ == "__main__":
    main()

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
import tempfile
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


def make_elevel_history_plot_data(source_path, run_type, homo_level, points):
    series_points = [
        _FlowMapping({"iteration": int(itr), "level": int(lvl), "energy_ev": float(ene)})
        for itr, lvl, ene in points
    ]
    return {
        "schema": "pdf-plot/1",
        "kind": "elevel_history",
        "title": "energy level history",
        "source": source_path,
        "run_type": run_type,
        "HOMO_level": int(homo_level) if homo_level is not None else None,
        "axes": {
            "x": {
                "label": "iteration",
                "unit": None,
            },
            "y": {
                "label": "energy / eV",
                "unit": "eV",
            },
        },
        "series": [
            {
                "name": "eigvals_hist.%s" % run_type,
                "points": series_points,
            }
        ],
    }


def save_elevel_history_yaml(plot_data, output_path):
    data = _prepare_plot_data_for_yaml(plot_data)
    content = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)


def plot_energy_level_history(points, homo_level, output_path):
    graph = pdf.DfEnergyLevelHistoryGraphH()
    if homo_level is not None:
        graph.set_HOMO_level(homo_level)

    with tempfile.TemporaryDirectory() as tmpdir:
        tmp_dat = os.path.join(tmpdir, "eigvals_hist.dat")
        with open(tmp_dat, "w", encoding="utf-8") as f:
            for itr, level, e in points:
                f.write("%d, %d, % 16.10f\n" % (itr, level, e))
        graph.load_data(tmp_dat)

    graph.save(output_path)


def read_energy_levels(pdfparam, run_type):
    """
    Read energy levels for run_type across iterations 1 to pdfparam.iterations.

    Returns:
        points: list of (iteration, level, energy_ev)
    """
    iterations = pdfparam.iterations
    points = []
    for itr in range(1, iterations + 1):
        try:
            eigvals = pdfparam.get_energy_level(run_type, itr)
        except Exception as e:
            logger.warning(
                "run_type %s, iteration %d: failed to read energy levels: %s, skipping",
                run_type,
                itr,
                e,
            )
            continue
        if eigvals is None or len(eigvals) == 0:
            logger.warning(
                "run_type %s, iteration %d: energy levels are not available, skipping",
                run_type,
                itr,
            )
            continue
        for level, e in enumerate(eigvals):
            energy_ev = float(e) * 27.2116
            points.append((itr, level, energy_ev))
    return points


def _get_opt_val(val):
    if isinstance(val, list):
        return val[0] if len(val) > 0 else None
    return val


def main():
    parser = argparse.ArgumentParser(description="plot energy level history")
    parser.add_argument(
        "pdfresults_db",
        nargs="?",
        default="pdfresults.h5",
        help="ProteinDF results file (default: pdfresults.h5)",
    )
    parser.add_argument(
        "-r",
        "--run_type",
        nargs=1,
        default=None,
        help="orbital run_type (default: all orbital run_types)",
    )
    parser.add_argument(
        "-o",
        "--output",
        nargs=1,
        default=None,
        help="output graph path",
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
    run_type_arg = _get_opt_val(args.run_type)
    output_arg = _get_opt_val(args.output)
    yaml_arg = _get_opt_val(args.yaml)

    pdfparam = pdf.PdfParam_H5()
    pdfparam.open(db_path)

    if run_type_arg:
        run_types = [run_type_arg]
    else:
        run_types = pdfparam.orbital_run_types()

    if not run_types:
        logger.warning("no orbital run_types found")
        return

    is_multi = len(run_types) > 1

    for rt in run_types:
        if yaml_arg:
            if is_multi:
                stem, ext = os.path.splitext(yaml_arg)
                yaml_path = "%s.%s%s" % (stem, rt, ext)
            else:
                yaml_path = yaml_arg
        else:
            yaml_path = None

        if output_arg:
            if is_multi:
                stem, ext = os.path.splitext(output_arg)
                png_path = "%s.%s%s" % (stem, rt, ext)
            else:
                png_path = output_arg
        elif yaml_arg is not None:
            png_path = None
        else:
            png_path = "eigvals_hist.%s.png" % rt

        if args.verbose:
            print("run_type: %s" % rt)
            if png_path:
                print("output: %s" % png_path)
            if yaml_path:
                print("yaml: %s" % yaml_path)

        homo_level = pdfparam.get_HOMO_level(rt)
        points = read_energy_levels(pdfparam, rt)

        if yaml_path:
            plot_data = make_elevel_history_plot_data(db_path, rt, homo_level, points)
            save_elevel_history_yaml(plot_data, yaml_path)

        if png_path:
            plot_energy_level_history(points, homo_level, png_path)


if __name__ == "__main__":
    main()

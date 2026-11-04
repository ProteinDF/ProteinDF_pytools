#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
from types import *
import argparse
import yaml

import proteindf_bridge as bridge
import proteindf_tools as pdf

import logging


class _FlowMapping(dict):
    """Mapping serialized in YAML inline flow style."""

    pass


class _FlowSequence(list):
    """Sequence serialized in YAML inline flow style."""

    pass


yaml.add_representer(
    _FlowMapping,
    lambda dumper, data: dumper.represent_mapping(
        "tag:yaml.org,2002:map", data, flow_style=True
    ),
    Dumper=yaml.SafeDumper,
)

yaml.add_representer(
    _FlowSequence,
    lambda dumper, data: dumper.represent_sequence(
        "tag:yaml.org,2002:seq", data, flow_style=True
    ),
    Dumper=yaml.SafeDumper,
)


def _prepare_plot_data_for_yaml(plot_data):
    if not isinstance(plot_data, dict) or "series" not in plot_data:
        return plot_data
    series_list = plot_data.get("series")
    if not isinstance(series_list, list):
        return plot_data
    new_data = dict(plot_data)
    new_series = []
    for s in series_list:
        if not isinstance(s, dict):
            new_series.append(s)
            continue
        new_s = dict(s)
        if "points" in new_s and isinstance(new_s["points"], list):
            new_s["points"] = [
                p if isinstance(p, _FlowMapping) else _FlowMapping(p)
                for p in new_s["points"]
            ]
        if "values" in new_s and isinstance(new_s["values"], list):
            new_s["values"] = [
                (
                    row
                    if isinstance(row, _FlowSequence)
                    else (_FlowSequence(row) if isinstance(row, list) else row)
                )
                for row in new_s["values"]
            ]
        new_series.append(new_s)
    new_data["series"] = new_series
    return new_data


def make_mat_plot_data(mat_path, mat, title="matrix value"):
    rows = int(mat.rows)
    cols = int(mat.cols)
    values = [
        _FlowSequence([float(mat.get(i, j)) for j in range(cols)])
        for i in range(rows)
    ]
    return {
        "schema": "pdf-plot/1",
        "kind": "mat",
        "title": title,
        "source": {
            "tool": "pdf-plot-mat.py",
            "matrix": mat_path,
            "rows": rows,
            "cols": cols,
        },
        "axes": {
            "x": {
                "label": "column",
                "unit": None,
            },
            "y": {
                "label": "row",
                "unit": None,
            },
        },
        "series": [
            {
                "name": "matrix",
                "rows": rows,
                "cols": cols,
                "values": values,
            }
        ],
    }


def save_mat_yaml(plot_data, output_path):
    data = _prepare_plot_data_for_yaml(plot_data)
    content = yaml.safe_dump(
        data, sort_keys=False, allow_unicode=True, width=float("inf")
    )
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)


def main():
    # parse args
    parser = argparse.ArgumentParser(description="plot matrix")
    parser.add_argument("matrix", nargs=1, help="matrix path")
    parser.add_argument("-o", "--output", nargs=1, default=None, help="output graph path")
    parser.add_argument("--yaml", nargs=1, default=None, help="output yaml path")
    parser.add_argument("--width", type=float, default=20, help="graph width (inch)")
    parser.add_argument("--height", type=float, default=20, help="graph height (inch)")
    parser.add_argument("--vmax", type=float)
    parser.add_argument("--vmin", type=float)
    parser.add_argument("--title", nargs=1, default=["matrix value"], help="graph title")
    parser.add_argument("--no-x-tick-labels", action="store_true", default=False, help="draw x-tick labels")
    parser.add_argument("--no-y-tick-labels", action="store_true", default=False, help="draw y-tick labels")
    parser.add_argument("--cmap", nargs=1, default=["bwr"], help="color map name")
    parser.add_argument("--write-values", action="store_true", default=False)
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("-D", "--debug", action="store_true", default=False)
    args = parser.parse_args()

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    mat_path = args.matrix[0]
    if verbose:
        print("matrix path: {}".format(mat_path))

    yaml_path = args.yaml[0] if args.yaml else None
    output_path = args.output[0] if args.output else (None if yaml_path else "mat.png")
    if verbose:
        if output_path:
            print("output path: {}".format(output_path))
        if yaml_path:
            print("yaml path: {}".format(yaml_path))

    width = args.width
    height = args.height
    if verbose:
        print("width: {}".format(width))
        print("height: {}".format(height))

    vmax = None
    if args.vmax != None:
        vmax = args.vmax
        if verbose:
            print("vmax: {}".format(vmax))
    vmin = None
    if args.vmin != None:
        vmin = args.vmin
        if verbose:
            print("vmin: {}".format(vmin))

    title = args.title[0]
    x_tick_labels = not args.no_x_tick_labels
    y_tick_labels = not args.no_y_tick_labels
    if verbose:
        print("draw x tick labels: {}".format(x_tick_labels))
        print("draw y tick labels: {}".format(y_tick_labels))

    cmap = args.cmap[0]
    if verbose:
        print("cmap: {}".format(cmap))

    write_values = args.write_values
    if verbose:
        print("write values: {}".format(write_values))

    # setup matrix
    mat = get_matrix(mat_path, verbose)
    rows = mat.rows
    cols = mat.cols

    if yaml_path:
        plot_data = make_mat_plot_data(mat_path, mat, title=title)
        save_mat_yaml(plot_data, yaml_path)

    # matrix elements
    if output_path:
        mat_plot = pdf.DfMatrixGraph(mat, figsize=(width, height))
        if vmax:
            mat_plot.vmax = vmax
        if vmin:
            mat_plot.vmin = vmin

        mat_plot.xticks = []
        mat_plot.xticklabels = []
        if x_tick_labels:
            xticks = [x for x in range(0, rows)]
            xticklabels = ["{}".format(x) for x in range(1, rows + 1)]
            mat_plot.xticks = xticks
            mat_plot.xticklabels = xticklabels

        mat_plot.yticks = []
        mat_plot.yticklabels = []
        if y_tick_labels:
            yticks = [y for y in range(0, cols)]
            yticklabels = ["{}".format(y) for y in range(1, cols + 1)]
            mat_plot.yticks = yticks
            mat_plot.yticklabels = yticklabels
        mat_plot.title = title
        mat_plot.cmap = cmap
        mat_plot.should_write_values = write_values

        mat_plot.save(output_path)


def get_matrix(path, verbose=False):
    if verbose:
        sys.stderr.write("load matrix: {}\n".format(path))
    mat = None

    if pdf.Matrix.is_loadable(path):
        mat = pdf.Matrix()
        mat.load(path)
    elif pdf.SymmetricMatrix.is_loadable(path):
        mat = pdf.SymmetricMatrix()
        mat.load(path)
        mat = mat.get_general_matrix()
    else:
        print("matrix file cannnot load.")
        raise ValueError("cannot load the file: {}".format(path))

    return mat


def distance_matrix(atomgroup, verbose=False):
    """
    create distance matrix

    support flat-atomgroup only
    """
    atomlist = atomgroup.get_atom_list()
    num_of_atoms = len(atomlist)
    if verbose:
        sys.stderr.write("num of atoms: {}n".format(num_of_atoms))
    mat = bridge.SymmetricMatrix(num_of_atoms)

    for i in range(num_of_atoms):
        for j in range(i):
            distance = atomlist[i].xyz.distance_from(atomlist[j].xyz)
            mat.set(i, j, distance)

    return mat


def dist_vs_matvalue(orb_info, distance_mat, mat, path, verbose=False):
    if verbose:
        sys.stderr.write("matrix size: {} x {}\n".format(mat.rows, mat.cols))

    f = open(path, "w")
    for i in range(mat.rows):
        atom_i = orb_info.get_atom_id(i)
        for j in range(i):
            atom_j = orb_info.get_atom_id(j)

            d = distance_mat.get(atom_i, atom_j)
            v = mat.get(i, j)
            data_str = "{}, {}\n".format(d, v)
            f.write(data_str)
    f.close()


if __name__ == "__main__":
    main()

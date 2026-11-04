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


def make_vector_plot_data(vct_path, vct, plot_abs_values=False):
    return {
        "schema": "pdf-plot/1",
        "kind": "vector",
        "title": "",
        "source": {
            "tool": "pdf-plot-vector.py",
            "vector": vct_path,
            "abs": bool(plot_abs_values),
        },
        "axes": {
            "x": {
                "label": "index",
                "unit": None,
            },
            "y": {
                "label": "value",
                "unit": None,
            },
        },
        "series": [
            {
                "name": "vector",
                "values": [float(v) for v in vct],
            }
        ],
    }


def save_vector_yaml(plot_data, output_path):
    data = _prepare_plot_data_for_yaml(plot_data)
    content = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)



def main():
    # parse args
    parser = argparse.ArgumentParser(description="plot vector")
    parser.add_argument("vector", nargs=1, help="vector path")
    parser.add_argument("-o", "--output", nargs=1, default=None, help="output graph path")
    parser.add_argument("--yaml", nargs=1, default=None, help="output yaml path")
    parser.add_argument("--abs", action="store_true", default=False, help="plot absolute values")
    parser.add_argument("--draw-grid", action="store_true", default=False, help="draw grid")
    parser.add_argument("--yscale", type=str)
    parser.add_argument("--ymax", type=float)
    parser.add_argument("--ymin", type=float)
    # parser.add_argument("--title", nargs=1, default=["matrix value"], help="graph title")
    # parser.add_argument("--no-x-tick-labels", action="store_true", default=False, help="draw x-tick labels")
    # parser.add_argument("--no-y-tick-labels", action="store_true", default=False, help="draw y-tick labels")
    # parser.add_argument("--cmap", nargs=1, default=["bwr"], help="color map name")
    # parser.add_argument("--write-values", action="store_true", default=False)
    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("-D", "--debug", action="store_true", default=False)
    args = parser.parse_args()

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    vct_path = args.vector[0]
    if verbose:
        print("vector path: {}".format(vct_path))

    yaml_path = args.yaml[0] if args.yaml else None
    output_path = args.output[0] if args.output else (None if yaml_path else "vct.png")
    if verbose:
        if output_path:
            print("output path: {}".format(output_path))
        if yaml_path:
            print("yaml path: {}".format(yaml_path))

    plot_abs_values = args.abs

    is_draw_grid = args.draw_grid

    yscale = None
    if args.yscale != None:
        yscale = args.yscale
        if verbose:
            print("yscale: {}".format(yscale))
    ymax = None
    if args.ymax != None:
        ymax = args.ymax
        if verbose:
            print("ymax: {}".format(ymax))
    ymin = None
    if args.ymin != None:
        ymin = args.ymin
        if verbose:
            print("ymin: {}".format(ymin))

    # title = args.title[0]
    # x_tick_labels = not args.no_x_tick_labels
    # y_tick_labels = not args.no_y_tick_labels
    # if verbose:
    #     print("draw x tick labels: {}".format(x_tick_labels))
    #     print("draw y tick labels: {}".format(y_tick_labels))

    # cmap = args.cmap[0]
    # if verbose:
    #     print("cmap: {}".format(cmap))

    # write_values = args.write_values
    # if verbose:
    #     print("write values: {}".format(write_values))

    # setup vector
    vct = get_vector(vct_path, verbose)
    if plot_abs_values:
        vct = vct.abs()

    if yaml_path:
        plot_data = make_vector_plot_data(vct_path, vct, plot_abs_values=plot_abs_values)
        save_vector_yaml(plot_data, yaml_path)

    # vector elements
    if output_path:
        vct_plot = pdf.DfVectorGraph(vct, figsize=(20, 20))

        vct_plot.is_draw_grid = is_draw_grid
        vct_plot.yscale = yscale
        if ymax:
            vct_plot.ymax = ymax
        if ymin:
            vct_plot.ymin = ymin

        # mat_plot.xticks = []
        # mat_plot.xticklabels = []
        # if x_tick_labels:
        #     xticks = [x for x in range(0, rows)]
        #     xticklabels = ["{}".format(x) for x in range(1, rows + 1)]
        #     mat_plot.xticks = xticks
        #     mat_plot.xticklabels = xticklabels

        # mat_plot.yticks = []
        # mat_plot.yticklabels = []
        # if y_tick_labels:
        #     yticks = [y for y in range(0, cols)]
        #     yticklabels = ["{}".format(y) for y in range(1, cols + 1)]
        #     mat_plot.yticks = yticks
        #     mat_plot.yticklabels = yticklabels
        # mat_plot.title = title
        # mat_plot.cmap = cmap
        # mat_plot.should_write_values = write_values

        vct_plot.save(output_path)


def get_vector(path, verbose=False):
    if verbose:
        sys.stderr.write("load matrix: {}\n".format(path))
    vct = None

    if pdf.Vector.is_loadable(path):
        vct = pdf.Vector()
        vct.load(path)
    else:
        print("vector file cannnot load.")
        raise ValueError("cannot load the file: {}".format(path))

    return vct


# def distance_matrix(atomgroup, verbose=False):
#     """
#     create distance matrix

#     support flat-atomgroup only
#     """
#     atomlist = atomgroup.get_atom_list()
#     num_of_atoms = len(atomlist)
#     if verbose:
#         sys.stderr.write("num of atoms: {}n".format(num_of_atoms))
#     mat = bridge.SymmetricMatrix(num_of_atoms)

#     for i in range(num_of_atoms):
#         for j in range(i):
#             distance = atomlist[i].xyz.distance_from(atomlist[j].xyz)
#             mat.set(i, j, distance)

#     return mat


# def dist_vs_matvalue(orb_info, distance_mat, mat, path, verbose=False):
#     if verbose:
#         sys.stderr.write("matrix size: {} x {}\n".format(mat.rows, mat.cols))

#     f = open(path, "w")
#     for i in range(mat.rows):
#         atom_i = orb_info.get_atom_id(i)
#         for j in range(i):
#             atom_j = orb_info.get_atom_id(j)

#             d = distance_mat.get(atom_i, atom_j)
#             v = mat.get(i, j)
#             data_str = "{}, {}\n".format(d, v)
#             f.write(data_str)
#     f.close()


if __name__ == "__main__":
    main()

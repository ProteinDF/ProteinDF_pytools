#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
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


def make_decaymat_plot_data(
    orb_info, distance_mat, mat, mat_path, db_path=None, param_path=None
):
    points = []
    for i in range(mat.rows):
        atom_i = orb_info.get_atom_id(i)
        for j in range(i):
            atom_j = orb_info.get_atom_id(j)
            d = float(distance_mat.get(atom_i, atom_j))
            v = float(mat.get(i, j))
            points.append(_FlowMapping({"x": d, "y": math.fabs(v)}))

    source = {
        "tool": "pdf-plot-decaymat.py",
        "matrix": mat_path,
        "rows": int(mat.rows),
        "cols": int(mat.cols),
    }
    if db_path is not None:
        source["result"] = db_path
    if param_path is not None:
        source["param"] = param_path

    return {
        "schema": "pdf-plot/1",
        "kind": "decaymat",
        "title": "",
        "source": source,
        "axes": {
            "x": {
                "label": "distance",
                "unit": None,
            },
            "y": {
                "label": "magnitude",
                "unit": None,
            },
        },
        "series": [
            {
                "name": "decaymat",
                "points": points,
            }
        ],
    }


def save_decaymat_yaml(plot_data, output_path):
    data = _prepare_plot_data_for_yaml(plot_data)
    content = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)


def main():
    # parse args
    parser = argparse.ArgumentParser(description='plot matrix')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "-d",
        "--db",
        nargs="?",
        action="store",
        const="pdfresults.h5",
        help="ProteinDF results file (.h5 は HDF5、それ以外は SQLite(従来))",
    )
    group.add_argument('-p', '--param',
                       nargs='?',
                       action='store',
                       const='pdfparam.mpac',
                       help='ProteinDF parameter file')
    parser.add_argument('matrix',
                        nargs=1,
                        help='matrix path')
    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=None,
                        help='output graph path')
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
    mat_path = args.matrix[0]
    yaml_path = args.yaml[0] if args.yaml else None
    output_path = args.output[0] if args.output else (None if yaml_path else 'mat.png')

    # setup orbinfo
    if verbose:
        sys.stderr.write("prepare orbital information\n")
    orb_info = None
    molecule = None
    if args.db:
        entry = pdf.open_archive(args.db)
        orb_info = pdf.OrbInfo(entry)
        molecule = entry.get_molecule()
    elif args.param:
        pdfparam = pdf.load_pdfparam(args.param)
        orb_info = pdf.OrbInfo(pdfparam)
        molecule = pdfparam.molecule

    # molecule
    if verbose:
        sys.stderr.write("prepare molecule information\n")
    distance_mat = distance_matrix(molecule)

    # setup matrix
    mat = get_matrix(mat_path, verbose)
    if mat == None:
        logging.critical("cannnot load")
        abort()

    if yaml_path:
        plot_data = make_decaymat_plot_data(
            orb_info,
            distance_mat,
            mat,
            mat_path,
            db_path=args.db,
            param_path=args.param,
        )
        save_decaymat_yaml(plot_data, yaml_path)

    # distance v.s. matrix elements
    if output_path:
        dist_vs_matvalue(orb_info, distance_mat, mat, 'mat.dat')
        distgraph = pdf.DfDistanceVsElementGraph()
        distgraph.load('mat.dat')
        distgraph.save(output_path)


def get_matrix(path, verbose=False):
    if verbose:
        sys.stderr.write("load matrix: {}\n".format(path))
    mat = None

    if pdf.SymmetricMatrix.is_loadable(path):
        mat = pdf.SymmetricMatrix()
        mat.load(path)
        mat = mat.get_general_matrix()
    elif pdf.Matrix.is_loadable(path):
        mat = pdf.Matrix()
        mat.load(path)

    return mat


def distance_matrix(atomgroup, verbose=False):
    '''
    create distance matrix

    support flat-atomgroup only
    '''
    atomlist = atomgroup.get_atom_list()
    num_of_atoms = len(atomlist)
    if verbose:
        sys.stderr.write("num of atoms: {}\n".format(num_of_atoms))
    mat = bridge.SymmetricMatrix(num_of_atoms)

    for i in range(num_of_atoms):
        for j in range(i):
            distance = atomlist[i].xyz.distance_from(atomlist[j].xyz)
            mat.set(i, j, distance)

    return mat


def dist_vs_matvalue(orb_info, distance_mat, mat, path, verbose=False):
    if verbose:
        sys.stderr.write("matrix size: {} x {}\n".format(mat.rows, mat.cols))

    f = open(path, 'w')
    for i in range(mat.rows):
        atom_i = orb_info.get_atom_id(i)
        for j in range(i):
            atom_j = orb_info.get_atom_id(j)

            d = distance_mat.get(atom_i, atom_j)
            v = mat.get(i, j)
            data_str = "{}, {}\n".format(d, v)
            f.write(data_str)
    f.close()


if __name__ == '__main__':
    main()

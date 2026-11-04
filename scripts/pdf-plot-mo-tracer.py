#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import logging
import os
import sys
import yaml

import proteindf_tools as pdf

AU2EV = 27.2116


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


def make_mo_tracer_plot_data(
    path1,
    path2,
    db_path1,
    db_path2,
    csc_path,
    runtype,
    itr1,
    itr2,
    pickup_ratio,
    elevel1_list,
    elevel2_list,
    pairs,
):
    return {
        "schema": "pdf-plot/1",
        "kind": "mo-tracer",
        "title": "",
        "source": {
            "tool": "pdf-plot-mo-tracer.py",
            "path1": path1,
            "path2": path2,
            "result1": db_path1,
            "result2": db_path2,
            "csc_path": csc_path,
            "run_type": runtype,
            "iteration1": int(itr1),
            "iteration2": int(itr2),
            "pickup_ratio": float(pickup_ratio),
        },
        "axes": {
            "x": {
                "label": "",
                "unit": None,
            },
            "y": {
                "label": "energy / eV",
                "unit": "eV",
            },
        },
        "series": [
            {
                "name": "elevel1",
                "values": [float(v) for v in elevel1_list],
            },
            {
                "name": "elevel2",
                "values": [float(v) for v in elevel2_list],
            },
            {
                "name": "pairs",
                "points": [
                    _FlowMapping({"mo1": int(r), "mo2": int(idx), "overlap": float(val)})
                    for r, idx, val in pairs
                ],
            },
        ],
    }


def save_mo_tracer_yaml(plot_data, output_path):
    data = _prepare_plot_data_for_yaml(plot_data)
    content = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)



def find_pdfparam_path(pdf_path, pdfparam_filename="pdfparam.mpac"):
    pdfparam_path = None
    path = os.path.join(pdf_path, pdfparam_filename)
    if os.path.isfile(path):
        pdfparam_path = path

    return pdfparam_path


def find_db(pdf_path, h5_filename="pdfresults.h5"):
    h5_path = None
    path = os.path.join(pdf_path, h5_filename)
    if os.path.isfile(path):
        h5_path = path

    return h5_path


def get_entry(pdfparam_path, h5_path, verbose=False):
    entry = None
    if h5_path:
        if verbose:
            print("load: {}".format(h5_path))
        entry = pdf.PdfParam_H5(h5_path)
    elif pdfparam_path:
        if verbose:
            print("load: {}".format(pdfparam_path))
        entry = pdf.load_pdfparam(pdfparam_path)
    else:
        print("cannot load parameter files")

    return entry


def make_mo_overlap_matrix(path1, path2, pdfparam_path1, pdfparam_path2, CSC_mat_path):
    pdfparam1 = pdf.load_pdfparam(pdfparam_path1)
    cmat_path1 = pdfparam1.get_c_mat_path()

    pdfparam2 = pdf.load_pdfparam(pdfparam_path2)
    cmat_path2 = pdfparam2.get_c_mat_path()

    args = [
        "info-mo-overlap",
        "-v",
        "-o",
        CSC_mat_path,
        pdfparam_path1,
        pdfparam_path2,
        os.path.join(path1, cmat_path1),
        os.path.join(path2, cmat_path2),
    ]
    # print(args)

    pdf.run_pdf(args)


def make_pair(CSC_mat, pickup_ratio=0.1, verbose=False):
    pickup_ratio = float(pickup_ratio)
    print("pickup ratio: {}".format(pickup_ratio))
    pairs = []

    if verbose:
        print("### (data1 MO, data2 MO): overlap ratio ###")
    rows = CSC_mat.rows
    cols = CSC_mat.cols
    for r in range(rows):
        col_vec = CSC_mat.get_row_vector(r)
        norm_col_vec = pdf.Vector(cols)
        norm = 0.0
        for c in range(cols):
            v = abs(col_vec[c])
            norm_col_vec[c] = v
            norm += v
        norm_col_vec *= 1.0 / norm
        col_index = norm_col_vec.argsort().flip()
        assert len(col_index) == cols

        # pickup by ratio
        for rank in range(cols):
            index = col_index[rank]
            value = norm_col_vec[index]

            if value > pickup_ratio:
                pair = [r, index, value]
                pairs.append(pair)
                if verbose:
                    print("{:6d}, {:6d}, {:8.5f}".format(r, index, value))
            else:
                break

        # rank = 0
        # ratio = 0.0
        # while ratio < 0.8:
        #     index = col_index[rank]
        #     # value = abs(norm_col_vec[index])
        #     value = norm_col_vec[index]
        #     if value > pickup_ratio:
        #         pair = [r, index, value]
        #         pairs.append(pair)
        #         if verbose:
        #             print("{}, {}: {:8.5f}".format(r, index, value))
        #             # print("#{}: {} th MO overlap {:5.2f}%".format(rank, index, value * 100.0))

        #     rank += 1
        #     ratio += value

    return pairs


def main():
    # parse args
    parser = argparse.ArgumentParser(description="plot energy levels")

    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("-D", "--debug", action="store_true", default=False)

    parser.add_argument("-m", "--csc-path", default="CSC.mat", type=str)

    parser.add_argument("--pickup-ratio", default=0.1, type=float)

    parser.add_argument("-o", "--output", default=None, type=str)
    parser.add_argument("--yaml", default=None, type=str, help="output yaml path")

    parser.add_argument("--ymin", default=-20, type=float)  # for graph setting (ymin)
    parser.add_argument("--ymax", default=10, type=float)  # for graph setting (ymax)

    parser.add_argument("pdf_path1", help="ProteinDF path1")
    parser.add_argument("pdf_path2", help="ProteinDF path2")
    args = parser.parse_args()
    # print(args)

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    pickup_ratio = args.pickup_ratio
    yaml_path = args.yaml
    output_path = args.output if args.output else (None if yaml_path else "mo-tracer.png")

    path1 = args.pdf_path1
    path2 = args.pdf_path2
    if verbose:
        print("ProteinDF path1: {}".format(path1))
        print("ProteinDF path2: {}".format(path2))

    ymin = args.ymin
    ymax = args.ymax
    if verbose:
        print("ymin: {}".format(ymin))
        print("ymax: {}".format(ymax))

    pdfparam_path1 = find_pdfparam_path(path1)
    pdfparam_path2 = find_pdfparam_path(path2)

    # entry1 = get_entry(pdfparam_path1, None, verbose)
    # entry2 = get_entry(pdfparam_path2, None, verbose)

    # itr1 = entry1.iterations
    # itr2 = entry2.iterations
    # print(itr1, itr2)

    csc_path = args.csc_path
    if verbose:
        print(">>>> make CSC matrix")
        print("CSC matrix path: {}".format(csc_path))
    make_mo_overlap_matrix(path1, path2, pdfparam_path1, pdfparam_path2, csc_path)

    db_path1 = find_db(path1)
    db_path2 = find_db(path2)

    # load DBs
    if verbose:
        print(">>>> read ProteinDF results")
    db1 = pdf.PdfParam_H5(db_path1)
    db2 = pdf.PdfParam_H5(db_path2)

    runtype = "rks"
    itr1 = db1.iterations
    if verbose:
        print("ProteinDF data1 iteration: {}".format(itr1))
    elevel1 = db1.get_energy_level(runtype, itr1)
    elevel1 *= AU2EV
    if verbose:
        print("ProteinDF data1 size of eigvals: {}".format(len(elevel1)))

    itr2 = db2.iterations
    if verbose:
        print("ProteinDF data2 iteration: {}".format(itr2))
    elevel2 = db2.get_energy_level(runtype, itr2)
    elevel2 *= AU2EV
    if verbose:
        print("ProteinDF data2 size of eigvals: {}".format(len(elevel2)))

    # MO overlap
    if verbose:
        print(">>>> read CSC matrix")
    CSC_mat = pdf.Matrix()
    CSC_mat.load(csc_path)
    if verbose:
        print("CSC matrix size: {} x {}".format(CSC_mat.rows, CSC_mat.cols))

    # make pair
    if verbose:
        print(">>>> check MO pairs")
    pairs = make_pair(CSC_mat, pickup_ratio, verbose)

    elevel1_list = elevel1.to_list()
    elevel2_list = elevel2.to_list()

    if yaml_path:
        if verbose:
            print("yaml: {}".format(yaml_path))
        plot_data = make_mo_tracer_plot_data(
            path1,
            path2,
            db_path1,
            db_path2,
            csc_path,
            runtype,
            itr1,
            itr2,
            pickup_ratio,
            elevel1_list,
            elevel2_list,
            pairs,
        )
        save_mo_tracer_yaml(plot_data, yaml_path)

    # graph
    if output_path:
        if verbose:
            print(">>>> make graph")
        graph = pdf.DfEnergyLevelTraceGraph()
        graph.ymin = ymin
        graph.ymax = ymax
        graph.xticks = []
        graph.set_data(elevel1_list, elevel2_list, pairs)

        if verbose:
            print("output: {}".format(output_path))
        graph.save(output_path)


if __name__ == "__main__":
    main()

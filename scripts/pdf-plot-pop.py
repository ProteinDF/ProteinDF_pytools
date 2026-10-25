#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import json
import tempfile
import warnings

import proteindf_bridge as bridge
import proteindf_tools as pdf

import logging


def make_pop_data(pop_vtr, output_path):
    num_of_items = len(pop_vtr)
    with open(output_path, "w") as f:
        for index in range(num_of_items):
            charge = pop_vtr[index]
            f.write("{}, {}\n".format(index + 1, charge))


def make_pop_data_diff(pop_vtr1, pop_vtr2, output_path):
    num_of_items = len(pop_vtr1)
    assert(num_of_items == len(pop_vtr2))

    with open(output_path, "w") as f:
        for index in range(num_of_items):
            charge = pop_vtr1[index] - pop_vtr2[index]
            f.write("{}, {}\n".format(index + 1, charge))


def get_atom_elements(pdfparam):
    if hasattr(pdfparam, "_h5root") and pdfparam._h5root is not None:
        if "molecule/xyzc" in pdfparam._h5root:
            ds = pdfparam._h5root["molecule/xyzc"]
            elements = []
            for row in ds:
                sym = row[0].decode("utf-8") if isinstance(row[0], bytes) else str(row[0])
                elements.append(sym)
            if len(elements) > 0:
                return elements
    if hasattr(pdfparam, "molecule") and pdfparam.molecule is not None:
        atoms = pdfparam.molecule.get_atom_list()
        if len(atoms) > 0:
            warnings.warn(
                "molecule/xyzc not found; falling back to pdfparam.molecule.get_atom_list()",
                UserWarning,
            )
            return [atom.symbol for atom in atoms]
    raise ValueError("Molecule geometry / atom list not found in results file")


def make_pop_diff_plot_data(
    h5_path,
    last_iteration,
    elements,
    pop_vtr_last,
    pop_vtr_ref,
    run_type="rks",
    reference_iteration=1,
):
    num_of_items = len(pop_vtr_last)
    if len(pop_vtr_ref) != num_of_items:
        raise ValueError(
            "Population vectors length mismatch: {} vs {}".format(
                len(pop_vtr_last), len(pop_vtr_ref)
            )
        )
    if len(elements) != num_of_items:
        raise ValueError(
            "Number of elements ({}) does not match population vector length ({})".format(
                len(elements), num_of_items
            )
        )

    atoms = []
    points = []
    for index in range(num_of_items):
        diff = float(pop_vtr_last[index] - pop_vtr_ref[index])
        atom_idx = index + 1
        atoms.append(
            {
                "index": atom_idx,
                "element": elements[index],
                "label": "",
                "chain": None,
                "residue_number": None,
                "residue_name": None,
                "atom_name": None,
            }
        )
        points.append(
            {
                "x": atom_idx,
                "y": diff,
                "atom": atom_idx,
            }
        )

    return {
        "schema": "pdf-plot/1",
        "kind": "scatter",
        "title": "Mulliken charge difference",
        "source": {
            "tool": "pdf-plot-pop.py",
            "result": h5_path,
            "run_type": run_type,
            "iteration": last_iteration,
            "reference_iteration": reference_iteration,
        },
        "axes": {
            "x": {
                "label": "atom index",
                "unit": None,
            },
            "y": {
                "label": "Mulliken charge difference (iteration {} \u2212 iteration {})".format(
                    last_iteration, reference_iteration
                ),
                "unit": "e",
            },
        },
        "atoms": atoms,
        "series": [
            {
                "name": "Mulliken (diff)",
                "points": points,
            }
        ],
    }


def save_pop_diff_json(plot_data, output_path):
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(plot_data, f, indent=2, ensure_ascii=False)
        f.write("\n")


def main():
    # parse args
    parser = argparse.ArgumentParser(description='plot population')
    group = parser.add_mutually_exclusive_group(required=True)
    # group.add_argument('-e', '--elevel_vector',
    #                   nargs=1,
    #                   help='ProteinDF energy level vector file')
    group.add_argument('--h5',
                       nargs='?',
                       action='store',
                       const='pdfresults.h5',
                       help='ProteinDF results file (HDF5 format)')

    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=['pop.png'],
                        help='output graph path')
    parser.add_argument('--json',
                        nargs=1,
                        default=None,
                        help='output json path')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    # setting
    verbose = args.verbose
    pdfparam_h5_path = args.h5
    output_path = args.output[0]
    json_path = args.json[0] if args.json else None

    #
    if verbose:
        print("H5 path: {}".format(pdfparam_h5_path))
        print("output path: {}".format(output_path))
        if json_path:
            print("json path: {}".format(json_path))

    #
    pdfparam = pdf.PdfParam_H5()
    pdfparam.open(pdfparam_h5_path)

    run_type = "rks"
    itr = pdfparam.iterations

    pop_vtr1 = pdfparam.get_pop_mulliken_atom(run_type, itr)
    if pop_vtr1 is None or len(pop_vtr1) == 0:
        raise ValueError(
            "Population data (Mulliken) for iteration {} not found in {}".format(
                itr, pdfparam_h5_path
            )
        )

    pop_vtr2 = pdfparam.get_pop_mulliken_atom(run_type, 1)
    if pop_vtr2 is None or len(pop_vtr2) == 0:
        raise ValueError(
            "Population data (Mulliken) for iteration 1 not found in {}".format(
                pdfparam_h5_path
            )
        )

    if json_path:
        elements = get_atom_elements(pdfparam)
        plot_data = make_pop_diff_plot_data(
            pdfparam_h5_path,
            itr,
            elements,
            pop_vtr1,
            pop_vtr2,
            run_type=run_type,
            reference_iteration=1,
        )
        save_pop_diff_json(plot_data, json_path)

        with tempfile.NamedTemporaryFile(mode="w", suffix=".txt", delete=False) as tmp:
            tmp_path = tmp.name
        try:
            make_pop_data_diff(pop_vtr1, pop_vtr2, tmp_path)
            pop_graph = pdf.DfPopulationGraph()
            pop_graph.load_data(tmp_path)
            pop_graph.save(output_path)
        finally:
            if os.path.exists(tmp_path):
                os.remove(tmp_path)
    else:
        data_path = "pop.mulliken.txt"
        # make_pop_data(pop_vtr, data_path)
        make_pop_data_diff(pop_vtr1, pop_vtr2, data_path)

        # draw
        pop_graph = pdf.DfPopulationGraph()
        pop_graph.load_data(data_path)
        pop_graph.save(output_path)


if __name__ == '__main__':
    main()

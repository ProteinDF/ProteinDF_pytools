#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import argparse
import json
import subprocess
import tempfile
import warnings

import h5py
import proteindf_bridge as bridge
import proteindf_tools as pdf
from proteindf_tools.poputils import compute_mulliken_pop, parse_mulliken_pop

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


def open_pdfparam_h5(filepath):
    """Open HDF5 file with PdfParam_H5, falling back to read-only mode if not writable.

    Parameters
    ----------
    filepath : str
        Path to the HDF5 file.

    Returns
    -------
    tuple of (PdfParam_H5, bool)
        The opened PdfParam_H5 instance and a boolean flag indicating whether the file is writable.
    """
    pdfparam = pdf.PdfParam_H5()
    if os.path.exists(filepath) and not os.access(filepath, os.W_OK):
        pdfparam._state["h5filepath"] = filepath
        pdfparam._h5file = h5py.File(filepath, "r")
        pdfparam._h5root = pdfparam._h5file["/"]
        pdfparam.load_model()
        TEs = pdfparam.get_TEs()
        if TEs is not None:
            for itr in range(pdfparam.iterations):
                pdfparam.set_total_energy(itr + 1, TEs[itr])
        return pdfparam, False

    try:
        pdfparam.open(filepath)
        return pdfparam, True
    except (PermissionError, OSError):
        pdfparam._state["h5filepath"] = filepath
        pdfparam._h5file = h5py.File(filepath, "r")
        pdfparam._h5root = pdfparam._h5file["/"]
        pdfparam.load_model()
        TEs = pdfparam.get_TEs()
        if TEs is not None:
            for itr in range(pdfparam.iterations):
                pdfparam.set_total_energy(itr + 1, TEs[itr])
        return pdfparam, False


def save_pop_to_h5(pdfparam, run_type, iteration, pop, is_writable):
    """Save population to HDF5 if writable, otherwise emit a warning.

    Parameters
    ----------
    pdfparam : PdfParam_H5
        ProteinDF parameter object.
    run_type : str
        Calculation run type (e.g. "rks").
    iteration : int
        Iteration number.
    pop : pdf.Vector
        Population vector.
    is_writable : bool
        Whether the file is writable.
    """
    if not is_writable:
        warnings.warn(
            "Cannot save population to {}: file is read-only".format(pdfparam.filepath),
            UserWarning,
        )
        return
    try:
        pdfparam.save_pop_mulliken_atom(run_type, iteration, pop)
        if hasattr(pdfparam, "_h5file") and pdfparam._h5file is not None:
            pdfparam._h5file.flush()
    except (PermissionError, OSError) as e:
        warnings.warn(
            "Failed to save population to {}: {}".format(pdfparam.filepath, e),
            UserWarning,
        )


def get_density_matrix_safe(pdfparam, run_type, iteration):
    try:
        return pdfparam.get_density_matrix(run_type, iteration)
    except Exception:
        return None


def get_input_density_matrix_safe(pdfparam, run_type, iteration):
    try:
        return pdfparam.get_input_density_matrix(run_type, iteration)
    except Exception:
        return None


def get_s_matrix_safe(pdfparam):
    try:
        return pdfparam.get_s_matrix()
    except Exception:
        return None


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
    pdfparam, is_writable = open_pdfparam_h5(pdfparam_h5_path)

    run_type = "rks"
    itr = pdfparam.iterations

    need_save_last = False
    need_save_1 = False

    pop_vtr1 = pdfparam.get_pop_mulliken_atom(run_type, itr)
    if pop_vtr1 is None or len(pop_vtr1) == 0:
        P_last = get_density_matrix_safe(pdfparam, run_type, itr)
        S = get_s_matrix_safe(pdfparam)
        if P_last is None or S is None:
            raise ValueError(
                "Population data (Mulliken) for iteration {} not found in {}".format(
                    itr, pdfparam_h5_path
                )
            )
        pop_vtr1 = compute_mulliken_pop(pdfparam, P_last, S)
        need_save_last = True

    pop_vtr2 = pdfparam.get_pop_mulliken_atom(run_type, 1)
    if pop_vtr2 is None or len(pop_vtr2) == 0:
        P_in = get_input_density_matrix_safe(pdfparam, run_type, 1)
        S = get_s_matrix_safe(pdfparam)
        if P_in is None or S is None:
            raise ValueError(
                "Population data (Mulliken) for iteration 1 not found in {}".format(
                    pdfparam_h5_path
                )
            )
        pop_vtr2 = compute_mulliken_pop(pdfparam, P_in, S)
        need_save_1 = True

    if need_save_last:
        save_pop_to_h5(pdfparam, run_type, itr, pop_vtr1, is_writable)
    if need_save_1:
        save_pop_to_h5(pdfparam, run_type, 1, pop_vtr2, is_writable)

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

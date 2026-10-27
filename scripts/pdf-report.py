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

import os
import sys
import argparse

import proteindf_bridge as bridge
import proteindf_tools as pdf

import logging
import logging.config

logger = logging.getLogger(__name__)


def main():
    # parse args
    parser = argparse.ArgumentParser(description="make ProteiDF report")
    parser.add_argument(
        "pdfresults_db",
        nargs="?",
        default="pdfresults.db",
        help="ProteinDF results file",
    )
    parser.add_argument(
        "output_dir", nargs="?", default="report", help="output directory"
    )

    parser.add_argument("-v", "--verbose", action="store_true", default=False)
    parser.add_argument("-d", "--debug", action="store_true", default=False)
    parser.add_argument(
        "-p", "--profile", nargs="?", const="pdf-report.stat", help="do profiling"
    )

    args = parser.parse_args()

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    do_profile = args.profile

    output_dir = args.output_dir
    entry = pdf.open_archive(args.pdfresults_db)

    # run
    if do_profile != None:
        import cProfile

        cProfile.runctx("do_report(output_dir, entry)", globals(), locals(), do_profile)
    else:
        do_report(output_dir, entry)


def do_report(output_dir, entry):
    # make output dir
    if not os.path.isdir(output_dir):
        os.mkdir(output_dir)

    # ------------------
    iterations = entry.iterations
    logging.debug("iterations=%d" % (iterations))

    plot_convergence_check(entry, output_dir)
    plot_convergence_energy_level(entry, output_dir)

    # plot_energy_level(pdfdata,
    #                  workdir + '/level.png')

    # plot_elapsed_time(entry, output_dir)

    if entry.scf_converged:
        contents = get_rst(entry)
        print(contents)


def plot_convergence_check(entry, output_dir):
    """
    TEのヒストリ
    """
    assert isinstance(entry, (pdf.PdfArchive, pdf.PdfParam_H5))
    itr = entry.iterations
    iterations = range(1, itr + 1)

    # TE
    TEs = [entry.get_total_energy(i) for i in iterations]
    data_path = os.path.join(output_dir, "TE_hist.dat")
    graph_path = os.path.join(output_dir, "TE_hist.png")
    dat = open(data_path, "w")
    for itr, TE in enumerate(TEs):
        if TE:
            dat.write("%d, % 16.10f\n" % (itr + 1, TE))
    dat.close()

    graph = pdf.DfTotalEnergyHistGraph()
    graph.load_data(data_path)
    graph.save(graph_path)


def plot_convergence_energy_level(entry, output_dir):
    """
    EnergyLevelのヒストリ
    """
    itr = entry.iterations
    method = entry.method
    HOMO_level = entry.get_HOMO_level(method)

    data_path = os.path.join(output_dir, "eigvals_hist.dat")
    dat = open(data_path, "w")
    for itr in range(1, entry.iterations + 1):
        eigvals = entry.get_energylevel(method, itr)
        if eigvals:
            for level, e in enumerate(eigvals):
                e *= 27.2116
                dat.write("%d, %d, % 16.10f\n" % (itr, level, e))
    dat.close()

    graphH_path = os.path.join(output_dir, "eigvals_hist.png")
    graphH = pdf.DfEnergyLevelHistoryGraphH()
    graphH.set_HOMO_level(HOMO_level)  # option base 0
    graphH.load_data(data_path)
    graphH.save(graphH_path)

    # lastのグラフを作成する
    graphV_path = os.path.join(output_dir, "eigvals_last.png")
    graphV = pdf.DfEnergyLevelHistoryGraphV()
    graphV.set_HOMO_level(HOMO_level)  # option base 0
    graphV.set_LUMO_level(HOMO_level + 1)  # option base 0
    graphV.load_data(data_path)
    if entry.scf_converged:
        graphV.select_iterations([itr])
    else:
        graphV.select_iterations([itr - 1])
    graphV.ylabel = ""
    graphV.is_draw_grid = False
    graphV.y_ticks = [-1]
    graphV.y_ticklabels = [""]
    graphV.save(graphV_path)


# def plot_energy_level(entry, output_path):
#    graph = EnergyLevelSingle()
#
#    HOMO_level = entry.get_HOMO_level('ALPHA') # TODO
#
#    data_path = output_dir + '/eigvals_last.dat'
#    graph_path = output_dir + '/eigvals_last.png'
#
#    dat = open(data_path, 'w')
#    last_it = entry.iterations
#    eigvals = entry.get_energylevel('ALPHA', itr) # TODO
#    if eigvals:
#        for level, e in enumerate(eigvals):
#            e *= 27.2116
#            dat.write('%d, %d, % 16.10f\n' % (itr, level, e))
#    dat.close()
#
#    graph = pdf.DfEigValsHistGraph()
#    graph.set_HOMO_level(HOMO_level) # option base 0
#    graph.load_data(data_path)
#    graph.save(graph_path)


# def get_eigenvalues(iterations):
#    eigval_vtr_path = 'fl_Work/eigenvalues.rks%d.vtr' % (iterations) # TODO
#    v = pdf.Vector()
#    if v.is_loadable(eigval_vtr_path):
#        v.load(eigval_vtr_path)
#    else:
#        logging.warning('cannot load: %s\n' % (eigval_vtr_path))
#    return v


def get_rst(pdfdata):
    assert isinstance(pdfdata, (pdf.PdfArchive, pdf.PdfParam_H5))
    contents = """
    =========================
    ProteinDF report:
    =========================

    ======================== ================
    item                     number
    ======================== ================
    the number of AOs        {num_of_AOs:>16d}
    the number of MOs        {num_of_MOs:>16d}
    the number of iterations {num_of_iterations:>16d}
    Total energy             {total_energy:>16.8f}
    ======================== ================
    """
    num_of_iterations = pdfdata.iterations
    return contents.format(
        num_of_AOs=pdfdata.num_of_AOs,
        num_of_MOs=pdfdata.num_of_MOs,
        num_of_iterations=num_of_iterations,
        total_energy=pdfdata.get_total_energy(num_of_iterations),
    )


if __name__ == "__main__":
    if os.path.exists("config.ini"):
        logging.config.fileConfig("config.ini", disable_existing_loggers=False)
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import textwrap

from .pdfgraph import DfTotalEnergyHistGraph, DfEnergyLevelHistoryGraphH, DfEnergyLevelHistoryGraphV, DfPopulationGraph
from .pdfparam_object import PdfParamObject

import logging
logger = logging.getLogger(__name__)

# , DfEigValsHistGraph


class PdfReport(object):
    def __init__(self, pdfparam, workdir="./report"):
        assert(isinstance(pdfparam, PdfParamObject))
        self._pdfparam = pdfparam
        self._workdir = workdir

    def _make_work_dir(self):
        if not os.path.isdir(self._workdir):
            os.mkdir(self._workdir)

    def report(self, iteration=0):
        if iteration == 0:
            iteration = self._pdfparam.iterations
        logger.debug('report iteration={}'.format(iteration))

        self._plot_convergence_check(iteration)
        
        for run_type in self._pdfparam.run_types():
            logger.info("make report [run_type: {}]".format(run_type))
            self._plot_convergence_energy_level(run_type, iteration)

        # plot_energy_level(pdfdata,
        #                  workdir + '/level.png')

        # plot_elapsed_time(pdfdata,
        #                  workdir + '/elapsed_time.png')

        # if self._pdfparam.scf_converged:
        # self._plot_pop_mulliken(iteration)

        contents = self._get_rst()
        print(contents)

    def _get_total_energy_vector(self):
        """
        Total Energyのヒストリ
        """
        num_of_iterations = self._pdfparam.iterations
        TEs = [None] * num_of_iterations
        for itr in range(num_of_iterations):
            TEs[itr] = self._pdfparam.get_total_energy(itr +1)
        return TEs

    def _get_rst(self):
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
        contents = textwrap.dedent(contents).strip()

        num_of_iterations = self._pdfparam.iterations
        TEs = self._get_total_energy_vector()
        answer = contents.format(
            num_of_AOs=self._pdfparam.num_of_AOs,
            num_of_MOs=self._pdfparam.num_of_MOs,
            num_of_iterations=num_of_iterations,
            total_energy=TEs[num_of_iterations -1])

        return answer

    def _plot_convergence_check(self, iteration):
        """
        TEのヒストリ
        """
        iterations = self._pdfparam.iterations

        # TE
        TEs = self._get_total_energy_vector()

        data_path = os.path.join(self._workdir, "TE_hist.dat")
        self._make_work_dir()
        with open(data_path, 'w') as dat:
            for i in range(iterations):
                dat.write('%d, % 16.10f\n' % (i +1, TEs[i]))

        graph = DfTotalEnergyHistGraph()
        graph.load_data(data_path)
        graph_path = os.path.join(self._workdir, "TE_hist.png")
        graph.save(graph_path)

    def _plot_convergence_energy_level(self, run_type, iteration):
        """
        EnergyLevelのヒストリ
        """
        HOMO_level = self._pdfparam.get_HOMO_level(run_type)

        data_path = os.path.join(self._workdir, "eigvals_hist.{}.dat".format(run_type))
        self._make_work_dir()
        with open(data_path, 'w') as dat:
            for itr in range(1, iteration + 1):
                eigvals = self._pdfparam.get_energy_level(run_type, itr)
                if eigvals:
                    for level, e in enumerate(eigvals):
                        e *= 27.2116
                        dat.write('%d, %d, % 16.10f\n' % (itr, level, e))

        graphH = DfEnergyLevelHistoryGraphH()
        graphH.set_HOMO_level(HOMO_level)  # option base 0
        graphH.load_data(data_path)
        graphH_path = os.path.join(self._workdir, "eigvals_hist.{}.png".format(run_type))
        graphH.save(graphH_path)

        # lastのグラフを作成する
        graphV = DfEnergyLevelHistoryGraphV()
        graphV.set_HOMO_level(HOMO_level)  # option base 0
        graphV.set_LUMO_level(HOMO_level + 1)  # option base 0
        graphV.load_data(data_path)
        graphV.select_iterations([iteration])
        graphV.ylabel = ''
        graphV.is_draw_grid = False
        graphV.y_ticks = [-1]
        graphV.y_ticklabels = [""]
        graphV_path = os.path.join(self._workdir, "eigvals_last.{}.png".format(run_type))
        graphV.save(graphV_path)

    def _plot_pop_mulliken(self, iteration):
        method = self._pdfparam.method
        run_type = "rks"

        pop = self._pdfparam.get_pop_mulliken_atom(run_type, iteration)
        pop_data_path = os.path.join(self._workdir, "mulliken.dat")
        self._make_plot_pop_mulliken_data(pop, pop_data_path)

        graph1 = DfPopulationGraph()
        graph1.load_data(pop_data_path)
        graph1_path = os.path.join(self._workdir, "mulliken.png")
        graph1.save(graph1_path)

    def _make_plot_pop_mulliken_data(self, pop_vtr, output_path):
        num_of_items = len(pop_vtr)
        with open(output_path, "w") as f:
            for index in range(num_of_items):
                charge = pop_vtr[index]
                f.write("{}, {}\n".format(index + 1, charge))

    # def plot_energy_level(self):
    #     HOMO_level = entry.get_HOMO_level('ALPHA') # TODO
    #
    #     data_path = os.path.join(self._workdir, "eigvals_last.dat")
    #     graph_path = os.path.join(self._workdir, "eigvals_last.png")
    #
    #     with open(data_path, 'w') as dat:
    #         last_it = entry.iterations
    #         eigvals = entry.get_energylevel('ALPHA', itr) # TODO
    #         if eigvals:
    #             for level, e in enumerate(eigvals):
    #                 e *= 27.2116
    #                 dat.write('%d, %d, % 16.10f\n' % (itr, level, e))
    #
    #     graph = DfEigValsHistGraph()
    #     graph.set_HOMO_level(HOMO_level) # option base 0
    #     graph.load_data(data_path)
    #     graph.save(graph_path)


# def plot_convergence_check0(pdfdata, output_path):
#     graph = pdfex.GraphConvergenceCheck()
#
#     graph.set_xlabel('SCF convergence step')
#     graph.set_ylabel('difference / a.u.')
#     graph.draw_grid(True)
#     graph.draw_legend(True)
#
#     itr = pdfdata.get_num_of_iterations()
#     # total energy
#     TE_history = range(itr +1)
#     TE_history[0] = None
#     # delta TE
#     deltaTE_history = range(itr +1)
#     deltaTE_history[0] = None
#     # delta Density Matrix
#     deltaDensMat_history = range(itr +1)
#     deltaDensMat_history[0] = None
#     for itr in range(1, iterations +1):
#         TE_history[itr] = dfdata.get_total_energy(itr)
#         info = dfdata.get_convergence_info(itr)
#         deltaTE_history[itr] = info.get('max_deviation_of_total_energy', None)
#         deltaDensMat_history[itr] = info.get('max_deviation_of_density_matrix', None)
#
#     #graph.plot(TE_history, "Total Energy")
#     graph.plotLog(deltaTE_history, "delta Total Energy")
#     graph.plotLog(deltaDensMat_history, "delta Density Matrix")
#
#     graph.prepare()
#     graph.file_out(output_path)
#
#
# def plot_elapsed_time(pdfdata, output_path):
#     graph = pdfex.BarGraph()
#
#     graph.set_xlabel('SCF convergence step')
#     graph.set_ylabel('elapsed time / sec')
#
#     # prepare data
#     num_of_iterations = pdfdata.get_num_of_iterations()
#     data = [None] * num_of_iterations
#     scf_sections = pdfdata.get_stat_scf_sections()
#     for itr in range(num_of_iterations):
#         contents = [None] * len(scf_sections)
#         for entry, scf_section in enumerate(scf_sections):
#             contents[entry] = pdfdata.get_elapsed_time(scf_section, itr +1)
#         data[itr] = contents
#
#     graph.plot(data, scf_sections)
#
#     graph.prepare()
#     graph.file_out(output_path)

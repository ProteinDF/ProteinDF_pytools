#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import textwrap
import yaml

from .pdfgraph import DfTotalEnergyHistGraph, DfEnergyLevelHistoryGraphH, DfEnergyLevelHistoryGraphV, DfPopulationGraph
from .pdfparam_object import PdfParamObject

import logging
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


def _prepare_report_data_for_yaml(data):
    if not isinstance(data, dict):
        return data
    te_hist = data.get("total_energy_history")
    if isinstance(te_hist, list):
        data["total_energy_history"] = [
            p if isinstance(p, _FlowMapping) else _FlowMapping(p)
            for p in te_hist
            if isinstance(p, dict)
        ]
    e_levels = data.get("energy_levels")
    if isinstance(e_levels, list):
        for el in e_levels:
            if isinstance(el, dict) and "history" in el:
                hist = el.get("history")
                if isinstance(hist, list):
                    el["history"] = [
                        p if isinstance(p, _FlowMapping) else _FlowMapping(p)
                        for p in hist
                        if isinstance(p, dict)
                    ]
    return data


def save_report_yaml(data, output_path):
    data = _prepare_report_data_for_yaml(data)
    content = yaml.safe_dump(data, sort_keys=False, allow_unicode=True)
    out_dir = os.path.dirname(output_path)
    if out_dir and not os.path.exists(out_dir):
        os.makedirs(out_dir, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        f.write(content)


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

        for run_type in self._pdfparam.orbital_run_types():
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

    def get_report_data(self, iteration=0, source=None):
        if iteration == 0:
            iteration = self._pdfparam.iterations

        num_of_iterations = self._pdfparam.iterations
        TEs = self._get_total_energy_vector()
        if (
            num_of_iterations > 0
            and len(TEs) >= num_of_iterations
            and TEs[num_of_iterations - 1] is not None
        ):
            total_energy = float(TEs[num_of_iterations - 1])
        else:
            total_energy = None

        total_energy_history = []
        for i in range(num_of_iterations):
            if TEs[i] is not None:
                total_energy_history.append(
                    _FlowMapping(
                        {
                            "iteration": i + 1,
                            "energy": float(TEs[i]),
                        }
                    )
                )

        energy_levels = []
        for run_type in self._pdfparam.orbital_run_types():
            homo_level = self._pdfparam.get_HOMO_level(run_type)
            history = []
            for itr in range(1, iteration + 1):
                eigvals = self._pdfparam.get_energy_level(run_type, itr)
                if eigvals:
                    for level, e in enumerate(eigvals):
                        history.append(
                            _FlowMapping(
                                {
                                    "iteration": itr,
                                    "level": level,
                                    "energy_ev": float(e) * 27.2116,
                                }
                            )
                        )

            last_eigvals = self._pdfparam.get_energy_level(run_type, iteration)
            if last_eigvals:
                last = [float(e) * 27.2116 for e in last_eigvals]
            else:
                last = []

            energy_levels.append(
                {
                    "run_type": run_type,
                    "HOMO_level": homo_level,
                    "history": history,
                    "last": last,
                }
            )

        if source is None:
            source = getattr(self._pdfparam, "filepath", "")

        return {
            "schema": "pdf-report/1",
            "kind": "report",
            "source": source,
            "iteration": iteration,
            "summary": {
                "num_of_AOs": int(self._pdfparam.num_of_AOs),
                "num_of_MOs": int(self._pdfparam.num_of_MOs),
                "num_of_iterations": num_of_iterations,
                "total_energy": total_energy,
            },
            "total_energy_history": total_energy_history,
            "energy_levels": energy_levels,
        }

    save_report_yaml = staticmethod(save_report_yaml)

    def _get_total_energy_vector(self):
        """
        History of the total energy.
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
        Total energy             {total_energy:>16}
        ======================== ================
        """
        contents = textwrap.dedent(contents).strip()

        num_of_iterations = self._pdfparam.iterations
        TEs = self._get_total_energy_vector()
        if num_of_iterations > 0 and len(TEs) >= num_of_iterations:
            total_energy_str = "{:>16.8f}".format(TEs[num_of_iterations - 1])
        else:
            total_energy_str = "{:>16}".format("n/a")

        answer = contents.format(
            num_of_AOs=self._pdfparam.num_of_AOs,
            num_of_MOs=self._pdfparam.num_of_MOs,
            num_of_iterations=num_of_iterations,
            total_energy=total_energy_str,
        )

        return answer

    def _plot_convergence_check(self, iteration):
        """
        History of TE.
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
        History of the energy level.
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

        # create the graph for the last iteration
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
        self._make_work_dir()
        for run_type in self._pdfparam.run_types():
            pop = self._pdfparam.get_pop_mulliken_atom(run_type, iteration)
            if pop is not None:
                pop_data_path = os.path.join(
                    self._workdir, "mulliken.{}.dat".format(run_type)
                )
                self._make_plot_pop_mulliken_data(pop, pop_data_path)

                graph1 = DfPopulationGraph()
                graph1.load_data(pop_data_path)
                graph1_path = os.path.join(
                    self._workdir, "mulliken.{}.png".format(run_type)
                )
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

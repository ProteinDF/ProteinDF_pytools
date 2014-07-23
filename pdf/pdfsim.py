#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import copy
import math
import logging

import bridge
import pdf

class PdfSim(object):
    """
    """
    def __init__(self):
        nullHandler = bridge.NullHandler()
        self._logger = logging.getLogger(__name__)
        self._logger.addHandler(nullHandler)

    def make_workdir(self, workdir ='.'):
        # make sub-directories
        dirs = ['fl_Work']
        for d in dirs:
            path = os.path.join(workdir, d)
            if not os.path.exists(path):
                os.mkdir(path)

        
    def setup(self, pdfparam =pdf.get_default_pdfparam(), workdir ="."):
        """
        setup to run ProteinDF
        """
        # make fl_Userinput
        assert(isinstance(pdfparam, pdf.PdfParam))

        if not os.path.exists(workdir):
            os.mkdir(workdir)
        else:
            self._logger.info('already exist: {}'.format(workdir))

        input_path = os.path.join(workdir, "fl_Userinput")
        f = open(input_path, 'w')
        f.write(pdfparam.get_inputfile_contents())
        f.close()

        self.make_workdir(workdir)

    def run_pdf(self, subcmd):
        """
        run ProteinDF command
        """
        if isinstance(subcmd, str):
            subcmd = [subcmd]
            
        cmd = os.path.join(pdf_home(), "bin", "pdf")
        cmdlist = [cmd]
        cmdlist.extend(subcmd)
        self._logger.debug("run: {0}".format(cmdlist))
        
        try:
            subprocess.check_call(cmdlist)
        except:
            print('-'*60)
            traceback.print_exc(file=sys.stdout)
            #print(traceback.format_exc())
            print('-'*60)


    # ------------------------------------------------------------------
    def sp(self, pdfparam, *args, **kwargs):
        """
        calc single-point

        Keyword arguments:
        pdfparam --- PdfParam object represented the calculation condition.
        workdir --- working directory.
        db_path --- ProteinDF ArchiveDB path (default: pdfresults.db)
        dry_run --- if True, the calculation is NOT carried out. Default is False.

        Returns:
        tuple of the number of iterations and the total energy.
        
        """
        workdir = kwargs.get('workdir', '')
        db_path = kwargs.get('db_path', 'pdfresults.db')
        dry_run = kwargs.get('dry_run', False)
        
        current_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or ".")

        pdf_workdir = os.path.join('.', workdir)
        if not os.path.exists(pdf_workdir):
            os.mkdir(pdf_workdir)
        else:
            self._logger.info('already exist {}'.format(pdf_workdir))

        self.setup(pdfparam, pdf_workdir)
        os.chdir(pdf_workdir)

        itr = None
        total_energy = None
        if not dry_run:
            self.run_pdf(['-o', 'pdf.log', 'serial'])
            self.run_pdf('archive')

            entry = pdf.PdfArchive(db_path)
            itr = entry.iterations
            total_energy = entry.get_total_energy(itr)

        os.chdir(current_dir)

        return (itr, total_energy)
    # ------------------------------------------------------------------





                

    def opt(self, pdfparam, workdir=".", max_cycle=100):
        current_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or ".")

        alpha = 0.5
        accuracy = 1.0E-3

        opt_cycle = 1
        while True:
            opt_workdir = os.path.join(current_dir, "opt_cycle_{}".format(opt_cycle))
            if not os.path.exists(opt_workdir):
                os.mkdir(opt_workdir)
            else:
                self._logger.info('already exist: {}'.format(opt_workdir))
            (grad_mat, rms) = self.numerical_grad(pdfparam, opt_workdir, accuracy)

            max_force = 0.0
            for i in range(len(grad_mat)):
                for j in range(3):
                    max_force = max(max_force, math.fabs(grad_mat[i][j]))

            if (max_force < 0.001 and rms < 0.001) or (opt_cycle > max_cycle):
                break

            # update
            i = 0
            for atom_id, atom in pdfparam.molecule.atoms():
                atom.shift_by(alpha * grad_mat[i][0],
                              alpha * grad_mat[i][1],
                              alpha * grad_mat[i][2])

            opt_cycle += 1

        self._logger.info("opt done")
        self._logger.info(str(pdfparam.molecule))

    def numerical_grad(self, pdfparam, workdir=".", accuracy=1.0E-3, delta=0.001):
        direction_str = ["x", "y", "z"]

        self._logger.debug(">>>> start numerical grad")
        molecule = pdfparam.molecule
        num_of_atoms = molecule.get_number_of_all_atoms()

        if pdfparam.convergence_threshold_energy < accuracy * 0.5:
            self._logger.warning('convergence_thresold_energy={} > 0.5*accuracy(={})'.format(
                pdfparam.convergence_threshold_energy,
                accuracy*0.5))

        grad_mat = [[0.0, 0.0, 0.0] for x in range(num_of_atoms)]
        index = 0
        for atom_id, atom in molecule.atoms():
            for direction in range(3):
                self._logger.debug("start numerical grad for {}".format(direction_str[direction]))
                h = delta

                delta_TE = 0.0
                while True:
                    atom1 = self._move(atom, direction, -h)
                    atom2 = self._move(atom, direction, +h)

                    pdfparam1 = pdfparam
                    pdfparam1.molecule[atom_id] = atom1
                    workdir1 = os.path.join(workdir, "{}_d{}_h{}_1".format(atom_id, direction,h))
                    (itr, total_energy1) = self.calc_pdf(pdfparam1, workdir1)

                    pdfparam2 = pdfparam
                    pdfparam2.molecule[atom_id] = atom2
                    workdir2 = os.path.join(workdir, "{}_d{}_h{}_2".format(atom_id, direction,h))
                    (itr, total_energy2) = self.calc_pdf(pdfparam2, workdir2)

                    delta_TE = total_energy2 - total_energy1
                    print("h={: e}, delta_TE={: e}, v={: e}".format(h, delta_TE, delta_TE / (2.0 * h)))
                    if (math.fabs(delta_TE) < accuracy) or (h < 1.0E-2):
                        self._logger.debug("numerical grad condition satisfied. value={} < {}".format(delta_TE, accuracy))
                        break

                    h *= 0.5
                    sys.stderr.flush()
                    sys.stdout.flush()

                # x,y,zのどれかが終了
                value = delta_TE / (2.0 * h)
                grad_mat[index][direction] = float(value)
                self._logger.debug("gradient value [{}][{}] = {}".format(atom_id,
                                                                         direction_str[direction],
                                                                         value))
                #self._show_grad_mat(grad_mat)

                sys.stderr.flush()
                sys.stdout.flush()
            index += 1

        self._logger.info("=== grad (accuracy={}) ===".format(accuracy))
        rms = 0.0
        index = 0
        for atom_id, atom in molecule.atoms():
            self._logger.info("[{:>8s}] {: 8.5f} {: 8.5f} {: 8.5f}".format(atom_id, grad_mat[index][0], grad_mat[index][1], grad_mat[index][2]))
            for i in range(3):
                rms += grad_mat[index][i] * grad_mat[index][i]
            index += 1
        rms /= num_of_atoms * 3
        rms = math.sqrt(rms)
        self._logger.info("RMS = {}".format(rms))
        self._logger.info("============\n")

        return (grad_mat, rms)

    def total_energy(self, pdfparam, workdir="."):
        (itr, total_energy) = self.calc_pdf(pdfparam, workdir)

        return (itr, total_energy)

    def _show_grad_mat(self, grad_mat):
        for i in range(len(grad_mat)):
            print("[{}] ({}, {}, {})".format(i, grad_mat[i][0], grad_mat[i][1], grad_mat[i][2]))

    def _move(self, atom, direction, delta):
        assert(0 <= direction < 3)
        #answer = copy.deepcopy(atom)
        answer = bridge.Atom(atom)

        if direction == 0:
            answer.xyz.x = float(answer.xyz.x) + delta
        elif direction == 1:
            answer.xyz.y = float(answer.xyz.y) + delta
        elif direction == 2:
            answer.xyz.z = float(answer.xyz.z) + delta
        else:
            self._logger.critical("program error.")
            exit(1)

        return answer

    def calc_pdf(self, pdfparam, workdir):
        current_dir = os.path.abspath(os.path.dirname(sys.argv[0]) or ".")

        pdf_workdir = os.path.join('.', workdir)
        if not os.path.exists(pdf_workdir):
            os.mkdir(pdf_workdir)
        else:
            self._logger.info('already exist {}'.format(pdf_workdir))

        self.setup(pdfparam, pdf_workdir)
        os.chdir(pdf_workdir)
        self.run_pdf(['-o', 'pdf.log', 'serial'])
        self.run_pdf('archive')

        entry = pdf.PdfArchive('pdfresults.db')
        itr = entry.iterations
        total_energy = entry.get_total_energy(itr)

        os.chdir(current_dir)

        return (itr, total_energy)


        

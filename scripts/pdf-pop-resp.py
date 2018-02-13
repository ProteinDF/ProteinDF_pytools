#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import numpy as np
import pandas as pd
try:
    import msgpack
except:
    import msgpack_pure as msgpack
import math

import proteindf_bridge as bridge
import proteindf_tools as pdf

class ESP_charge(object):
    def get_RRMS(self, mpac_path, atoms):
        grids, ESPs = self._load_ESPs(mpac_path)
        rrms = self._calcRRMS(grids, ESPs, atoms)

        return rrms

    def _load_ESPs(self, mpac_path):
        print('load: {}'.format(mpac_path))
        data = None
        with open(mpac_path, 'rb') as f:
            mpac_data = f.read()
            data = msgpack.unpackb(mpac_data)
            data = bridge.Utils.to_unicode_dict(data)

            grids = []
            for p in data['grids']:
                pos = bridge.Position(p[0], p[1], p[2])
                pos *= (1.0 / 0.5291772108) # Angstroam to a.u.
                grids.append(pos)

            ESPs = []
            for v in data['ESP']:
                ESPs.append(v)

            assert(len(grids) == len(ESPs))
        return (grids, ESPs)

    def _calcRRMS(self, grids, ESPs, atoms):
        sum_delta2 = 0.0
        sum_v2 = 0.0
        num_of_grids = len(grids)
        print('# of grids: {}'.format(num_of_grids))
        assert(len(ESPs) == num_of_grids)
        for i in range(num_of_grids):
            estimate_esp = self._calc_esp(grids[i], atoms)
            exact_esp = ESPs[i]
            delta = estimate_esp - exact_esp
            delta2 = delta * delta
            sum_delta2 += delta2
            sum_v2 += exact_esp * exact_esp
            #print("grid={} est={: 8.3e} exact={: 8.3e} >> delta2={: 8.3e}".format(grids[i], estimate_esp, exact_esp, delta2))
        rrms2 = sum_delta2 / sum_v2
        rrms = math.sqrt(rrms2)
        print('delta2={} sum_v2={} RRMS2={} RRMS={}'.format(sum_delta2, sum_v2, rrms2, rrms))

        return rrms

    def _calc_esp(self, pos, atoms):
        esp = 0.0
        num_of_atoms = len(atoms)
        for i in range(num_of_atoms):
            d = pos.distance_from(atoms[i].xyz)
            esp += atoms[i].charge / d;

        return esp



class Resp(object):
    def __init__(self, alpha=0.0005, beta=0.1, refs=None, max_iter=1000, tol=0.0001):
        self._alpha = float(alpha)
        self._beta = float(beta)
        self._refs = refs
        self._max_iter = int(max_iter)
        self._iter = 0

        self._conv_threshold = float(tol)
        self._num_of_converged = 0
        self._is_converged = False
        self._coef = None

    @property
    def iterations(self):
        return self._iter

    @property
    def is_converged(self):
        return self._is_converged

    @property
    def coef(self):
        return self._coef

    def _is_convergence(self, iter, x, prev_x):
        answer = False
        if iter > 1:
            dim = len(x)

            diff_x = x - prev_x
            max_diff = max(abs(diff_x.max), abs(diff_x.min))
            rms_diff = diff_x * diff_x / dim
            print("#{} MAX delta: {} MAX RMS: {}".format(iter, max_diff, rms_diff))

            if (max_diff < self._conv_threshold * 0.1):
                self._num_of_converged += 1
                if self._num_of_converged >= 2:
                    answer = True
            else:
                self._num_of_converged = 0

        return answer

    def fit_q(self, model_mat, predicted):
        print('>>> fit RESP charges')
        dim = model_mat.rows
        assert(dim == model_mat.cols)
        assert(predicted.size() == dim)

        # initial guess
        coef = bridge.Vector(dim)
        prev_coef = bridge.Vector(coef)
        if self._refs == None:
            self._refs = bridge.Vector(dim)

        for self._iter in range(1, self._max_iter +1):
            A = bridge.Matrix(model_mat)
            y = bridge.Vector(predicted)
            for i in range(dim -1):
                v = A.get(i, i)
                A.set(i, i, v - 2.0 * A.self._alpha * (coef[i] - self._refs[i]))

            invA = A.inverse()
            coef = invA * y

            if self._is_convergence(self._iter, coef, prev_coef):
                self._is_converged = True
                break
            prev_coef = bridge.Vector(coef)

        self._coef = coef


    def fit_h(self, model_mat, predicted):
        print('>>> fit RESP charges (hyperbolic)')
        dim = model_mat.rows
        assert(dim == model_mat.cols)
        assert(predicted.size() == dim)

        # initial guess
        coef = bridge.Vector(dim)
        prev_coef = bridge.Vector(coef)
        if self._refs == None:
            self._refs = bridge.Vector(dim)

        b = self._beta
        b2 = b * b
        for self._iter in range(1, self._max_iter +1):
            A = bridge.Matrix(model_mat)
            y = bridge.Vector(predicted)
            for i in range(dim -1):
                v = A.get(i, i)
                q = coef[i]
                q2 = q * q
                v += self._alpha * math.sqrt(q2 + b2)
                A.set(i, i, v)

            invA = A.inverse()
            coef = invA * y

            if self._is_convergence(self._iter, coef, prev_coef):
                self._is_converged = True
                break
            prev_coef = bridge.Vector(coef)

        self._coef = coef


def get_charge(x):
    charge = 0.0
    for i in range(len(x) -1):
        charge += x[i]
    return charge


def set_atomlist(atom_charges, atomlist):
    q_total = 0.0
    for i in range(len(atom_charges)):
        q = atom_charges[i];
        atomlist[i].charge = q;
        q_total += q
        print(atomlist[i])
    print('charge={:.3f}'.format(q_total))


def main():
    # parse args
    parser = argparse.ArgumentParser(description='lasso')

    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-d', '--db',
                       nargs='?',
                       action='store',
                       const='pdfresults.db',
                       help='ProteinDF results file')
    group.add_argument('-p', '--param',
                       nargs='?',
                       action='store',
                       const='pdfparam.mpac',
                       help='ProteinDF parameter file')

    parser.add_argument("-m", "--max_iterations",
                        type=int,
                        nargs=1,
                        action='store',
                        default=['1000'],
                        help="ridge alpha")
    parser.add_argument("-t", "--tolerance",
                        type=float,
                        nargs=1,
                        action='store',
                        default=['0.0001'],
                        help="ridge alpha")
    parser.add_argument("--alpha",
                        type=float,
                        nargs=1,
                        action='store',
                        default=['0.0005'],
                        help="RESP alpha")
    parser.add_argument("--beta",
                        type=float,
                        nargs=1,
                        action='store',
                        default=['0.1'],
                        help="ridge alpha")
    parser.add_argument("-e", "--espdat",
                        nargs=1,
                        action='store',
                        default=['grid-esp.mpac'],
                        help="ridge alpha")
    parser.add_argument("--csv",
                        type=str,
                        nargs=1,
                        action='store',
                        default=[''],
                        help="output csv path")

    parser.add_argument("design_matrix_path",
                        action='store',
                        help="design matrix path")
    parser.add_argument("target_vector_path",
                        action='store',
                        help="target vector path")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    # setting
    design_matrix_path = args.design_matrix_path
    target_vector_path = args.target_vector_path
    output_csv_path = args.csv[0]
    verbose = args.verbose
    max_iter = args.max_iterations[0]
    tol = args.tolerance[0] # 0.0001
    alpha = float(args.alpha[0])
    beta = float(args.beta[0])
    if verbose:
        print("alpha={}".format(alpha))
        print("beta={}".format(beta))
    esp_data_path = args.espdat[0]

    atomlist = []
    if args.db:
        entry = pdf.PdfArchive(args.db)
    elif args.param:
        pdfparam = pdf.load_pdfparam(args.param)
        atoms = pdfparam.molecule.get_atom_list()
        for atom in atoms:
            if atom.symbol != 'X':
                atomlist.append(atom)

    # load design matrix
    # design matrix size = (atom +1) * (atom +1): +1 for lagurange parameter of MK
    X = pdf.Matrix()
    if pdf.Matrix.is_loadable(design_matrix_path):
        X.load(design_matrix_path)
    elif pdf.SymmetricMatrix.is_loadable(design_matrix_path):
        tmp = pdf.SymmetricMatrix()
        tmp.load(design_matrix_path)
        X = tmp.get_general_matrix()
    else:
        raise
    y = pdf.Vector()
    y.load(target_vector_path)

    atom_charges = []
    resp = Resp(alpha=alpha, beta=beta, refs=None, max_iter=max_iter, tol=tol)
    resp.fit_h(X, y)
    print(resp.iterations)
    print(resp.is_converged)
    atom_charges = resp.coef[:-1]
    set_atomlist(atom_charges, atomlist)

    # RRMS
    esp_charge = ESP_charge()
    rrms = esp_charge.get_RRMS(esp_data_path, atomlist)
    print("RRMS={:.3f}".format(rrms))

    # output csv
    if output_csv_path != '':
        with open(output_csv_path, 'w') as f:
            for i in range(len(atomlist)):
                line = "{:2}, {: 8.3f}\n".format(atomlist[i].symbol,
                                                 atomlist[i].charge)
                f.write(line)

if __name__ == '__main__':
    main()

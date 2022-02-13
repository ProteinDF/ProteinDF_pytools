#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import math

import proteindf_bridge as bridge
import proteindf_tools as pdf

import logging
import logging.config
logger = logging.getLogger(__name__)


class Resp(object):
    def __init__(self, alpha=0.0005, beta=0.1, refs=None, max_iter=1000, threshold=1e-3, rms_threshold=1e-2):
        self._alpha = float(alpha)
        self._beta = float(beta)
        self._refs = refs
        self._max_iter = int(max_iter)
        self._iter = 0

        self._max_threshold = float(threshold)
        self._rms_threshold = float(rms_threshold)

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
            logger.info("#{} MAX delta: {} MAX RMS: {}".format(
                iter, max_diff, rms_diff))

            if (max_diff < self._max_threshold) and (rms_diff < self._rms_threshold):
                answer = True

        return answer

    def fit_q(self, model_mat, predicted):
        logger.info('>>> fit RESP charges (quadric) a={}'.format(self._alpha))
        dim = model_mat.rows
        assert(dim == model_mat.cols)
        assert(predicted.size() == dim)

        # initial guess
        coef = bridge.Vector(dim)
        prev_coef = bridge.Vector(coef)
        if self._refs == None:
            self._refs = bridge.Vector(dim)

        for self._iter in range(1, self._max_iter + 1):
            A = bridge.Matrix(model_mat)
            y = bridge.Vector(predicted)
            for i in range(dim - 1):
                v = A.get(i, i)
                rest = -2.0 * self._alpha * (self._refs[i] - coef[i])
                A.set(i, i, v + rest)
                y[i] += self._refs[i] * rest

            invA = A.inverse()
            coef = invA * y

            if self._is_convergence(self._iter, coef, prev_coef):
                self._is_converged = True
                break
            prev_coef = bridge.Vector(coef)

        self._coef = coef

    def fit_h(self, model_mat, predicted):
        logger.info('>>> fit RESP charges (hyperbolic) a={} b={}'.format(
            self._alpha, self._beta))
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
        for self._iter in range(1, self._max_iter + 1):
            A = bridge.Matrix(model_mat)
            y = bridge.Vector(predicted)
            for i in range(dim - 1):
                v = A.get(i, i)
                q = coef[i]
                q2 = q * q
                rest = self._alpha * q * (1.0 / math.sqrt(q2 + b2))
                A.set(i, i, v + rest)
                y[i] += self._refs[i] * rest

            invA = A.inverse()
            coef = invA * y

            if self._is_convergence(self._iter, coef, prev_coef):
                self._is_converged = True
                break
            prev_coef = bridge.Vector(coef)

        self._coef = coef


def main():
    # parse args
    parser = argparse.ArgumentParser(description='calculate RESP charge')

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
                        help="max iteration")

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

    parser.add_argument("--use-quadratic",
                        action='store_true',
                        default=False,
                        help="use quadratic function")

    parser.add_argument("--design_matrix_path",
                        nargs=1,
                        action='store',
                        default=['esp_design.mat'],
                        help="design matrix path (default: esp_design.mat), which is created by `pdf-pop-esp`.")
    parser.add_argument("--target_vector_path",
                        nargs=1,
                        action='store',
                        default=['esp_target.vtr'],
                        help="target vector path (default: esp_target.vtr), which is created by `pdf-pop-esp`.")
    parser.add_argument("--model_vector_path",
                        nargs=1,
                        action='store',
                        default=['resp_model.vtr'],
                        help="output model vector path (default: resp_model.vtr).")

    parser.add_argument("-e", "--espdat",
                        nargs=1,
                        action='store',
                        default=['grid-esp.mpac'],
                        help="ESP values on grids by msgpack format (default: grid-esp.mpac)")

    parser.add_argument("--csv",
                        type=str,
                        nargs=1,
                        action='store',
                        default=[''],
                        help="output csv path")

    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    args = parser.parse_args()

    # logging
    logging.basicConfig(level=logging.INFO)

    # setting
    design_matrix_path = args.design_matrix_path[0]
    target_vector_path = args.target_vector_path[0]
    model_vector_path = args.model_vector_path[0]
    output_csv_path = args.csv[0]
    verbose = args.verbose
    max_iter = args.max_iterations[0]
    # tol = args.tolerance[0] # 0.0001
    alpha = float(args.alpha[0])
    beta = float(args.beta[0])
    use_quadratic = args.use_quadratic
    print("use_quadratc: ", use_quadratic)
    if verbose:
        logger.info("alpha={}".format(alpha))
        logger.info("beta={}".format(beta))
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
    resp = Resp(alpha=alpha, beta=beta, refs=None, max_iter=max_iter)
    if use_quadratic:
        resp.fit_q(X, y)
    else:
        resp.fit_h(X, y)
    logger.info("iterations: {}".format(resp.iterations))
    if resp.is_converged:
        logger.debug("well converged.")
    else:
        logger.warning("NOT converged.")
    atom_charges = resp.coef[:-1]
    pdf.PopEspUtils.set_atomlist(atom_charges, atomlist)

    # output
    if len(model_vector_path) > 0:
        atom_charges = pdf.Vector(atom_charges)
        atom_charges.save(model_vector_path)

    # output csv
    if output_csv_path != '':
        with open(output_csv_path, 'w') as f:
            for i in range(len(atomlist)):
                line = "{:2}, {: 8.3f}\n".format(atomlist[i].symbol,
                                                 atomlist[i].charge)
                f.write(line)

    # RRMS
    espUtils = pdf.PopEspUtils()
    rrms = espUtils.get_RRMS(esp_data_path, atomlist)
    print("RRMS={: .5f}".format(rrms))


if __name__ == '__main__':
    main()

#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import logging
import logging.config
logger = logging.getLogger(__name__)

import proteindf_bridge as bridge
import proteindf_tools as pdf


class Ridge(object):
    def __init__(self, alpha=1.0):
        self._alpha = alpha
        self._coef = None

    def fit(self, X, y):
        '''
        X: design matrix (n x b)
        '''
        n = X.rows
        b = X.cols
        assert(y.size() == n)

        tX = pdf.Matrix(X)
        tX.transpose()
        tXX = tX * X

        rI = bridge.identity_matrix(b)
        rI *= self._alpha

        M = (tXX + rI)
        Minv = M.inverse();

        MinvX = Minv * tX
        self._coef = MinvX * y

    def fit_lagrangian(self, X, y):
        '''
        X: design matrix (n x b)
        '''
        print('>>> lagrangian')
        n = X.rows
        b = X.cols
        assert(y.size() == n)

        tX = pdf.Matrix(X)
        tX.transpose()
        tXX = tX * X

        rI = bridge.identity_matrix(b)
        rI.set(b -1, b -1, 0.0)
        rI *= self._alpha

        # 一般化逆行列を使う方法
        M = (tXX + rI)
        Minv = M.inverse();

        MinvX = Minv * tX
        self._coef = MinvX * y

    def fit_atomic_charges(self, X, y):
        '''
        X: design matrix (n x n)
        電荷計算用に特化したもの(最後のモデル変数が総電荷用lambda)
        '''
        print('>>> fit ridge atomic charges')
        n = X.rows
        assert(n == X.cols)
        assert(y.size() == n)

        rI = bridge.identity_matrix(n)
        rI.set(n -1, n -1, 0.0)
        rI *= self._alpha

        M = X + rI
        Minv = M.inverse()
        self._coef = Minv * y


    @property
    def coef(self):
        return self._coef


def main():
    # parse args
    parser = argparse.ArgumentParser(description='calculate ridge partial charges')

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

    parser.add_argument("--alpha",
                        type=float,
                        nargs=1,
                        action='store',
                        default=['1.0'],
                        help="ridge alpha")

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
                        default=['ridge_model.vtr'],
                        help="output model vector path (default: ridge_model.vtr).")

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
    tol = 1.0e-3
    alpha = float(args.alpha[0])
    if verbose:
        print("alpha={}".format(alpha))
    esp_data_path = args.espdat[0]

    # make atomlist
    atomlist = []
    if args.db:
        atomlist, charge = pdf.PopUtils.get_atomlist_by_db(args.db)
    elif args.param:
        atomlist, charge = pdf.PopUtils.get_atomlist_by_param(args.param)
    num_of_real_atoms = len(atomlist)

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

    # set constrained value (charge)
    y[num_of_real_atoms] = charge

    # solve ridge
    ridge = Ridge(alpha=alpha)
    ridge.fit_atomic_charges(X, y)
    atom_charges = ridge.coef[:-1]
    pdf.PopUtils.set_atomlist(atom_charges, atomlist)

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
    if os.path.exists("config.ini"):
        logging.config.fileConfig("config.ini",
                                  disable_existing_loggers=False)
    main()

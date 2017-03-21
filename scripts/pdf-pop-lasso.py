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
    
import pdfbridge as bridge
import pdfpytools as pdf

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

   

class Lasso(object):
    def __init__(self, alpha=1.0, max_iter=1000, tol=0.0001):
        self._alpha = alpha
        self._fit_intercept = False
        self._max_iter = max_iter
        self._iter = 0
        
        self._conv_threshold = tol
        self._num_of_converged = 0
        self._is_converged = False
        self._coef = None
        # self._intercept = None

    @property
    def iterations(self):
        return self._iter
        
    @property
    def is_converged(self):
        return self._is_converged
        
    @property
    def coef(self):
        return self._coef
        
    def fit(self, X, y):
        '''
        X: design matrix (n x b)
        '''
        assert(isinstance(X, bridge.Matrix))
        assert(isinstance(y, bridge.Vector))
        if self._fit_intercept:
            pass
        else:
            self._fit_without_intercept(X, y)

    def _fit_without_intercept(self, X, y):
        n = X.rows
        b = X.cols
        assert(y.size() == n)
        
        beta = pdf.Vector(b)
        prev_beta = pdf.Vector(b)
        for self._iter in range(1, self._max_iter +1):
            for j in range(b):
                beta_j = pdf.Vector(beta)
                beta_j[j] = 0.0

                r_j = y - X * beta_j
                X_j = X.get_col_vector(j)
                arg1 = X_j * r_j
                arg2 = self._alpha * n
                
                X_j2 = X_j * X_j
                beta[j] = self._soft_thresholding_operator(arg1, arg2) / X_j2

            if self._is_convergence(self._iter, beta, prev_beta):
                self._is_converged = True
                break

            #print(self._iter)
            #print(beta)
            prev_beta = pdf.Vector(beta)

        self._coef = beta

    def _fit_with_intercept(self, X, y):
        # to inprement
        raise
        
    def _soft_thresholding_operator(self, x, lambda_):
        if x > 0 and lambda_ < abs(x):
            return x - lambda_
        elif x < 0 and lambda_ < abs(x):
            return x + lambda_
        else:
            return 0
        #if abs(x) <= lambda_:
        #    return 0
        #elif x > lambda_:
        #    return x - lambda_
        #else:
        #    return x + lambda_
       
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

    def fit_atomic_charges(self, X, y):
        n = X.rows
        b = X.cols
        assert(y.size() == n)
        assert(n == b)

        beta = pdf.Vector(b)
        prev_beta = pdf.Vector(b)
        for self._iter in range(1, self._max_iter +1):
            delta = pdf.Matrix(n, b)
            for i in range(n -1):
                delta.set(i, i, abs(beta[i]))
            delta *= self._alpha

            M = X + delta
            Minv = M.inverse()

            prev_beta = pdf.Vector(beta)
            beta = Minv * y
            self._coef = beta
                
            if self._is_convergence(self._iter, beta, prev_beta):
                self._is_converged = True
                print("iter={}".format(self._iter))
                break

    def fit_TH(self, X, y):
        print('>>> fit lasso atomic charges')
        n = X.rows
        assert(n == X.cols)
        assert(y.size() == n)

        # initial guess
        coef = bridge.Vector(n)
        prev_coef = bridge.Vector(coef)

        for self._iter in range(1, self._max_iter +1):
            rI = bridge.SymmetricMatrix(n)
            for i in range(n -1):
                rI.set(i, i, abs(coef[i]))
            rI *= self._alpha

            M = X + rI
            Minv = M.inverse()
            coef = Minv * y

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
                        default=['1.0'],
                        help="ridge alpha")
    parser.add_argument("-e", "--espdat",
                        nargs=1,
                        action='store',
                        default=['grid-esp.mpac'],
                        help="ridge alpha")
    parser.add_argument("-s", "--use_scikit_learn",
                        action='store_true',
                        help="use scikit-learn module")
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
    use_scikit_learn = args.use_scikit_learn
    verbose = args.verbose
    max_iter = args.max_iterations[0]
    tol = args.tolerance[0] # 0.0001
    alpha = float(args.alpha[0])
    if verbose:
        print("alpha={}".format(alpha))
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
    if use_scikit_learn:
        print('use scikit-learn module')
        print('alpha={}'.format(alpha))
        import sklearn.linear_model as lm
        lasso = lm.Lasso(alpha=alpha, fit_intercept=False, max_iter=max_iter, tol=tol)
        X.resize(X.rows -1, X.cols -1)
        y.resize(len(y) -1)
        lasso.fit(X.get_ndarray(), y.get_ndarray())
        # print(lasso.intercept_)
        # atom_charges = lasso.coef_[:-1]
        atom_charges = lasso.coef_[:]
        set_atomlist(atom_charges, atomlist)
    else:
        lasso = Lasso(alpha=alpha, max_iter=max_iter, tol=tol)
        lasso.fit_TH(X=X, y=y)
        print(lasso.iterations)
        print(lasso.is_converged)
        atom_charges = lasso.coef[:-1]
        set_atomlist(atom_charges, atomlist)

        # debug
        #X.resize(X.rows -1, X.cols -1)
        #y.resize(len(y) -1)
        #lasso.fit(X=X, y=y) # same as scikit-learn
        #print(lasso.iterations)
        #print(lasso.is_converged)
        #atom_charges = lasso.coef[:]
        #set_atomlist(atom_charges, atomlist)
        
        
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
    

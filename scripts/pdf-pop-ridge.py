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
import copy

import pdfbridge as bridge
import pdfpytools as pdf

import sklearn.linear_model as lm

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
    parser = argparse.ArgumentParser(description='ridge')

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
    parser.add_argument("-e", "--espdat",
                        nargs=1,
                        action='store',
                        default=['grid-esp.mpac'],
                        help="ESP values on grids by msgpack format")
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
    tol = 1.0e-3
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
        ridge = lm.Ridge(alpha=alpha, fit_intercept=False)
        X.resize(X.rows -1, X.cols -1)
        y.resize(len(y) -1)
        ridge.fit(X.get_ndarray(), y.get_ndarray())
        # print(ridge.coef_)
        # print('charge={:.6f}'.format(get_charge(ridge.coef_)))
        # atom_charges = ridge.coef_[:-1]
        atom_charges = ridge.coef_[:]
        set_atomlist(atom_charges, atomlist)
    else:    
        ridge = Ridge(alpha=alpha)

        ridge.fit_atomic_charges(X, y)
        atom_charges = ridge.coef[:-1]
        set_atomlist(atom_charges, atomlist)
        
        # debug
        #X.resize(X.rows -1, X.cols -1)
        #y.resize(len(y) -1)
        #ridge.fit(X, y) # same as scikit-learn
        #atom_charges = ridge.coef[:]
        #set_atomlist(atom_charges, atomlist)
        
        #ridge.fit_atomic_charges(X, y)
        #atom_charges = ridge.coef[:-1]
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
    

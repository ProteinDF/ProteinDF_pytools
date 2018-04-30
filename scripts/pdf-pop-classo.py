#!/usr/bin/env python
# -*- coding: utf-8 -*-

import argparse
import os

import logging
import logging.config
logger = logging.getLogger(__name__)

import sklearn.linear_model as lm

import proteindf_bridge as bridge
import proteindf_tools as pdf


TOO_SMALL = 1.0E-5


class ConstrainedLasso(object):
    """
    [LASSO]: arg min(beta) 1/2 | Y - X * beta |_2^2 + lambda * |beta|_1
    X : p  x p
    Y : p
    beta: p

    [Constrained-LASSO]
    C : m x p
    b : m

    [solve]
    A  : n (n <= p)
    XA : p x n
    XA~: p x (p-n)
    CA : m x n
    CA~: m x (p-n)

    CA^-1: n x m
    Y* := Y - XA * CA^-1 * b  [eq. 10]
       :  p - (p x n)*(n x m)*m
       = p
    X* := XA~ - XA * CA^-1 * CA~  [eq. 10]
       :  (p x (p-n)) - (p x n)*(n x m)*(m x (p-n))
       = p x (p-n)
    theta_A: (p-n) [eq. 10]

    X^-1: p x p
    X~:= X*T X~ = I
      :  (p x (p-n))T * X~ = I
      :  ((p-n) x p) * X~ = I
      = p x (p-n)
    theta~As := CA^-1 * (b - CA~ * theta_As)
             :  (n x m) * (m - (m x (p-n)) * (p-n))   (∵theta_As = theta_A)
             :  n
    Y~:= Y* + lambda * X~ * (CA^-1 * CA~)^T * s  [eq.12]
      :  p  +          (p x (p-n))*((n x m)*(m x (p-n)))^T * n
      :  p  +          (p x (p-n))*(n x (p-n))^T * n
      :  p  +          (p x (p-n))*((p-n) x n) * n


    theta_bar: CA^-1 * (b - CAtilda * theta)
             : (n x m) * (m - (m x (p-n)) * (p-n))

    """

    def __init__(self):
        '''
        alpha: Constant that multiplies the L1 term.
        '''
        # self.alpha = float(alpha)
        #self.lambda_min = float(lambda_min)
        #self.lambda_max = float(lambda_max)
        self.max_itr = 100
        self.step_size = 1.0

        self.coef = [] # answer


    def fit(self, X, y, C, b):
        '''
        X: design matrix (n x b)
        '''
        assert(isinstance(X, bridge.Matrix))
        assert(isinstance(y, bridge.Vector))
        assert(isinstance(C, bridge.Matrix))
        assert(isinstance(b, bridge.Vector))
        self.X = X
        self.y = y
        self.C = C
        self.b = b

        self._is_converged = False
        #self._beta = []
        self._A = [None]
        self._s = [None]
        #self._theta = bridge.Vector()
        #self._theta_bar = bridge.Vector()

        self.step1()
        step_size_was_too_large = False
        while ((self._is_converged == False) and (self.itr < self.max_itr)):
            if step_size_was_too_large:
                self.step5()
            else:
                self.step2()
            logging.debug("itr: {}".format(self.itr))

            self.step3()
            step_size_is_too_large = self.step4()


        # solved!
        answer = bridge.Vector(len(self._theta_bar) + len(self._theta))
        for i in range(len(self._theta_bar)):
            answer[i] = self._theta_bar[i]
        for i in range(len(self._theta)):
            answer[i + len(self._theta_bar)] = self._theta[i]
        self.coef = answer


    def is_converged(self):
        return self._is_converged


    def _do_lasso(self, X, y, alpha):
        assert(isinstance(X, bridge.Matrix))
        assert(isinstance(y, bridge.Vector))
        assert(isinstance(alpha, float))

        logger.debug("begin LASSO...")
        max_iter=100000 # default: 1000
        tol=0.001 # default: 0.0001
        lasso = lm.Lasso(alpha=alpha, fit_intercept=False, max_iter=max_iter, tol=tol)

        lasso.fit(X.get_ndarray(), y.get_ndarray())
        beta = lasso.coef_[:]
        logger.debug("end LASSO...")
        logger.debug("LASSO results: itr={}\n{}".format(lasso.n_iter_, beta))
        beta = pdf.Vector(beta)
        return beta

    def _find_abs_max_elements(self, beta):
        max_val = abs(beta[0])
        for i in range(len(beta)):
            max_val = max(max_val, abs(beta[i]))

        A = []
        s = []
        for i in range(len(beta)):
            if abs(max_val - abs(beta[i])) < TOO_SMALL:
                A.append(i)
                if beta[i] > 0:
                    s.append(1)
                else:
                    s.append(-1)

        return (A, s)


    def _select_by_A(self, in_mat, A):
        '''
        see sec.3
        '''
        assert(isinstance(in_mat, bridge.Matrix))
        assert(isinstance(A, (bridge.Vector, list)))

        m = len(A)
        coef_mat = bridge.Matrix(in_mat.cols, m)
        for i in range(m):
            row = A[i]
            coef_mat.set(row, i, 1.0)

        out_mat = in_mat * coef_mat
        return out_mat


    def _select_by_Abar(self, in_mat, A):
        assert(isinstance(in_mat, bridge.Matrix))
        assert(isinstance(A, (bridge.Vector, list)))

        cols = in_mat.cols
        m = len(A)

        Abar = pdf.Vector(cols - m)
        Abar_index = 0
        for i in range(cols):
            if i not in A:
                Abar[Abar_index] = i
                Abar_index += 1
        return self._select_by_A(in_mat, Abar)


    def _get_Xstar(self, A):
        '''
        eq. 10
        '''
        X = self.X
        C = self.C

        XA = self._select_by_A(X, A)
        XAtilda = self._select_by_Abar(X, A)

        CA = self._select_by_A(C, A)
        logger.debug("X*: CA: {}".format(str(CA)))
        CAinv = CA.inverse()

        CAtilda = self._select_by_Abar(C, A)

        Xstar = XAtilda - XA * CAinv * CAtilda
        return Xstar


    def _get_Ystar(self, A):
        '''
        calculate Y*
        see eq. 10
        '''
        assert(isinstance(A, (bridge.Vector, list)))

        X = self.X
        y = self.y
        C = self.C
        b = self.b

        XA = self._select_by_A(X, A)
        assert(XA.rows == X.rows)

        CA = self._select_by_A(C, A)
        assert(CA.rows == C.rows)

        # CAinv means CA^-1
        CAinv = CA.inverse()

        # Y*
        Ystar = y - XA * CAinv * b

        return Ystar


    def _get_Xtilda(self, A):
        '''
        calculate X~
        X^*T * X~ = I
        '''
        assert(isinstance(A, (pdf.Vector, list)))

        X = self.X
        Xstar = self._get_Xstar(A)

        # transpose
        XstarT = bridge.Matrix(Xstar)
        XstarT.transpose()
        assert(XstarT.rows == (X.rows - len(A)))
        assert(XstarT.cols == X.rows)

        Xtilda = XstarT.pseudo_inverse()
        assert(Xtilda.rows == X.rows)
        assert(Xtilda.cols == (X.rows - len(A)))

        return Xtilda


    def _get_Ytilda(self, lambda_, A, s):
        '''
        calculate Y~
        see eq.12
        '''
        assert(isinstance(lambda_, float))
        assert(isinstance(A, (bridge.Vector, list)))
        assert(isinstance(s, (bridge.Vector, list)))
        assert(len(A) == len(s))

        X = self.X
        C = self.C

        Xinv = X.inverse()
        Ystar = self._get_Ystar(A)

        CA = self._select_by_A(C, A)
        CAinv = CA.inverse()
        CAbar = self._select_by_Abar(C, A)

        Xtilda = self._get_Xtilda(A)

        CC = CAinv * CAbar
        # transpose
        CCt = bridge.Matrix(CC)
        CCt.transpose()

        itr = self.itr
        lambda_val = self._lambda[itr]

        # Ytilda = Ystar + lambda_val * Xtilda * CCt * s
        mat_s = bridge.Matrix(len(s), 1)
        for i in range(len(s)):
            mat_s.set(i, 0, s[i])
        term = lambda_val * Xtilda * CCt * mat_s

        term_vtr = bridge.Vector(term.rows)
        for i in range(term.rows):
            term_vtr[i] = term.get(i, 0)

        Ytilda = Ystar + term_vtr
        return Ytilda


    def _get_theta_bar(self, A, theta):
        assert(isinstance(theta, pdf.Vector))

        b = self.b
        C = self.C

        CA = self._select_by_A(C, A)
        CAinv = CA.inverse()

        CAbar = self._select_by_Abar(C, A)

        theta_bar = CAinv * (b - CAbar * theta)
        assert(len(theta_bar) == len(A))

        return theta_bar


    def _get_s(self, a):
        assert(isinstance(a, bridge.Vector))

        d = len(a)
        s = pdf.Vector(d)
        for i in range(d):
            if a[i] >= 0.0:
                s[i] = +1.0
            else:
                s[i] = -1.0

        return s

    # ------------------------------------------------------------------
    # solve
    # ------------------------------------------------------------------
    def step1(self):
        # initialize beta_0 by solving eq.9 using lambda_0 = lambda_max
        self.itr = 0
        self._lambda = [self.lambda_max]
        beta = self._do_lasso(self.X, self.y, alpha=self._lambda[0])
        logger.debug("step1: beta={}".format(beta))
        self._beta = [beta]


    def step2(self):
        '''
        k: the number of step
        '''
        self.itr += 1

        itr = self.itr
        logger.debug("step2: itr={}".format(itr))
        logger.debug("step2: len(beta)={}".format(len(self._beta)))
        (A, s) = self._find_abs_max_elements(self._beta[itr -1])
        logging.debug("step2: A=\n{}".format(str(A)))
        logging.debug("step2: s=\n{}".format(str(s)))
        self._A.append(A)
        self._s.append(s)
        assert(len(self._A) == self.itr +1)
        assert(len(self._s) == self.itr +1)

        lambda_val = pow(10.0, - self.step_size) * self._lambda[itr -1]
        self._lambda.append(lambda_val)
        assert(len(self._lambda) == self.itr +1)
        logger.debug("new lambda: {} at {}-th".format(lambda_val, self.itr))

        if self._lambda[itr] < self.lambda_min:
            self._is_converged = True


    def step3(self):
        '''
        '''
        itr = self.itr
        # eq.12
        Xstar = self._get_Xstar(self._A[itr])
        lambda_val = self._lambda[itr]
        Ytilda = self._get_Ytilda(lambda_val, self._A[itr], self._s[itr])

        self._theta = self._do_lasso(Xstar, Ytilda, alpha=self._lambda[itr])
        self._theta_bar = self._get_theta_bar(self._A[itr], self._theta)


    def step4(self):
        itr = self.itr
        logger.debug("step4: len(beta)=\n{}".format(len(self._beta)))

        s = self._get_s(self._theta_bar)
        step_size_was_too_large = (self._s[itr] != s)
        if step_size_was_too_large == False:
            # update beta
            beta = pdf.Vector(len(self._theta_bar) + len(self._theta))
            for i in range(len(self._theta_bar)):
                beta[i] = theta_bar[i]
            for i in range(len(self._theta)):
                beta[len(theta_bar) + i] = theta[i]
            self._beta.append(beta)
        else:
            self._beta.append(self._beta[itr -1])

        logger.debug("step4: len(beta)=\n{}".format(len(self._beta)))
        return step_size_was_too_large


    def step5(self):
        itr = self.itr

        lambda_k = self._lambda[itr]
        lambda_k1 = self._lambda_[itr -1]
        self._lambda[itr] = lambda_k1 -0.5*(lambda_k1 - lambda_k)


def main():
    # parse args
    parser = argparse.ArgumentParser(description='calculate Constrained LASSO partial atomic charges')

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

    parser.add_argument("--max-alpha",
                        type=float,
                        nargs=1,
                        action='store',
                        default=['1.0'],
                        help="max alpha value")
    parser.add_argument("--min-alpha",
                        type=float,
                        nargs=1,
                        action='store',
                        default=['1.0'],
                        help="min alpha value")

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
                        default=['classo_model.vtr'],
                        help="output model vector path (default: classo_model.vtr).")

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
    max_alpha = float(args.max_alpha[0])
    min_alpha = float(args.min_alpha[0])
    if verbose:
        print("max alpha={}".format(max_alpha))
        print("min alpha={}".format(min_alpha))
    esp_data_path = args.espdat[0]

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
    X.resize(num_of_real_atoms, num_of_real_atoms)

    y = pdf.Vector()
    y.load(target_vector_path)
    y.resize(num_of_real_atoms)

    # make Constrained Matrix C
    C = bridge.Matrix(1, num_of_real_atoms)
    for i in range(num_of_real_atoms):
        C.set(0, i, 1.0)

    # make b
    b = bridge.Vector(1)
    b.set(0, charge)

    # solve
    classo = ConstrainedLasso()
    classo.step_size = 1.0
    classo.lambda_max = max_alpha
    classo.lambda_min = min_alpha

    classo.fit(X, y, C, b)
    logger.info("itr: {}".format(classo.itr))
    atom_charges = classo.coef
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

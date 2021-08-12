#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import math
from types import *
import argparse

from matplotlib import use
import numpy 

import proteindf_bridge as bridge
import proteindf_tools as pdf

import logging

def dosGaussian(eps, omega = 0.005, eLevels=[]):
    answer = 0.0
    coef1 = omega * math.sqrt(0.5 * math.pi)
    coef2 = -2.0 / (omega * omega)

    for el in eLevels:
        x = eps - el
        x2 = x * x
        answer += math.exp(coef2 * x2)

    answer *= coef1
    return answer

def dosLorentzian(eps, omega = 0.005, eLevels=[]):
    answer = 0.0
    omega2 = omega * omega

    for el in eLevels:
        x = eps - el
        x2 = x * x
        answer += 1.0 / (omega2 + 4.0 * x2)

    answer *= (2.0 * omega)/(math.pi)
    return answer


def main():
    # parse args
    parser = argparse.ArgumentParser(description='plot density of state')
    parser.add_argument('-d', '--db',
                       nargs=1,
                       default=['pdfresults.h5'],
                       help='ProteinDF results file')
    parser.add_argument('-r', '--resolution',
                       nargs=1,
                       default=[0.005],
                       help='resolution (default: 0.005 eV)')

    group_func = parser.add_mutually_exclusive_group()
    group_func.add_argument('--lorentzian',
                        action="store_true",
                        default=False,
                        help='use Lorentzian')
    group_func.add_argument('--gaussian',
                        action="store_true",
                        default=False,
                        help='use Gaussian')

    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=['elevel.png'],
                        help='output graph path')
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    parser.add_argument('-D', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()
    # print(args)

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    output_path = args.output[0]
    pdfresults_db = args.db[0]

    use_func = "gaussian"
    if args.lorentzian:
        use_func = "lorentzian"

    AU2EV = 27.2116

    # load DB
    pdfparam = pdf.PdfParam_H5()
    pdfparam.open(pdfresults_db)

    # eLevels
    method = pdfparam.method
    run_type = "rks"
    HOMO_level = pdfparam.get_HOMO_level(run_type)

    itr = pdfparam.iterations
    eigvals = pdfparam.get_energy_level(method, itr)

    # make data
    occLevels = [None] * HOMO_level
    for i in range(HOMO_level):
        occLevels[i] = eigvals[i] * AU2EV

    minLevel = -20.0
    maxLevel = 5.0
    omega = 0.005
    steps = int((maxLevel - minLevel) / omega)
    data_path = os.path.join(".", "dos.dat")
    with open(data_path, 'w') as dat:
        for step in range(steps +1):
            e = minLevel + omega * step
            if use_func == "lorentzian":
                intensity = dosLorentzian(e, omega, occLevels)
            else:
                intensity = dosGaussian(e, omega, occLevels)
            dat.write('% 16.10f, % 16.10f\n' % (e, intensity))

    # plot data
    graph = pdf.DfGraph2D()
    graph.load_data("dos.dat")
    graph.title = "Density of State"
    graph.xmax = maxLevel
    graph.xmin = minLevel
    graph.xlabel = "energy level / eV"
    graph.ylabel = "Intensity"
    graph.save("dos.png")


if __name__ == '__main__':
    main()

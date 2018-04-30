#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from jinja2 import Environment, FileSystemLoader

import pdfbridge as bridge
import pdfpytools as pdf


pdf_templ_str = """
>>>>MAIN
        step-control    = [create integral guess scf]
        cut-value       = 1.0e-10
        charge-extrapolate-number = 0
        scf-start-guess = harris
        orbital_independence_threshold = 0.00
        orbital_independence_threshold/canonical = 0.00
        orbital_independence_threshold/lowdin = 0.00
        orbital-overlap-correspondence  = off
        orbital-overlap-correspondence-first    = off
        max-iteration   = 100
        method  = rks
        method/rks/electron-number      = {{ "%d"|format(rks_electrons) }}
        method/rks/occlevel     = [ 1 - {{ "%d"|format(rks_electrons / 2) }} ]
        convergence/type        = density
        convergence/threshold   = 1e-4
        convergence/threshold-energy    = 1e-5
        scf-acceleration        = damping
        scf-acceleration/damping/damping-factor = 0.67
        scf-acceleration/damping/damping-type = density_matrix
        xc-potential    = HF
        J_engine = CD
        K_engine = CD
        XC_engine = grid
        CDAM_tau = 1.0E-10
        CD_epsilon = 1.0E-4

>>>>MOLECULE
        geometry/cartesian/unit = angstrom
        geometry/cartesian/input        = {
        {%- for atom in atomlist %}
            {{ atom.symbol }}    {{ "% 8.5f"|format(atom.x) }} {{ "% 8.5f"|format(atom.y) }} {{ "% 8.5f"|format(atom.z) }}
        {%- endfor %}
        }end

        basis-set/orbital       = {
                H = "O-DZVP2.H"
                O = "O-DZVP2.O"
                C = "O-DZVP2.C"
                N = "O-DZVP2.N"
        }end

        basis-set/density-auxiliary     = {
                H = "A-DZVP2.H"
                O = "A-DZVP2.O"
                C = "A-DZVP2.C"
                N = "A-DZVP2.N"
        }end

        basis-set/exchange-auxiliary    = {
                H = "A-DZVP2.H"
                O = "A-DZVP2.O"
                C = "A-DZVP2.C"
                N = "A-DZVP2.N"
        }end

"""


def main():
    # parse args
    parser = argparse.ArgumentParser(description='make ProteinDF input from xyz file.')
    parser.add_argument("-t", "--template",
                        action="store",
                        default=[""],
                        help="ProteinDF input template (Jinja2 format)")
    parser.add_argument("-o", "--output",
                        action="store",
                        default=["fl_Userinput"],
                        help="output path")
    parser.add_argument("-v", "--verbose",
                        action="store_true",
                        default=False)
    parser.add_argument('-D', '--debug',
                        action='store_true',
                        default=False)
    parser.add_argument("xyz",
                        nargs=1,
                        action='store',
                        help='XYZ file')
    args = parser.parse_args()
        
    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    pdf_templ = args.template[0]
    pdf_output = args.output[0]
    xyz_path = args.xyz[0]
    
    if verbose:
        sys.stderr.write("template: {}\n".format(pdf_templ))
        sys.stderr.write("xyz file: {}\n".format(xyz_path))
        sys.stderr.write("output: {}\n".format(pdf_output))
        
    xyz = bridge.Xyz()
    xyz.load(xyz_path)
    atomgroup = xyz.get_atom_group()

    pdfparam = pdf.PdfParam()
    pdfparam.molecule = atomgroup
    rks_electrons = pdfparam.num_of_electrons
    atomlist = []
    for i in atomgroup.get_atom_list():
        atom = {}
        atom['symbol'] = i.symbol
        atom['x'] = i.xyz.x
        atom['y'] = i.xyz.y
        atom['z'] = i.xyz.z
        atomlist.append(atom)

    # make input
    env = Environment(loader=FileSystemLoader('./', encoding='utf8'))
    tpl = env.from_string(pdf_templ_str)
    if len(pdf_templ) > 0:
        tpl = env.get_template(pdf_templ)
    

    pdf_input = tpl.render(rks_electrons=rks_electrons,
                           atomlist=atomlist)

    # output
    with open(pdf_output, "w") as f:
        f.write(pdf_input)

        
if __name__ == '__main__':
    main()
        

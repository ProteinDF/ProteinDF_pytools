#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import argparse
from jinja2 import Environment, FileSystemLoader

import pdfbridge as bridge
import pdfpytools as pdf


templ_str = """
#p {{ method }}/{{ basisset }} sp 5D GFprint GFinput NoSymm
SCF(tight, maxcycle=500)

{{ title }} {{ method }}/{{ basisset }}

0 1
{%- for atom in atomlist %}
{{ atom.symbol }}    {{ "% 8.5f"|format(atom.x) }} {{ "% 8.5f"|format(atom.y) }} {{ "% 8.5f"|format(atom.z) }}
{%- endfor %}
"""


def main():
    # parse args
    parser = argparse.ArgumentParser(description='make Gaussian input from xyz file.')
    parser.add_argument("-t", "--template",
                        action="store",
                        default=[""],
                        help="input template (Jinja2 format)")
    parser.add_argument("-o", "--output",
                        action="store",
                        default=["sample.gau"],
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

    templ_path = args.template[0]
    output_path = args.output[0]
    xyz_path = args.xyz[0]
    
    if verbose:
        sys.stderr.write("template: {}\n".format(templ_path))
        sys.stderr.write("xyz file: {}\n".format(xyz_path))
        sys.stderr.write("output: {}\n".format(output_path))
        
    xyz = bridge.Xyz()
    xyz.load(xyz_path)
    atomgroup = xyz.get_atom_group()

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
    tpl = env.from_string(g09_templ_str)
    if len(templ_path) > 0:
        tpl = env.get_template(templ_path)
    

    input_content = tpl.render(atomlist=atomlist)

    
    # output
    with open(output_path, "w") as f:
        f.write(input_content)

        
if __name__ == '__main__':
    main()
        

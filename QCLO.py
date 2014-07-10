#!/usr/bin/env python
# -*- coding: utf-8 -*-

# Copyright (C) 2002-2014 The ProteinDF project
# see also AUTHORS and README.
#
# This file is part of ProteinDF.
#
# ProteinDF is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# ProteinDF is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with ProteinDF.  If not, see <http://www.gnu.org/licenses/>.

import os
import sys
import argparse
import logging
import logging.config
import math
try:
    import msgpack
except:
    import msgpack_pure as msgpack

import bridge
import pdf


class QcAtom(bridge.Atom):
    def __init__(self, *args, **kwargs):
        super(QcAtom, self).__init__(*args, **kwargs)
        suffix = "." + self.symbol
        self.basisset = kwargs.get('basisset', 'O-DZVP2' + suffix)
        self.basisset_j = kwargs.get('basisset_j', 'A-DZVP2' + suffix)
        self.basisset_xc = kwargs.get('basisset_xc', 'A-DZVP2' + suffix)
        self.basisset_gridfree = kwargs.get('basisset_gridfree', 'O-DZVP2' + suffix)

    #
    def get_number_of_AOs(self):
        basis2 = pdf.Basis2()
        bs = basis2.get_basisset(self.basisset)
        return len(bs)
        
    # basisset -------------------------------------------------------
    def _get_basisset(self):
        return self._basisset

    def _set_basisset(self, name):
        self._basisset = bridge.Utils.byte2str(name)

    basisset = property(_get_basisset, _set_basisset)

    # basisset_j ------------------------------------------------------
    def _get_basisset_j(self):
        return self._basisset_j

    def _set_basisset_j(self, name):
        self._basisset_j = bridge.Utils.byte2str(name)

    basisset_j = property(_get_basisset_j, _set_basisset_j)

    # basisset_xc -----------------------------------------------------
    def _get_basisset_xc(self):
        return self._basisset_xc

    def _set_basisset_xc(self, name):
        self._basisset_xc = bridge.Utils.byte2str(name)

    basisset_xc = property(_get_basisset_xc, _set_basisset_xc)

    # basisset_gridfree -----------------------------------------------
    def _get_basisset_gridfree(self):
        return self._basisset_gridfree

    def _set_basisset_gridfree(self, name):
        self._basisset_gridfree = bridge.Utils.byte2str(name)

    basisset_gridfree = property(_get_basisset_gridfree, _set_basisset_gridfree)

    # atom label ------------------------------------------------------
    def _get_atomlabel(self):
        symbol = self.symbol
        label = self.label
        kind = symbol
        if len(label) > 0:
            kind = '{}@{}'.format(symbol, label)
        return kind

    atomlabel = property(_get_atomlabel)


class QcFragment(bridge.AtomGroup):
    def __init__(self, *args, **kwargs):
        super(QcFragment, self).__init__(*args, **kwargs)
        self._lo = None

        if len(args) > 0:
            if len(args) == 1:
                rhs = args[0]
                if isinstance(rhs, QcFragment):
                    self._lo = rhs._lo

    #
    def get_number_of_AOs(self):
        def get_number_of_AOs_sub(ag):
            AOs = 0
            for key, subgrp in ag.groups():
                AOs += get_number_of_AOs_sub(subgrp)
            for key, atm in ag.atoms():
                AOs += atm.get_number_of_AOs()
            return AOs

        return get_number_of_AOs_sub(self)

    # basisset ---------------------------------------------------------
    def set_basisset(self, pdfparam):
        for key, frg in self.groups():
            frg.set_basisset(pdfparam)
        for key, atm in self.atoms():
            symbol = atm.symbol
            if symbol != 'X':
                atomlabel = atm.atomlabel
                bsname = atm.basisset
                bsname_j = atm.basisset_j
                bsname_xc = atm.basisset_xc
                bsname_gridfree = atm.basisset_gridfree
                pdfparam.set_basisset_name(atomlabel, bsname)
                pdfparam.set_basisset_j_name(atomlabel, bsname_j)
                pdfparam.set_basisset_xc_name(atomlabel, bsname_xc)
                pdfparam.set_basisset_gridfree_name(atomlabel, bsname_gridfree)

    # ------------------------------------------------------------------
    def set_atom(self, key, value):
        key = bridge.Utils.byte2str(key)
        self._atoms[key] = QcAtom(value, parent = self)
        
    # ------------------------------------------------------------------
    def __str__(self):
        answer = super(QcFragment, self).__str__()
        return answer 
    

class QcFrame(object):
    _CALCD_SP = 1
    _CALCD_LO = 2

    def __init__(self, name):
        nullHandler = bridge.NullHandler()
        self._logger = logging.getLogger(__name__)
        self._logger.addHandler(nullHandler)

        self._name = str(name)
        self._fragments = {}
        self._state = 0
        self._frame_molecule = None
        
        self.make_workdir()

    # ------------------------------------------------------------------
    def _get_frame_molecule(self):
        if self._frame_molecule == None:
            frame_molecule = bridge.AtomGroup()
            for frg_name, frg in self._fragments.items():
                frame_molecule[frg_name] = bridge.AtomGroup(frg)
            self._frame_molecule = frame_molecule
        return self._frame_molecule

    frame_molecule = property(_get_frame_molecule)
    # ------------------------------------------------------------------
    def make_workdir(self):
        workdir = self.name
        self._logger.debug('make workdir: {}'.format(workdir))
        if not os.path.exists(workdir):
            os.mkdir(workdir)
        else:
            self._logger.warning('already exit: {}'.format(workdir))

    #def is_neutralize(self):
    #    charge = self.frame_molecule.charge
    #    return (math.fabs(charge) < 1.0E-5)

    def neutralize(self):
        pass
        
    def calc_sp(self):
        pdfsim = pdf.PdfSim()
        pdfparam = pdf.get_default_pdfparam()

        for frg_name, frg in self._fragments.items():
            frg.set_basisset(pdfparam)
        pdfparam.molecule = self.frame_molecule

        # for debug
        #self.output_xyz("{}/model.xyz".format(self.name))

        #if not self.is_neutralize():
        #    self._logger.warning('The frame-molecule is NOT neutralized.')
        # 

        pdfsim.sp(pdfparam,
                  workdir = self.name,
                  dry_run = True)


    def calc_lo(self):
        for frg_name, frg in self._fragments.items():
            frg_AOs = frg.get_number_of_AOs()
            print('frg={}, AOs={}'.format(frg_name, frg_AOs))

    def pickup_lo(self):
        pass

    # name -------------------------------------------------------------
    def _get_name(self):
        return self._name

    name = property(_get_name)

    # basisset ---------------------------------------------------------
    def _set_basisset(self, pdfparam):
        for fragment_name, fragment in self._fragments.items():
            fragment.set_basisset(pdfparam)
    
    # outout XYZ -------------------------------------------------------
    def output_xyz(self, file_path):
        xyz = bridge.Xyz(self.frame_molecule)
        xyz.save(file_path)

    # operator[] -------------------------------------------------------
    def __getitem__(self, fagment_name):
        fragment_name = str(fragment_name)
        return self._fragments.getdefault(fragment_name, None)

    def __setitem__(self, fragment_name, fragment):
        self._frame_molecule = None # clear molecule cache
        fragment_name = str(fragment_name)
        fragment = QcFragment(fragment)
        self._fragments[fragment_name] = fragment

    # operator str -----------------------------------------------------
    def __str__(self):
        answer = ""
        for key, fragment in self._fragments.items():
            answer += '>> fragment: {}\n'.format(key)
            answer += str(fragment)
            answer += '\n'
        return answer


def main():
    parser = argparse.ArgumentParser(description='QCLO program')
    parser.add_argument('brdfile',
                        nargs=1,
                        help='Bridge file')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        default=False)
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        default=False)
    args = parser.parse_args()

    # setting
    verbose = args.verbose
    if args.debug:
        logging.basicConfig(level=logging.DEBUG)
    brdfile_path = args.brdfile[0]

    # load
    if verbose:
        print("reading: {}".format(brdfile_path))
    brdfile = open(brdfile_path, 'rb')
    brddata = msgpack.unpackb(brdfile.read())
    brdfile.close()

    # all models
    models = bridge.AtomGroup(brddata)
    #print('>>>> models')
    #print(models)
    #print('----')

    # model
    model = models["model_1"]
    # print(model)

    #
    N_TERM = 0 # N-termがNMEなら1
    C_TERM = 0 # C-termがACEなら1
    modeling = bridge.Modeling()
    for chain_name, chain in model.groups():
        max_resid = chain.get_number_of_groups()
        logging.info("chain {}: max_resid={}".format(chain_name, max_resid))
        #logging.info(str(chain))
        
        for resid, res in chain.groups():
            resid = int(resid)
            print(">>>> resid: {}".format(resid))
            res_model = bridge.AtomGroup(res)

            logging.debug("resid = {}".format(resid))

            # create base fragment
            fragment_res = QcFragment(res_model)

            ACE = QcFragment()
            if resid > 1 + N_TERM:
                #logging.debug("add ACE")

                prev = chain[resid -1]
                #logging.debug(">>> prev")
                #logging.debug("\n" + str(prev))
                # ag_ACE = modeling.get_ACE(res_model, prev)
                ag_ACE = modeling.get_ACE_simple(prev)
                ACE = QcFragment(ag_ACE)

            NME = QcFragment()
            if resid + C_TERM < max_resid:
                #logging.debug("add NME")

                ahead = chain[resid +1]
                #logging.debug(">>> ahead")
                #logging.debug("\n" + str(ahead))
                # ag_NME = modeling.get_NME(res_model, ahead)
                ag_NME = modeling.get_NME_simple(ahead)
                NME = QcFragment(ag_NME)

            frame = QcFrame("res_{}-{}".format(resid, resid))
            frame[resid] = fragment_res
            frame["ACE"] = ACE
            frame["NME"] = NME
            
            #print(">>>> resid: {}".format(resid))
            #print(frame)
            #ag = frame.atomgroup
            #print(str(ag))

            frame.calc_sp()
            frame.calc_lo()

            #break

    exit()

    #======
    for chain in protein:
        chain

    #=======
    frame1_3 = QcFrame()
    frame2_4 = QcFrame()
    frame3_5 = QcFrame()
    frame1_3["body"] = chain.select_seq(1, 3)
    frame2_4["body"] = chain.select_seq(2, 4)
    frame3_5["body"] = chain.select_seq(3, 5)

    frame1_3.add_ACE() # add ACE submolecule and fragment "ACE"
    frame1_3.add_NME() # add "NME"
    for fragment in frame1_3.fragments():
        print(fragment)

    frame2_4.add_ACE() # add ACE submolecule and fragment "ACE"
    frame2_4.add_NME() # add "NME"
    frame3_5.add_ACE() # add ACE submolecule and fragment "ACE"
    frame3_5.add_NME() # add "NME"

    frame1_3.calc_sp()
    frame2_4.calc_sp()
    frame3_5.calc_sp()

    frame1_5 = QcFrame()
    frame1_5["ACE"] = frame1_3["ACE"] # .get_fragment("ACE")
    frame1_5["body"] = frame1_3["body"]
    frame1_5["body"] += frame2_4["body"]
    frame1_5["body"] += frame3_5["body"]
    frame1_5["NME"] = frame1_3["NME"]

    frame1_5.make_guess(method="QCLO")
    frame1_5.calc_sp()

    #=====


if __name__ == '__main__':
    logging.config.fileConfig('logconfig.ini')
    main()

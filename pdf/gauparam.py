#!/usr/bin/env python
# -*- coding: utf-8 -*-

import pdf

class GaussianParam(pdf.QmSim):
    """
    Gaussianの計算条件を設定する
    """
    def __init__(self):
        self._data = {}
        self._init_gen_data()

    def _init_gen_data(self):
        self._gen_data = {}

        # DZ
        self._gen_data['DZ'] = {}
        self._gen_data['DZ']['H']= """
 H 0
 S 3
      1.92406000E+01     3.28280000E-02
      2.89920000E+00     2.31208000E-01
      6.53400000E-01     8.17238000E-01
 S 1
      1.77600000E-01     1.00000000E+00
****"""
        self._gen_data['DZ']['N'] = """
 N 0
 S 6
       5.90944000E+03     2.00400000E-03
       8.87451000E+02     1.53100000E-02
       2.04749000E+02     7.42930000E-02
       5.98376000E+01     2.53364000E-01
       1.99981000E+01     6.00576000E-01
       2.68600000E+00     2.45111000E-01
 S 1
       7.19270000E+00     1.00000000E+00
 S 1
       7.00000000E-01     1.00000000E+00
 S 1
       2.13300000E-01     1.00000000E+00
 P 4
       2.67860000E+01     1.82570000E-02
       5.95640000E+00     1.16407000E-01
       1.70740000E+00     3.90111000E-01
       5.31400000E-01     6.37221000E-01
 P 1
       1.65400000E-01     1.00000000E+00
****"""
        # DZP
        self._gen_data['DZP'] = {}
        self._gen_data['DZP']['H'] = """
 H 0
 S 3
      1.92406000E+01     3.28280000E-02
      2.89920000E+00     2.31208000E-01
      6.53400000E-01     8.17238000E-01
 S 1
      1.77600000E-01     1.00000000E+00
 P 1
      1.00000000E+00     1.00000000E+00
****"""
        
        self._gen_data['DZP']['N'] = """
N 0
 S 6
      5.90944000E+03     2.00400000E-03
      8.87451000E+02     1.53100000E-02
      2.04749000E+02     7.42930000E-02
      5.98376000E+01     2.53364000E-01
      1.99981000E+01     6.00576000E-01
      2.68600000E+00     2.45111000E-01
 S 1
      7.19270000E+00     1.00000000E+00
 S 1
      7.00000000E-01     1.00000000E+00
 S 1
      2.13300000E-01     1.00000000E+00
 P 4
      2.67860000E+01     1.82570000E-02
      5.95640000E+00     1.16407000E-01
      1.70740000E+00     3.90111000E-01
      5.31400000E-01     6.37221000E-01
 P 1
      1.65400000E-01     1.00000000E+00
 D 1
      8.00000000E-01     1.00000000E+00
****"""

        self._gen_data['DZVP'] = {}
        self._gen_data['DZVP']['H'] = """
 H 0
 S 4
      5.09991780E+01     9.66050000E-03
      7.48321810E+00     7.37289000E-02
      1.77746760E+00     2.95858100E-01
      5.19329500E-01     7.15905300E-01
 S 1
      1.54110000E-01     1.00000000E+00
****"""
        self._gen_data['DZVP']['N'] = """
 N 0
 S 6
      3.84541490E+03     2.01860000E-03
      5.77533230E+02     1.54078000E-02
      1.31319830E+02     7.53714000E-02
      3.68237810E+01     2.48212200E-01
      1.16701150E+01     4.79827400E-01
      3.85426040E+00     3.31801200E-01
 S 2
      7.82956110E+00    -7.76669000E-02
      6.87735100E-01     5.65459800E-01
 S 1
      2.04038800E-01     1.00000000E+00
 P 4
      2.68098410E+01     1.54663000E-02
      6.06815400E+00     9.64397000E-02
      1.76762560E+00     3.08361000E-01
      5.46672700E-01     4.91159700E-01
 P 1
      1.58728900E-01     1.00000000E+00
 D 1
      7.00000000E-01     1.00000000E+00
****"""

        self._gen_data['DZVP2'] = {}
        self._gen_data['DZVP2']['H'] = """
 H 0
 S 4
      5.09991780E+01     9.66050000E-03
      7.48321810E+00     7.37289000E-02
      1.77746760E+00     2.95858100E-01
      5.19329500E-01     7.15905300E-01
 S 1
      1.54110000E-01     1.00000000E+00
 P 1
      7.50000000E-01     1.00000000E+00
****"""
        self._gen_data['DZVP2']['N'] = """
 N 0
 S 7
      8.10417610E+03     7.96900000E-04
      1.21731380E+03     6.12890000E-03
      2.77739930E+02     3.10471000E-02
      7.88475980E+01     1.15368200E-01
      2.55371610E+01     3.02573800E-01
      9.00457110E+00     4.55791300E-01
      3.28352780E+00     2.43020800E-01
 S 2
      7.84935730E+00    -7.76364000E-02
      6.86223900E-01     5.67981500E-01
 S 1
      2.03502600E-01     1.00000000E+00
 P 5
      4.90146080E+01    -5.90070000E-03
      1.13166710E+01    -4.16444000E-02
      3.40340530E+00    -1.61024900E-01
      1.16111070E+00    -3.58353800E-01
      3.95335800E-01    -4.48841500E-01
 P 1
      1.26898100E-01     1.00000000E+00
 D 1
      7.00000000E-01     1.00000000E+00
****"""

    # title
    def _get_title(self):
        return self._data.get('title', '')

    def _set_title(self, title):
        self._data['title'] = str(title)

    title = property(_get_title, _set_title)
    
    # jobtype
    def _get_jobtype(self):
        return self._data.get('jobtype', 'sp')

    def _set_jobtype(self, jobtype):
        self._data['jobtype'] = str(jobtype)

    jobtype = property(_get_jobtype, _set_jobtype)

    # basisset
    def _get_basisset(self):
        return self._data.get('basisset', '6-31G')

    def _set_basisset(self, basisset):
        self._data['basisset'] = str(basisset)

    basisset = property(_get_basisset, _set_basisset)

    # ext_basisset for gen
    def _get_ext_basisset(self):
        return self._data.get('ext_basisset', '')

    def _set_ext_basisset(self, basisset):
        self._data['ext_basisset'] = str(basisset)

    ext_basisset = property(_get_ext_basisset, _set_ext_basisset)

    # charge
    def _get_charge(self):
        return self._data.get('charge', 0)

    def _set_charge(self, charge):
        self._data['charge'] = int(charge)

    charge = property(_get_charge, _set_charge)
    
    # 
    def get_inputfile_contents(self):
        data = ''
        data += "#P {method}/{basisset} {jobtype} \n".format(
            method = self.method,
            basisset = self.basisset,
            jobtype = self.jobtype
            )
        data += "5D GFPrint GFInput NoSymm \n"
        data += "Integral(SG1Grid, BWeights, NoJEngine) NoRaff \n"
        data += "\n"
        data += "{title} \n".format(title = self.title)
        data += "\n"
        data += "{charge} 1 \n".format(charge = self.charge)
        data += self._get_inputfile_contents_geometry()
        if self.basisset == 'gen':
            data += self._get_inputfile_contents_basisset_gen()
        
        return data

    def _get_inputfile_contents_geometry(self, atom_group = None):
        output = ""
        if (atom_group == None):
            output += self._get_inputfile_contents_geometry(self.molecule)
        else:
            for key, group in atom_group.groups():
                output += self._get_inputfile_contents_geometry(group)
            for key, atom in atom_group.atoms():
                output += "{} {} {} {}\n".format(
                    atom.symbol,
                    atom.xyz.x, atom.xyz.y, atom.xyz.z)
        return output

    def _get_inputfile_contents_basisset_gen(self):
        data = ""
        basisset_name = self.ext_basisset
        kinds = self.molecule.get_atom_kinds()
        for atom in kinds:
            data += self._gen_data[basisset_name][atom]
        data += "\n\n"
        return data

    

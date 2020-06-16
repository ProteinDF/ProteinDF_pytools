#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import codecs
import argparse
import requests
import re
import pprint
import pickle
import copy
from bs4 import BeautifulSoup
try:
    set
except NameError:
    from sets import Set as set

import proteindf_tools as pdf
import proteindf_bridge as bridge

import logging.config
from logging import getLogger, DEBUG, INFO
logger = getLogger(__name__)


class MakeBasis2(object):
    MAX_SHELL_TYPE = 'spdfg'

    def __init__(self,
                 alias_path=None,
                 cache_path=None,
                 output_path=None,
                 uncontract_shells='spd',
                 force_download=False,
                 verbose=False,
                 debug=False):
        self._uncontract_shells = uncontract_shells
        self._force_download = force_download
        self._debug = debug

        self._load_alias(alias_path)

        self._cache = {}
        self._cache_path = cache_path

        self._calced_set = set()

        self._stdout = sys.stdout
        if output_path != None:
            self._stdout = open(output_path, 'w')
        # self._stderr = open(os.devnull, 'a')
        # if verbose:
        #    self._stderr = sys.stderr

        self._load_cache()

    def _alias2real(self, alias):
        real_name = alias
        if alias in self._alias.keys():
            real_name = self._alias[alias]
        return real_name

    def show_list(self):
        basisset_info = self._get_basissets_info()
        for basisset_name in basisset_info.keys():
            self._stdout.write(basisset_name + "\n")

    def run(self, basisset_names=[], header=""):
        # output header
        self._stdout.write(header + '\n')

        # make basis set
        basisset_info = self._get_basissets_info()

        self._cache.setdefault('Gaussian94', {})
        self._cache.setdefault('ProteinDF', {})
        self._cache.setdefault('ProteinDF-DFT', {})
        self._cache.setdefault('ProteinDF-MRD', {})

        if len(basisset_names) == 0:
            basisset_names = basisset_info.keys()

        num_of_items = len(basisset_names)
        count = 0
        for basisset_name in basisset_names:
            logger.info('### {0:>4}/{1:>4}'.format(count + 1, num_of_items))
            self._run_basisset(basisset_info, basisset_name)
            count += 1
        assert(num_of_items == count)

    def _run_basisset(self, basisset_info, alias):
        basisset_name = self._alias2real(alias)

        logger.debug("run getting basisset: {}".format(basisset_name))
        if basisset_name in basisset_info.keys():
            is_parsed = False
            if basisset_info[basisset_name]['status'] == 'published':
                if basisset_info[basisset_name]['type'] == 'orbital':
                    self._run_for_orbital(basisset_info, alias)
                    is_parsed = True
                elif basisset_info[basisset_name]['type'] == 'dftorb':
                    self._run_for_dftorb(basisset_info, alias)
                    is_parsed = True
                else:
                    self._run_for_others(basisset_info, alias)

            if not is_parsed:
                logger.info(
                    "# pass: {name} ({status}; {bs_type})".format(name=basisset_name,
                                                                  status=basisset_info[basisset_name]['status'],
                                                                  bs_type=basisset_info[basisset_name]['type']))
        else:
            logger.error(
                "basisset \"{}\" is not found. ".format(basisset_name))

    def _sort_atom_symbols(self, atom_symbols):
        atom_table = [None] * bridge.PeriodicTable.get_num_of_atoms()
        for atom_symbol in atom_symbols:
            atom_num = bridge.PeriodicTable.get_atomic_number(atom_symbol)
            atom_table[atom_num] = atom_symbol
        return [x for x in atom_table if x != None]

    def _run_for_orbital(self, basisset_info, alias):
        """
        """
        db_basisset_name = self._alias2real(alias)
        logger.debug("orb: [{}] -> [{}]".format(alias, db_basisset_name))

        bs_content = self._download_gaussian94(basisset_info, db_basisset_name)
        logger.info("# reg: {name} ({status}; {bs_type})".format(name=db_basisset_name,
                                                                 status=basisset_info[db_basisset_name]['status'],
                                                                 bs_type=basisset_info[db_basisset_name]['type']))
        if self._debug:
            filepath = '{}-orbital.txt'.format(alias)
            filepath = filepath.replace(os.sep, '_')
            with open(filepath, 'wb') as f:
                # f = codecs.EncodedFile(f, 'utf-8')
                f.write(bridge.Utils.to_bytes(bs_content))

        # make basissets & basis2
        bsp = pdf.BasisSetParser()
        (basissets, basissets_density, basissets_xc) = bsp.parse(bs_content)
        # basissets is dict (key="atom-symbol", value=BasisSet object)

        basis2 = ""
        if db_basisset_name in self._cache['ProteinDF']:
            basis2 = self._cache['ProteinDF'][db_basisset_name]
        else:
            basis2 = self._get_basis2(db_basisset_name, basissets)
            self._cache['ProteinDF'][db_basisset_name] = basis2
            self._save_cache()
        if len(basis2) > 0:
            self._stdout.write('# {} -------------\n'.format(alias))
            self._stdout.write(basis2 + '\n')

        # uncontract shell basis set for Gridfree MRD
        if len(self._uncontract_shells) > 0:
            for i in range(len(self._uncontract_shells)):
                uncontract_shell = self._uncontract_shells[0:i+1]
                db_basisset_name_uncontract = db_basisset_name + uncontract_shell

                basis2_MRD = ""
                if db_basisset_name_uncontract in self._cache["ProteinDF-MRD"]:
                    basis2_MRD = self._cache["ProteinDF-MRD"][db_basisset_name_uncontract]
                else:
                    basis2_MRD = self._get_basis2_MRD(
                        db_basisset_name, basissets, uncontract_shell)
                    self._cache["ProteinDF-MRD"][db_basisset_name_uncontract] = basis2_MRD
                    self._save_cache()

                if len(basis2_MRD) > 0:
                    self._stdout.write(
                        '# {} -------------\n'.format(alias + "-" + uncontract_shell))
                    self._stdout.write(basis2_MRD + '\n')

    def _run_for_dftorb(self, basisset_info, alias):
        """
        """
        db_basisset_name = self._alias2real(alias)
        logger.debug("dft: alias[{}] -> {}".format(alias, db_basisset_name))

        bs_content = self._download_gaussian94(basisset_info, db_basisset_name)
        logger.info("# reg: {name} ({status}; {bs_type})".format(name=db_basisset_name,
                                                                 status=basisset_info[db_basisset_name]['status'],
                                                                 bs_type=basisset_info[db_basisset_name]['type']))

        if self._debug:
            filepath = '{}-dftorb.txt'.format(alias)
            filepath = filepath.replace(os.sep, '_')
            with open(filepath, 'wb') as f:
                # f = codecs.EncodedFile(f, 'utf-8')
                f.write(bridge.Utils.to_bytes(bs_content))

        # make basissets & basis2
        bsp = pdf.BasisSetParser()
        (basissets, basissets_density, basissets_xc) = bsp.parse(bs_content)
        # basissets is dict (key="atom-symbol", value=BasisSet object)

        basis2 = ""
        if db_basisset_name in self._cache["ProteinDF"]:
            basis2 = self._cache["ProteinDF"][db_basisset_name]
        else:
            basis2 = self._get_basis2(db_basisset_name, basissets)
            self._cache["ProteinDF"][db_basisset_name] = basis2
            self._save_cache()
        if len(basis2) > 0:
            self._stdout.write('# {} -------------\n'.format(alias))
            self._stdout.write(basis2 + '\n')

        basis2_dft = ""
        if db_basisset_name in self._cache["ProteinDF-DFT"]:
            basis2_dft = self._cache["ProteinDF-DFT"][db_basisset_name]
        else:
            basis2_dft = self._get_basis2(
                db_basisset_name, basissets_density, basissets_xc)
            self._cache["ProteinDF-DFT"][db_basisset_name] = basis2_dft
            self._save_cache()
        if len(basis2_dft) > 0:
            self._stdout.write('# {} (DFT)---------\n'.format(alias))
            self._stdout.write(basis2_dft + '\n')

        # uncontract shell basis set for Gridfree MRD
        if len(self._uncontract_shells) > 0:
            for i in range(len(self._uncontract_shells)):
                uncontract_shell = self._uncontract_shells[0:i+1]
                db_basisset_name_uncontract = db_basisset_name + uncontract_shell

                basis2_MRD = ""
                if db_basisset_name_uncontract in self._cache["ProteinDF-MRD"]:
                    basis2_MRD = self._cache["ProteinDF-MRD"][db_basisset_name_uncontract]
                else:
                    basis2_MRD = self._get_basis2_MRD(
                        db_basisset_name, basissets, uncontract_shell)
                    self._cache["ProteinDF-MRD"][db_basisset_name_uncontract] = basis2_MRD
                    self._save_cache()

                if len(basis2_MRD) > 0:
                    self._stdout.write(
                        '# {} -------------\n'.format(alias + "-" + uncontract_shell))
                    self._stdout.write(basis2_MRD + '\n')

    def _run_for_others(self, basisset_info, alias):
        basisset_name = self._alias2real(alias)
        bs_content = self._download_gaussian94(basisset_info, basisset_name)
        logger.info(
            "# reg: {name} ({status}; {bs_type})".format(name=basisset_name,
                                                         status=basisset_info[basisset_name]['status'],
                                                         bs_type=basisset_info[basisset_name]['type']))

    def _download_gaussian94(self, basisset_info, basisset_name):
        if (self._force_download == True) or (basisset_name not in self._cache['Gaussian94']):
            logger.info(
                "# get: {name} ({status}; {bs_type})".format(name=basisset_name,
                                                             status=basisset_info[basisset_name]['status'],
                                                             bs_type=basisset_info[basisset_name]['type']))
            bs_content = self._get_bs_data(basisset_name,
                                           basisset_info[basisset_name]['url'],
                                           basisset_info[basisset_name]['elements'])
            self._cache['Gaussian94'][basisset_name] = bs_content
            self._save_cache()
        return self._cache['Gaussian94'][basisset_name]

    def _load_alias(self, path):
        self._alias = {}
        re_alias = re.compile("^\s*(\S.*?)\s*=\s*(\S.*?)\s*$")

        if path != None:
            with open(path, mode='r') as f:
                line = f.readline()
                while line:
                    line = line.rstrip()
                    if (len(line) > 0) and (line[0] != '#'):
                        re_result = re_alias.match(line)
                        if re_result != None:
                            alias_name, real_name = re_result.group(1, 2)
                            self._alias[alias_name] = real_name
                        else:
                            logger.error(
                                "illegal alias format: {}".format(line))
                    line = f.readline()

        # debug out
        for alias in self._alias.keys():
            logger.debug(
                "alias: [{}] -> [{}]".format(alias, self._alias[alias]))

    def _load_cache(self):
        if os.path.exists(self._cache_path):
            with open(self._cache_path, mode='rb') as f:
                self._cache = pickle.load(f)

        # cache clean
        self._cache["ProteinDF"] = {}
        self._cache["ProteinDF-DFT"] = {}
        self._cache["ProteinDF-MRD"] = {}

    def _save_cache(self):
        with open(self._cache_path, mode='wb') as f:
            pickle.dump(self._cache, f)

    def _get_content(self, url, params={}):
        cache_key = tuple((url, tuple(params.items())))
        if cache_key not in self._cache:
            response = requests.get(url, params=params)
            html = response.text.encode(response.encoding)
            html = bridge.Utils.to_unicode(html)
            self._cache[cache_key] = html
            self._save_cache()

        return self._cache[cache_key]

    def _get_basissets_info(self):
        '''download basisset information

        retval: basisset information dict.
        '''

        url = "https://bse.pnl.gov/bse/portal/user/anon/js_peid/11535052407933/panel/Main/template/content"
        index_html = self._get_content(url)

        db_info = self._parse_basissets_index(index_html)

        return db_info

    def _parse_basissets_index(self, html):
        '''
        0 - path to xml file
        1 - basis set name
        2 - categorization: "dftcfit", "dftorb", "dftxfit", "diffuse",
                             "ecporb","effective core potential", "orbital",
                             "polarization", "rydberg", or "tight"
        3 - parameterized elements by symbol e.g. '[H, He, B, C, N, O, F, Ne]'
        4 - curation status; only 'published' is trustworthy
        5 - boolean: has ECP
        6 - boolean: has spin
        '''
        answer = {}

        re_bs = re.compile(
            'new basisSet\("?(.+?)"?,\s*"(.+?)",\s*"?(.+?)"?,\s*"?\[(.+?)\]"?,\s*"?(.+?)"?,\s*"?(.+?)"?,\s*"?(.+?)"?,\s*"?(.+?)"?.*\)')
        lines = re_bs.findall(html)
        for l in lines:
            url = l[0]
            name = l[1]
            bs_type = l[2]
            elements = l[3]
            elements = filter(lambda w: len(w) > 0, re.split(', ', elements))
            status = l[4]
            has_ECP = (l[5][0].upper() == 'T')
            has_spin = (l[6][0].upper() == 'T')
            last_modified = l[7]

            answer[name] = {'url': url,
                            'type': bs_type,
                            'elements': elements,
                            'status': status,
                            'has_ESP': has_ECP,
                            'has_spin': has_spin,
                            'last_modified': last_modified
                            }

        return answer

    def _get_bs_data(self, name, bsurl, elements, text_format='Gaussian94'):
        contraction = 'True'
        url = "https://bse.pnl.gov:443/bse/portal/user/anon/js_peid/11535052407933/action/portlets.BasisSetAction/template/courier_content/panel/Main/"
        url += "/eventSubmit_doDownload/true"

        params = {'bsurl': bsurl,
                  'bsname': name,
                  'elts': " ".join(elements),
                  'format': text_format,
                  'minimize': contraction}
        html = self._get_content(url, params=params)

        soup = BeautifulSoup(html, "html.parser")
        text = soup.pre.get_text()

        return text

    def _get_basis2(self, basisset_name, basissets, basissets2=None):
        assert isinstance(basissets, dict)
        max_shell_type = self.MAX_SHELL_TYPE

        # escape and remove strings for basis set name
        basisset_name = self._escape_basisset_name(basisset_name)
        if basissets2 == None:
            basisset_name = "O-" + basisset_name
        else:
            basisset_name = "A-" + basisset_name

        basis2 = ''
        for atom_symbol in self._sort_atom_symbols(basissets.keys()):
            bs = basissets[atom_symbol]
            assert(isinstance(bs, pdf.BasisSet))
            bs.name = basisset_name + '.' + atom_symbol

            if (basissets2 != None) and (atom_symbol in basissets2):
                # for aux basis
                bs2 = basissets2[atom_symbol]
                assert(isinstance(bs2, pdf.BasisSet))
                if ((bs.max_shell_type in max_shell_type) and
                        (bs2.max_shell_type in max_shell_type)):
                    basis2 += str(bs) + '\n' + str(bs2) + '\n'
                else:
                    logger.warning('# {basisset_name} is not supported: (max shell type: "{shell_type}")'.format(basisset_name=bs.name,
                                                                                                                 shell_type=bs.max_shell_type))
            else:
                # for orbital basis
                if bs.max_shell_type in max_shell_type:
                    basis2 += str(bs) + '\n'
                else:
                    logger.warning('# {basisset_name} is not supported: (max shell type: "{shell_type}")'.format(basisset_name=bs.name,
                                                                                                                 shell_type=bs.max_shell_type))

        return basis2

    def _get_basis2_MRD(self, basisset_name, basissets, uncontract_shell):
        assert isinstance(basissets, dict)
        max_shell_type = self.MAX_SHELL_TYPE

        # escape and remove strings for basis set name
        basisset_name = "O-" + self._escape_basisset_name(basisset_name)

        basis2 = ''
        for atom_symbol in self._sort_atom_symbols(basissets.keys()):
            bs = basissets[atom_symbol]
            assert(isinstance(bs, pdf.BasisSet))
            bs_MRD = self._makeBasisSet_MRD(bs, uncontract_shell)
            bs_MRD.name = basisset_name + \
                '-{}'.format(uncontract_shell) + '.' + atom_symbol

            # for orbital basis
            if bs_MRD.max_shell_type in max_shell_type:
                basis2 += str(bs_MRD) + '\n'
            else:
                logger.warning('# {basisset_name} is not supported: (max shell type: "{shell_type}")'.format(basisset_name=bs_MRD.name,
                                                                                                             shell_type=bs_MRD.max_shell_type))

        return basis2

    def _makeBasisSet_MRD(self, basisset, uncontract_shell='spd'):
        assert isinstance(basisset, pdf.BasisSet)

        new_cgtos = []
        for cgto in basisset:
            shell_type = cgto.shell_type
            if shell_type in uncontract_shell:
                for pgto in cgto:
                    cgto1 = pdf.ContractedGTO(shell_type, 1)
                    cgto1[0] = pdf.PrimitiveGTO(pgto.exp, 1.0)
                    if cgto1 not in new_cgtos:
                        new_cgtos.append(cgto1)
            else:
                new_cgtos.append(pdf.ContractedGTO(cgto))

        # reverse replace!
        new_name = '-{}.'.format(uncontract_shell).join(
            basisset.name.rsplit('.', 1))
        answer = pdf.BasisSet(new_name, len(new_cgtos))
        for i in range(len(new_cgtos)):
            answer[i] = new_cgtos[i]

        return answer

    def _escape_basisset_name(self, name):
        # escape and remove strings for basis set name
        name = copy.copy(name)
        name_old = copy.copy(name)
        name = name.rstrip()
        name = name.replace(' (DFT Orbital)', '')
        name = re.sub(' \(Dunning.*\)', '', name)
        name = name.replace(' + ', '+')
        name = re.sub(' +$', '', name)
        name = name.rstrip()
        # name = name.replace(' ', '_')
        logger.debug("_get_basis2(): [{}] => [{}]".format(name_old, name))

        return name


def main():
    default_cache_path = os.path.join(os.environ['PDF_HOME'],
                                      'data',
                                      'basis2.cache')

    parser = argparse.ArgumentParser(description='make basis2 file')
    parser.add_argument('--list',
                        action='store_true',
                        help='show list')
    parser.add_argument('-a', '--alias',
                        nargs=1,
                        default=[None],
                        help='text file of alias list of basisset name')
    parser.add_argument('-u', '--uncontract',
                        nargs=1,
                        default=['spd'],
                        help='uncontract shells for gridfree')
    parser.add_argument('-c', '--cache',
                        nargs=1,
                        default=[default_cache_path],
                        help='cache file: default={}'.format(default_cache_path))
    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=[None],
                        help='output file path')
    parser.add_argument('-f', '--force-download',
                        action='store_true')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='show verbosely')
    parser.add_argument('-d', '--debug',
                        action='store_true',
                        help='debug use')
    parser.add_argument('basisset_name',
                        nargs='*',
                        help="basisset name")
    args = parser.parse_args()

    is_show_list = args.list
    basisset_names = args.basisset_name
    alias_path = args.alias[0]
    cache_path = args.cache[0]
    output_path = args.output[0]
    uncontract_shells = args.uncontract[0]
    force_download = args.force_download
    verbose = args.verbose
    debug = args.debug

    if verbose:
        logger.info('cache: {}'.format(cache_path))
        logger.info('alias: {}'.format(alias_path))
        logger.info('output: {}'.format(output_path))
        logger.info('uncontract_shells: {}'.format(uncontract_shells))
        logger.info('debug: {}'.format(debug))
        if len(basisset_names) > 0:
            logger.info('basisset: ')
            for basisset_name in basisset_names:
                logger.info('{} '.format(basisset_name))

    header = ""
    # header += '##### basis2 file created by {}. #####\n'.format(parser.prog)
    header += '#' * 80
    header += '\n'
    bs2 = MakeBasis2(alias_path=alias_path,
                     cache_path=cache_path,
                     output_path=output_path,
                     uncontract_shells=uncontract_shells,
                     force_download=force_download,
                     verbose=verbose,
                     debug=debug)

    if is_show_list:
        bs2.show_list()
    else:
        bs2.run(basisset_names=basisset_names,
                header=header)


if __name__ == '__main__':
    if os.path.exists("config.ini"):
        logging.config.fileConfig("config.ini",
                                  disable_existing_loggers=False)
    main()

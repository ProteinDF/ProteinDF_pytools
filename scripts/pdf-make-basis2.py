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

import logging.config
from logging import getLogger, DEBUG, INFO
logger = getLogger(__name__)

from bs4 import BeautifulSoup
import pdfbridge
import pdfpytools as pdf


class MakeBasis2(object):
    def __init__(self,
                 alias_path=None,
                 cache_path=None,
                 output_path=None,
                 uncontract_shells='spd',
                 force_download = False,
                 verbose=False,
                 debug=False):
        self._uncontract_shells = uncontract_shells
        self._force_download = force_download
        self._debug = debug

        self._load_alias(alias_path)
        
        self._cache = {}
        self._cache_path = cache_path

        self._stdout = sys.stdout
        if output_path != None:
            self._stdout = open(output_path, 'a')
        # self._stderr = open(os.devnull, 'a')
        #if verbose:
        #    self._stderr = sys.stderr

        self._load_cache()

    def _alias2real(self, alias):
        real_name = alias
        if alias in self._alias.keys():
            real_name = self._alias[alias]
        return real_name
        
    def show_list(self):
        bs_index = self._get_basissets_index()
        for bs in bs_index.keys():
            self._stdout.write(bs + "\n")
        
    def run(self, basisset_list=[], header=""):
        # output header
        self._stdout.write(header + '\n')

        # make basis set
        bs_index = self._get_basissets_index()

        self._cache.setdefault('Gaussian94', {})
        self._cache.setdefault('ProteinDF', {})

        if len(basisset_list) == 0:
            basisset_list = bs_index.keys()

        num_of_items = len(basisset_list)
        count = 0
        for bs_name in basisset_list:
            logger.info('### {0:>4}/{1:>4}'.format(count +1, num_of_items))
            self._run_basisset(bs_index, bs_name)
            count += 1
        assert(num_of_items == count)


    def _run_basisset(self, bs_index, alias):
        bs_name = self._alias2real(alias)

        logger.debug("run getting basisset: {}".format(bs_name))
        if bs_name in bs_index:
            is_parsed = False
            if bs_index[bs_name]['status'] == 'published':
                if bs_index[bs_name]['type'] == 'orbital':
                    self._run_for_orbital(bs_index, alias)
                    is_parsed = True
                elif bs_index[bs_name]['type'] == 'dftorb':
                    self._run_for_dftorb(bs_index, alias)
                    is_parsed = True
                else:
                    self._run_for_others(bs_index, alias)
        
            if not is_parsed:
                logger.info(
                    "# pass: {name} ({status}; {bs_type})".format(name=bs_name,
                                                                  status=bs_index[bs_name]['status'],
                                                                    bs_type=bs_index[bs_name]['type']))
        else:
            logger.error("basisset \"{}\" is not found. ".format(bs_name))
        
    
    def _run_for_orbital(self, bs_index, alias):
        bs_name = self._alias2real(alias)
        logger.debug("orb: [{}] -> [{}]".format(alias, bs_name))

        bs_content = self._download_gaussian94(bs_index, bs_name)
        logger.info("# reg: {name} ({status}; {bs_type})".format(name=bs_name,
                                                                 status=bs_index[bs_name]['status'],
                                                                 bs_type=bs_index[bs_name]['type']))
        if self._debug:
            filepath = '{}-orbital.txt'.format(alias)
            filepath = filepath.replace(os.sep, '_')
            with open(filepath, 'wb') as f:
                # f = codecs.EncodedFile(f, 'utf-8')
                f.write(pdfbridge.Utils.to_bytes(bs_content))

        # make basissets & basis2
        bsp = pdf.BasisSetParser()
        (basissets, basissets_density, basissets_xc) = bsp.parse(bs_content)
        basis2 = self._get_basis2('O-' + bs_name, basissets)
        self._cache['ProteinDF'][bs_name] = basis2
        self._save_cache()
        self._stdout.write('# {} -------------\n'.format(alias))
        if len(basis2) > 0:
            self._stdout.write(basis2 + '\n')

        # uncontract shell basis set for Gridfree MRD
        if len(self._uncontract_shells) > 0:
            for i in range(len(self._uncontract_shells)):
                uncontract_shell = self._uncontract_shells[0:i+1]
                uncontract_basissets = {}
                for atom_symbol, bs in basissets.items():
                    uncontract_bs = self._bs2MRD(bs, uncontract_shell)
                    uncontract_basissets[atom_symbol] = uncontract_bs
                    uncontract_basis_name = 'O-' + bs_name + '-{}'.format(uncontract_shell)
                    uncontract_basis2 = self._get_basis2(uncontract_basis_name, uncontract_basissets)
                    self._cache['ProteinDF'][uncontract_basis_name] = uncontract_basis2
                    self._save_cache()
                    if len(uncontract_basis2) > 0:
                        self._stdout.write(uncontract_basis2 + '\n')
                
    def _run_for_dftorb(self, bs_index, alias):
        bs_name = self._alias2real(alias)
        logger.debug("dft: alias[{}] -> {}".format(alias, bs_name))
        
        bs_content = self._download_gaussian94(bs_index, bs_name)
        logger.info(
            "# reg: {name} ({status}; {bs_type})".format(name=bs_name,
                                                         status=bs_index[bs_name]['status'],
                                                         bs_type=bs_index[bs_name]['type']))
        
        if self._debug:
            filepath = '{}-dftorb.txt'.format(alias)
            filepath = filepath.replace(os.sep, '_')
            with open(filepath, 'wb') as f:
                # f = codecs.EncodedFile(f, 'utf-8')
                f.write(pdfbridge.Utils.to_bytes(bs_content))
        
        # make basissets & basis2
        bsp = pdf.BasisSetParser()
        (basissets, basissets_density, basissets_xc) = bsp.parse(bs_content)
        basis2 = self._get_basis2('O-' + bs_name, basissets)
        basis2 += self._get_basis2('A-' + bs_name, basissets_density, basissets_xc)
        self._cache['ProteinDF'][bs_name] = basis2
        self._save_cache()
        self._stdout.write('# {} -------------\n'.format(alias))
        self._stdout.write(basis2 + '\n')

        # uncontract shell basis set for Gridfree MRD
        if len(self._uncontract_shells) > 0:
            for i in range(len(self._uncontract_shells)):
                uncontract_shell = self._uncontract_shells[0:i+1]
                uncontract_basissets = {}
                for atom_symbol, bs in basissets.items():
                    uncontract_bs = self._bs2MRD(bs, uncontract_shell)
                    uncontract_basissets[atom_symbol] = uncontract_bs
                    uncontract_basis_name = 'O-' + bs_name + '-{}'.format(uncontract_shell)
                    uncontract_basis2 = self._get_basis2(uncontract_basis_name, uncontract_basissets)
                    self._cache['ProteinDF'][uncontract_basis_name] = uncontract_basis2
                    self._save_cache()
                    if len(uncontract_basis2) > 0:
                        self._stdout.write(uncontract_basis2 + '\n')
                
        
    def _run_for_others(self, bs_index, alias):
        bs_name = self._alias2real(alias)
        bs_content = self._download_gaussian94(bs_index, bs_name)
        logger.info(
            "# reg: {name} ({status}; {bs_type})".format(name=bs_name,
                                                         status=bs_index[bs_name]['status'],
                                                         bs_type=bs_index[bs_name]['type']))

        
    def _download_gaussian94(self, bs_index, bs_name):
        if (self._force_download == True) or (bs_name not in self._cache['Gaussian94']):
            logger.info(
                "# get: {name} ({status}; {bs_type})".format(name=bs_name,
                                                             status=bs_index[bs_name]['status'],
                                                             bs_type=bs_index[bs_name]['type']))
            bs_content = self._get_bs_data(bs_name,
                                           bs_index[bs_name]['url'],
                                           bs_index[bs_name]['elements'])
            self._cache['Gaussian94'][bs_name] = bs_content
            self._save_cache()
        return self._cache['Gaussian94'][bs_name]
        

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
                            logger.error("illegal alias format: {}".format(line))
                    line = f.readline()

        # debug out
        for alias in self._alias.keys():
            logger.debug("alias: [{}] -> [{}]".format(alias, self._alias[alias]))
    
    
    def _load_cache(self):
        if os.path.exists(self._cache_path):
            with open(self._cache_path, mode='rb') as f:
                self._cache = pickle.load(f)
            
    def _save_cache(self):
        with open(self._cache_path, mode='wb') as f:
            pickle.dump(self._cache, f)

            
    def _get_content(self, url, params={}):
        cache_key = tuple((url, tuple(params.items())))
        if cache_key not in self._cache:
            response = requests.get(url, params=params)
            html = response.text.encode(response.encoding)
            html = pdfbridge.Utils.to_unicode(html)
            self._cache[cache_key] = html
            self._save_cache()
        
        return self._cache[cache_key]

    
    def _get_basissets_index(self):
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

        re_bs = re.compile('new basisSet\("?(.+?)"?,\s*"(.+?)",\s*"?(.+?)"?,\s*"?\[(.+?)\]"?,\s*"?(.+?)"?,\s*"?(.+?)"?,\s*"?(.+?)"?,\s*"?(.+?)"?.*\)')
        lines = re_bs.findall(html)
        for l in lines:
            url = l[0]
            name = l[1]
            bs_type = l[2]
            elements = l[3]
            elements = filter(lambda w: len(w)>0, re.split(', ', elements))
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


    def _bs2MRD(self, basisset, uncontract_shell='spd'):
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

        new_name = '-{}.'.format(uncontract_shell).join(basisset.name.rsplit('.', 1)) # reverse replace!
        answer = pdf.BasisSet(new_name, len(new_cgtos))
        for i in range(len(new_cgtos)):
            answer[i] = new_cgtos[i]

        return answer

    def _get_basis2(self, basissets_name, basissets, basissets2=None):
        assert isinstance(basissets, dict)
        max_shell_type = 'spdfg'

        # escape and remove strings for basis set name
        basissets_name_old = basissets_name
        basissets_name = basissets_name.rstrip()
        basissets_name = basissets_name.replace(' (DFT Orbital)', '')
        basissets_name = re.sub(' \(Dunning.*\)', '', basissets_name)
        basissets_name = basissets_name.replace(' + ', '+')
        basissets_name = re.sub(' +$', '', basissets_name)
        basissets_name = basissets_name.rstrip()
        #basissets_name = basissets_name.replace(' ', '_')
        logger.debug("_get_basis2(): [{}] => [{}]".format(basissets_name_old, basissets_name))
        
        basis2 = ''
        for atom_symbol, bs in basissets.items():
            bs.name = basissets_name + '.' + atom_symbol
            if (basissets2 != None) and (atom_symbol in basissets2):
                # for aux basis
                bs2 = basissets2[atom_symbol]
                if ((bs.max_shell_type in max_shell_type) and
                    (bs2.max_shell_type in max_shell_type)):
                    basis2 += str(bs) + '\n' + str(bs2) + '\n'
                else:
                    logger.warning('# {bs_name} is not supported: (max shell type: "{shell_type}")'.format(bs_name=bs.name,
                                                                                                           shell_type=bs.max_shell_type))
            else:
                # for orbital basis
                if bs.max_shell_type in max_shell_type:
                    basis2 += str(bs) + '\n'
                else:
                    logger.warning('# {bs_name} is not supported: (max shell type: "{shell_type}")'.format(bs_name=bs.name,
                                                                                                           shell_type=bs.max_shell_type))

        return basis2

    
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
    basisset_list = args.basisset_name
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
        if len(basisset_list) > 0:
            logger.info('basisset: ')
            for bs_name in basisset_list:
                logger.info('{} '.format(bs_name))
        

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
        bs2.run(basisset_list=basisset_list,
                header=header)
    
    
if __name__ == '__main__':
    if os.path.exists("config.ini"):
        logging.config.fileConfig("config.ini",
                                  disable_existing_loggers=False)
    main()
        


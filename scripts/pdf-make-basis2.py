#!/usr/bin/env python
# -*- coding: utf-8 -*-

import os
import sys
import argparse
import requests
import re
import pprint
import pickle

from bs4 import BeautifulSoup
import pdfpytools as pdf


class MakeBasis2(object):
    def __init__(self, cache_path=None, output_path=None, verbose=False, header=""):
        self._cache = {}
        self._cache_path = cache_path
        
        self._load_cache()
        self.run(output_path, verbose, header)
        
    def run(self, output_path=None, verbose=False, header=""):
        stdout = sys.stdout
        if output_path != None:
            stdout = open(output_path, 'w')

        stderr = open(os.devnull, 'w')
        if verbose:
            stderr = sys.stderr

        # output header
        stdout.write(header + '\n')

        # make basis set
        bs_index = self._get_basissets_index()

        self._cache.setdefault('Gaussian94', {})
        self._cache.setdefault('ProteinDF', {})
        num_of_items = len(bs_index.keys())
        count = 1
        for bs_name in bs_index.keys():
            if (bs_index[bs_name]['status'] == 'published' and
                bs_index[bs_name]['type'] == 'orbital') :
                if bs_name not in self._cache['Gaussian94']:
                    stderr.write("# {count}/{max_count}> get: {name} ({status}; {bs_type})\n".format(count=count,
                                                                                                     max_count=num_of_items,
                                                                                                     name=bs_name,
                                                                                                     status=bs_index[bs_name]['status'],
                                                                                                     bs_type=bs_index[bs_name]['type']))
                    bs_content = self._get_bs_data(bs_name,
                                                   bs_index[bs_name]['url'],
                                                   bs_index[bs_name]['elements'])
                    self._cache['Gaussian94'][bs_name] = bs_content
                    self._save_cache()

                stderr.write("# {count}/{max_count}> reg: {name} ({status}; {bs_type})\n".format(count=count,
                                                                                                 max_count=num_of_items,
                                                                                                 name=bs_name,
                                                                                                 status=bs_index[bs_name]['status'],
                                                                                                 bs_type=bs_index[bs_name]['type']))
                bs_content = self._cache['Gaussian94'][bs_name]

                bsp = pdf.BasisSetParser()
                basissets = bsp.parse(bs_content)

                basis2 = ''
                for atom_symbol, bs in basissets.items():
                    bs.name = bs_name + '.' + atom_symbol
                    if bs.max_shell_type in ('s', 'p', 'd', 'f', 'g'):
                        basis2 += str(bs) + '\n'
                    else:
                        stderr.write('# {bs_name} is not supported, because the max shell type is "{shell_type}"\n'.format(bs_name=bs.name,
                                                                                                                         shell_type=bs.max_shell_type))
                    
                self._cache['ProteinDF'][bs_name] = basis2
                self._save_cache()

                stdout.write('# {} -------------\n'.format(bs_name))
                stdout.write(basis2 + '\n')
                
            else:
                stderr.write("# {count}/{max_count}> pass: {name} ({status}; {bs_type})\n".format(count=count,
                                                                                                max_count=num_of_items,
                                                                                                name=bs_name,
                                                                                                status=bs_index[bs_name]['status'],
                                                                                                bs_type=bs_index[bs_name]['type']))
                
            count += 1

            
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
            self._cache[cache_key] = html
            self._save_cache()
        
        return self._cache[cache_key]

    
    def _get_basissets_index(self):
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

    
def main():
    #default_cache_path = os.path.join(os.path.dirname(sys.modules['pdfpytools'].__file__),
    #                                  'basis2.cache')
    default_cache_path = os.path.join(os.environ['PDF_HOME'],
                                      'data',
                                      'basis2.cache')
    
    parser = argparse.ArgumentParser(description='make basis2 file')
    parser.add_argument('-c', '--cache',
                        nargs=1,
                        default=[default_cache_path],
                        help='cache file: default={}'.format(default_cache_path))
    parser.add_argument('-o', '--output',
                        nargs=1,
                        default=[None],
                        help='output file path')
    parser.add_argument('-v', '--verbose',
                        action='store_true',
                        help='show verbosely')
    args = parser.parse_args()

    cache_path = args.cache[0]
    output_path = args.output[0]
    verbose = args.verbose
    
    header = '# basis2 file created by {}.\n'.format(parser.prog)
    header += '#' * 80
    header += '\n'
    bs2 = MakeBasis2(cache_path=cache_path,
                     output_path=output_path,
                     verbose=verbose,
                     header=header)
    
if __name__ == '__main__':
    main()
        


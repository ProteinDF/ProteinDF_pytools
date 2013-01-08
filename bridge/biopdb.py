#!/usr/bin/env python
# -*- coding: utf-8 -*-

import sys
import os
import optparse
import re

import bridge
#from atomgroup import *
#from position import *

class Pdb(object):
    """
    """
    def __init__(self, file_path = None):
        """
        create empty PDB object
        """
        self._data = {}
        if (file_path):
            self.load(file_path)

    def __get_debug(self):
        if not '_debug' in self.__dict__:
            self._debug = False
        return self._debug

    def __set_debug(self, yn):
        self._debug = yn

    debug = property(__get_debug, __set_debug)
            
    #def get_number_of_items(self):
    #    return len(self._data)

    #def get_name(self, index):
    #    answer = None
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        answer = self._data[index].get('name', None)
    #    return answer

    #def set_name(self, index, name):
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        self._data[index]['name'] = name

    #def get_element(self, index):
    #    answer = None
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        answer = self._data[index].get('element', None)
    #    return answer

    #def set_element(self, index, element):
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        self._data[index]['element'] = element

    #def get_posision(self, index):
    #    answer = None
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        answer = self._data[index].get('coord', None)
    #    return answer

    #def set_position(self, index, position):
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        self._data[index]['coord'] = position

    #def get_charge(self, index):
    #    answer = 0
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        answer = self._data[index].get('charge', 0)
    #    return answer

    #def set_charge(self, index, charge):
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        self._data[index]['charge'] = charge

    #def get_occupancy(self, index):
    #    answer = None
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        answer = self._data[index].get('occupancy', None)
    #    return answer

    #def get_temp_factor(self, index):
    #    answer = None
    #    if ((0 <= index) and (index < self.get_number_of_items())):
    #        answer = self._data[index].get('temp_factor', None)
    #    return answer
        
    #def set_temp_factor(self, serial, temp_factor):
    #    serial = int(serial)
    #    self._data[serial]['temp_factor'] = temp_factor

    def renumber(self):
        for model_serial, model in self._data.iteritems():
            for index in range(len(model)):
                model[index]['serial'] = index +1

    def load(self, file_path):
        if (os.path.isfile(file_path) != True):
            return

        model_serial = 1
        self._data.setdefault(model_serial, [])

        fin = open(file_path, 'r')
        while True:
            line = fin.readline()
            if (len(line) == 0):
                break
            line = line.rstrip('\n')
            
            #if (len(line) != 80):
            #    continue
            if (self.debug == True):
                print(line)

            record_name = line[0:6]
            if ((record_name == 'ATOM  ') or (record_name == 'HETATM')):

                if (len(line) < 80):
                    line = line + (' ' * (80 - len(line)))

                serial = int(line[6:11])
                name = line[12:16]
                alt_loc = line[16]
                res_name = line[17:20]
                chain_id = line[21]
                res_seq = line[22:26]
                i_code = line[26]
                coord_x = line[30:38]
                coord_y = line[38:46]
                coord_z = line[46:54]
                occupancy = line[54:60].strip()
                temp_factor = line[60:66].strip()
                element = line[76:78].strip()
                charge = line[78:80].strip()
            
                item = {}
                item['serial'] = serial
                item['record_name'] = record_name
                item['name'] = name
                item['alt_loc'] = alt_loc
                item['res_name'] = res_name
                item['chain_id'] = chain_id
                item['res_seq'] = int(res_seq)
                item['i_code'] = i_code
                item['coord'] = [float(coord_x), float(coord_y), float(coord_z)]

                if (len(occupancy) != 0):
                    item['occupancy'] = float(occupancy)
                else:
                    item['occupancy'] = 1.0
                if (len(temp_factor) != 0):
                    item['temp_factor'] = float(temp_factor)
                else:
                    item['temp_factor'] = 0.0

                if (len(element) != 0):
                    item['element'] = element
                else:
                    # TODO: テーブルを持って変換するように変更
                    if (name == ' MG '):
                        element = 'Mg'
                    elif (name == 'FE  '):
                        element = 'Fe'
                    else:
                        element = name[1]
                    item['element'] = element

                if (len(charge) != 0):
                    item['charge'] = charge
                else:
                    item['charge'] = "  "

                self._data[model_serial].append(item)
                continue
            elif (record_name == 'MODEL '):
                serial = int(line[10:14])
                model_serial = serial
                self._data.setdefault(model_serial, [])
                continue
            elif (record_name == 'TER   '):
                if (len(line) < 27):
                    line = line + (' ' * (27 - len(line)))
                serial = int(line[6:11]) if line[6:11].isdigit() else 0
                resname = line[17:20]
                chain_id = line[21]
                res_seq = line[22:26]
                i_code = line[26]
                item = {}
                item['serial'] = serial
                item['record_name'] = record_name
                item['res_name'] = res_name
                item['chain_id'] = chain_id
                item['res_seq'] = res_seq
                item['i_code'] = i_code

                self._data[model_serial].append(item)
                continue

    def get_atom_group(self):
        """
        return AtomGroup object
        """
        root = bridge.AtomGroup()

        for model_serial, model_items in self._data.iteritems():
            model = bridge.AtomGroup()
            model_name = "model_%d" % (model_serial)
            model.name = model_name
            
            for index in range(len(model_items)):
                item = model_items[index]
                record_name = item['record_name']
                serial = item['serial']

                if ((record_name == 'ATOM  ') or (record_name == 'HETATM')):
                    name = item['name'].strip()
                    alt_loc = item['alt_loc']
                    res_name = item['res_name']
                    chain_id = item['chain_id']
                    res_seq = str(item['res_seq'])
                    i_code = item['i_code']
                    coord = item['coord']
                    occupancy = item.get('occupancy', 1.0)
                    temp_factor = item.get('temp_factor', 0.0)
                    element = item.get('element', 'X')
                    charge = item.get('charge', 0.0)
                
                    if (chain_id == " "):
                        chain_id = "_"
                        
                    if (model.has_group(chain_id) == False):
                        chain = bridge.AtomGroup()
                        chain.name = chain_id
                        model.set_group(chain_id, chain)
                    if (model[chain_id].has_group(res_seq) == False):
                        residue = bridge.AtomGroup()
                        residue.name = res_name
                        model[chain_id].set_group(res_seq, residue)

                    atom = bridge.Atom()
                    atom.symbol = element
                    atom.position = bridge.Position(coord)
                    atom.name = name
                    model[chain_id][res_seq].set_atom(serial, atom)
                    
            root.set_group(model_name, model)
        return root

    def set_by_atomgroup(self, atomgroup, set_b_factor=None):
        assert(isinstance(atomgroup, AtomGroup))

        re_model_serial = re.compile("^model_(\d+)")
        self._data = {}
        item = {}
        model_serial = 1
        for model_key, model in atomgroup.groups():
            match_obj = re_model_serial.match(model_key)
            if (match_obj != None):
                model_serial = int(match_obj.group(1))
            self._data.setdefault(model_serial, [])
                
            for chain_id, chain in model.groups():
                if (chain_id != "_"):
                    item['chain_id'] = chain_id
                else:
                    item['chain_id'] = " "

                for res_seq, residue in chain.groups():
                    item["res_seq"] = res_seq
                    item["res_name"] = residue.name

                    for serial, atom in residue.atoms():
                        item["record_name"] = "ATOM  "
                        item["serial"] = serial
                        item["alt_loc"] = " "
                        item["i_code"] = " "
                        name = atom.name
                        if (len(name) != 4):
                            name = self.__match_name_table(atom.name,
                                                           atom.symbol)
                        item["name"] = name
                        item["coord"] = atom.position.xyz

                        if set_b_factor == 'charge':
                            item['temp_factor'] = atom.charge
                        
                        self._data[model_serial].append(copy.deepcopy(item))
        self._sort_by_serial()
        
    def _sort_by_serial(self):
        for model_serial, model in self._data.iteritems():
            model.sort(cmp = lambda x, y: cmp(int(x['serial']), int(y['serial'])))

    def __str__(self):
        occupancy = 1.0
        output = ""
        for model_serial, model in self._data.iteritems():
            output += "MODEL     %4d\n" % (model_serial)
            for index in range(len(model)):
                item = model[index]
                record_name = item['record_name']
                serial = int(item['serial'])
                if ((record_name == 'ATOM  ') or (record_name == 'HETATM')):
                    name = item['name']
                    alt_loc = item['alt_loc']
                    res_name = item['res_name']
                    chain_id = item['chain_id']
                    res_seq = int(item['res_seq'])
                    i_code = item['i_code']
                    coord = item['coord']
                    occupancy = item.setdefault('occupancy', 1.0)
                    temp_factor = item.setdefault('temp_factor', 1.0)
                    element = item.setdefault('element', '  ')
                    charge = item.setdefault('charge', '  ')
                    line = "ATOM  %5d %4s%c%3s %c%4d%c   %8.3f%8.3f%8.3f%6.2f%6.2f          %2s%2s\n" % (
                        serial, name, alt_loc,
                        res_name, chain_id, res_seq, i_code,
                        coord[0], coord[1], coord[2],
                        occupancy, temp_factor, element.upper(), charge)
                    output += line
                elif (record_name == 'TER   '):
                    line = "TER   %5d      %3s %c%4d%c\n" % (serial, res_name, chain_id, res_seq, i_code)
                    output += line
        return output

    def __match_name_table(self, name, symbol):
        assert(isinstance(name, str) == True)
        assert(isinstance(symbol, str) == True)
        capital_name = name.upper()
        capital_symbol = symbol.upper()
        if (capital_symbol == "NA"):
            name = "NA  "
        elif (capital_symbol == "CL"):
            name = "CL  "
        elif (capital_name[0] == capital_symbol[0]):
            name = " %s" % (name)
        name += " " * (4 - len(name))
        return name


def main():
    # initialize

    # parse args
    parser = optparse.OptionParser(usage = "%prog [options] PDB_FILE",
                                   version = "%prog 1.0")
    parser.add_option("-o", "--output", dest = "output_path",
                      help = "PDB output file", metavar = "FILE")
    parser.add_option("-v", "--verbose", dest = "verbose",
                      action="store_false", default = False,
                      help = "print message")
    (opts, args) = parser.parse_args()
        
    if (len(args) == 0):
        parser.print_help()
        sys.exit(1)

    # setting
    file_path = args[0]
    verbose = opts.verbose

    # 
    pdb_obj = Pdb(file_path)
    print(pdb_obj)

    # end

if __name__ == '__main__':
    main()


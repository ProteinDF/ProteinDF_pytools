#!/usr/bin/env python
# -*- coding: utf-8 -*-

class PeriodicTable(object):
    """
    >>> PeriodicTable.get_symbol(1)
    'H'
    >>> PeriodicTable.get_symbol(20)
    'Ca'
    >>> PeriodicTable.get_atomic_number('C')
    6
    >>> PeriodicTable.get_atomic_number('Cu')
    29
    """
    __table = [
        'X',
        'H', 'He',
        'Li', 'Be', 'B', 'C', 'N', 'O', 'F', 'Ne',
        'Na', 'Mg', 'Al', 'Si', 'P', 'S', 'Cl', 'Ar',
        'K', 'Ca', 'Sc', 'Ti', 'V', 'Cr', 'Mn', 'Fe', 'Co', 'Ni', 'Cu', 'Zn', 'Ga', 'Ge', 'As', 'Se', '', 'Kr',
        'Rb', 'Sr', 'Y', 'Zr', 'Nb', 'Mo', 'Tc', 'Ru', 'Rh', 'Pd', 'Ag', 'Cd', 'In', 'Sn', 'Sb', 'Te', 'I', 'Xe',
        'Cs', 'Ba',
        'La', 'Ce', 'Pr', 'Nd', 'Pm', 'Sm', 'Eu', 'Gd', 'Tb', 'Dy', 'Ho', 'Er', 'Tm', 'Yb', 'Lu',
        'Hf', 'Ta', 'W', 'Re', 'Os', 'Ir', 'Pt', 'Au', 'Hg', 'Tl', 'Pb', 'Bi', 'Po', 'At', 'Rn',
        'Fr', 'Ra',
        'Ac', 'Th', 'Pa', 'U', 'Np', 'Pu', 'Am', 'Cm', 'Bk', 'Cf', 'Es', 'Fm', 'Md', 'No', 'Lr',
        'Rf', 'Db', 'Sg', 'Bh', 'Hs', 'Mt', 'Ds', 'Rg', 'Cn', 'Uut', 'Uuq', 'Uup', 'Uuh', 'Uus', 'Uuo'
        ]


    @staticmethod
    def get_symbol(atomic_number):
        try:
            answer = PeriodicTable.__table[atomic_number]
            return answer
        except:
            print("ERROR @PeriodicTable::get_symbol(): not found input:%s." % (atomicNumber))
            raise


    @staticmethod
    def get_atomic_number(symbol):
        try:
            symbol = str(symbol)
        except:
            raise

        symbol = symbol.lower()
        symbol = symbol.capitalize()
        
        try:
            answer = PeriodicTable.__table.index(symbol)
            return answer
        except ValueError:
            print("ERROR @PeriodicTable::get_atomic_number(): not found symbol:%s." % (symbol))
            raise
        except:
            raise

    
if __name__ == "__main__":
    import doctest
    doctest.testmod()

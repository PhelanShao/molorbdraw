"""
Constants for molecular visualization and analysis.
"""
from enum import Enum

# Unit conversion constants
BOHR_TO_ANGSTROM = 0.529177249  # 1 Bohr = 0.529177249 Å

class BondType(Enum):
    """Bond types enumeration"""
    SINGLE = 1
    DOUBLE = 2
    TRIPLE = 3
    AROMATIC = 1.5

# Extended periodic table information: (symbol, radius, color)
ELEMENTS = {
    1: ('H', 0.31, (1.0, 1.0, 1.0)),      # Hydrogen
    2: ('He', 0.28, (0.85, 1.0, 1.0)),    # Helium
    3: ('Li', 1.28, (0.8, 0.5, 1.0)),     # Lithium
    4: ('Be', 0.96, (0.76, 1.0, 0.0)),    # Beryllium
    5: ('B', 0.84, (1.0, 0.71, 0.71)),    # Boron
    6: ('C', 0.76, (0.5, 0.5, 0.5)),      # Carbon
    7: ('N', 0.71, (0.05, 0.05, 1.0)),    # Nitrogen
    8: ('O', 0.66, (1.0, 0.05, 0.05)),    # Oxygen
    9: ('F', 0.57, (0.7, 1.0, 1.0)),      # Fluorine
    10: ('Ne', 0.58, (0.7, 0.89, 0.96)),  # Neon
    11: ('Na', 1.66, (0.67, 0.36, 0.95)), # Sodium
    12: ('Mg', 1.41, (0.54, 1.0, 0.0)),   # Magnesium
    13: ('Al', 1.21, (0.75, 0.65, 0.65)), # Aluminum
    14: ('Si', 1.11, (0.5, 0.6, 0.6)),    # Silicon
    15: ('P', 1.07, (1.0, 0.5, 0.0)),     # Phosphorus
    16: ('S', 1.05, (1.0, 1.0, 0.19)),    # Sulfur
    17: ('Cl', 1.02, (0.12, 0.94, 0.12)), # Chlorine
    18: ('Ar', 1.06, (0.5, 0.82, 0.89)),  # Argon
    26: ('Fe', 1.32, (0.71, 0.45, 0.13)), # Iron
    29: ('Cu', 1.28, (0.72, 0.45, 0.20)), # Copper
    30: ('Zn', 1.22, (0.49, 0.50, 0.69)), # Zinc
    35: ('Br', 1.20, (0.65, 0.16, 0.16)), # Bromine
    53: ('I', 1.39, (0.58, 0.0, 0.58)),   # Iodine
}

# Extended bond length ranges definition
BOND_RANGES = {
    ('C', 'C'): [
        (1.17, 1.25, BondType.TRIPLE, 3),
        (1.29, 1.35, BondType.DOUBLE, 2),
        (1.30, 1.45, BondType.AROMATIC, 1.5),
        (1.45, 1.70, BondType.SINGLE, 1),
    ],
    ('C', 'O'): [
        (1.13, 1.22, BondType.TRIPLE, 3),
        (1.20, 1.31, BondType.DOUBLE, 2),
        (1.31, 1.50, BondType.SINGLE, 1),
    ],
    ('C', 'N'): [
        (1.13, 1.19, BondType.TRIPLE, 3),
        (1.21, 1.30, BondType.DOUBLE, 2),
        (1.30, 1.45, BondType.SINGLE, 1),
    ],
    ('O', 'H'): [(0.85, 1.10, BondType.SINGLE, 1)],
    ('N', 'H'): [(0.95, 1.10, BondType.SINGLE, 1)],
    ('C', 'H'): [(1.06, 1.12, BondType.SINGLE, 1)],
    ('O', 'O'): [
        (1.20, 1.35, BondType.DOUBLE, 2),
        (1.35, 1.60, BondType.SINGLE, 1)
    ],
    ('N', 'N'): [
        (1.10, 1.25, BondType.TRIPLE, 3),
        (1.25, 1.35, BondType.DOUBLE, 2),
        (1.35, 1.45, BondType.SINGLE, 1)
    ],
    ('P', 'O'): [
        (1.35, 1.45, BondType.DOUBLE, 2),
        (1.45, 1.65, BondType.SINGLE, 1)
    ],
    ('S', 'O'): [
        (1.35, 1.45, BondType.DOUBLE, 2),
        (1.45, 1.65, BondType.SINGLE, 1)
    ],
    ('C', 'S'): [
        (1.55, 1.70, BondType.SINGLE, 1),
        (1.70, 1.85, BondType.DOUBLE, 2)
    ],
    ('C', 'P'): [(1.80, 1.85, BondType.SINGLE, 1)],
    ('C', 'Si'): [(1.85, 1.90, BondType.SINGLE, 1)],
    ('Si', 'O'): [(1.55, 1.65, BondType.SINGLE, 1)],
}

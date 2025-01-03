"""
Molecular structure analyzer module.
"""
import numpy as np
from collections import defaultdict
from typing import List, Dict, Tuple, Optional

from .constants import ELEMENTS, BOND_RANGES, BondType

class MoleculeAnalyzer:
    """Molecular structure analyzer class for analyzing molecular properties."""
    
    def __init__(self, atoms: List[Tuple]):
        """
        Initialize the analyzer with atomic coordinates.
        
        Args:
            atoms: List of tuples containing (atomic_number, x, y, z)
        """
        self.atoms = atoms
        self.bonds = []
        self.rings = []
        self.functional_groups = []
        
    def analyze_structure(self) -> Dict:
        """
        Perform complete structural analysis.
        
        Returns:
            Dict containing composition, bonds, rings, functional groups and topology
        """
        return {
            'composition': self.analyze_composition(),
            'bonds': self.analyze_bonds(),
            'rings': self.find_rings(),
            'functional_groups': self.identify_functional_groups(),
            'topology': self.analyze_topology()
        }
        
    def analyze_composition(self) -> Dict[str, int]:
        """
        Analyze atomic composition.
        
        Returns:
            Dict mapping element symbols to their counts
        """
        composition = defaultdict(int)
        for atomic_num, *_ in self.atoms:
            symbol = ELEMENTS[atomic_num][0]
            composition[symbol] += 1
        return dict(composition)
        
    def analyze_bonds(self) -> List[Dict]:
        """
        Analyze chemical bonds.
        
        Returns:
            List of bond information dictionaries
        """
        bonds = []
        for i in range(len(self.atoms)):
            for j in range(i + 1, len(self.atoms)):
                bond_info = self._check_bond(i, j)
                if bond_info:
                    bonds.append(bond_info)
        self.bonds = bonds
        return bonds
        
    def _check_bond(self, i: int, j: int) -> Optional[Dict]:
        """
        Check if two atoms form a bond.
        
        Args:
            i: Index of first atom
            j: Index of second atom
            
        Returns:
            Bond information dictionary if bond exists, None otherwise
        """
        pos1 = np.array(self.atoms[i][1:])
        pos2 = np.array(self.atoms[j][1:])
        dist = np.linalg.norm(pos2 - pos1)
        
        symbol1 = ELEMENTS[self.atoms[i][0]][0]
        symbol2 = ELEMENTS[self.atoms[j][0]][0]
        
        atoms_pair = tuple(sorted([symbol1, symbol2]))
        if atoms_pair in BOND_RANGES or atoms_pair[::-1] in BOND_RANGES:
            ranges = BOND_RANGES.get(atoms_pair) or BOND_RANGES.get(atoms_pair[::-1])
            for min_dist, max_dist, bond_type, bond_order in ranges:
                if min_dist <= dist <= max_dist:
                    return {
                        'atoms': (i, j),
                        'symbols': (symbol1, symbol2),
                        'distance': dist,
                        'type': bond_type,
                        'order': bond_order
                    }
        return None
        
    def find_rings(self) -> List[Dict]:
        """
        Identify ring systems.
        
        Returns:
            List of ring information dictionaries
        """
        def find_cycles(start: int, visited: set, path: List[int]) -> List[List[int]]:
            cycles = []
            neighbors = self._get_neighbors(start)
            
            for next_atom in neighbors:
                if next_atom in path[:-1]:
                    if len(path) >= 3 and next_atom == path[0]:
                        cycles.append(path[:])
                elif next_atom not in visited:
                    visited.add(next_atom)
                    path.append(next_atom)
                    cycles.extend(find_cycles(next_atom, visited, path))
                    path.pop()
                    visited.remove(next_atom)
            return cycles
        
        rings = []
        visited = set()
        for i in range(len(self.atoms)):
            if i not in visited:
                visited.add(i)
                cycles = find_cycles(i, visited, [i])
                for cycle in cycles:
                    if len(cycle) >= 3:  # Only consider 3+ membered rings
                        ring_info = self._analyze_ring(cycle)
                        if ring_info not in rings:
                            rings.append(ring_info)
        
        self.rings = rings
        return rings
        
    def _get_neighbors(self, atom_idx: int) -> List[int]:
        """
        Get indices of all atoms bonded to specified atom.
        
        Args:
            atom_idx: Index of atom to find neighbors for
            
        Returns:
            List of neighboring atom indices
        """
        neighbors = []
        for bond in self.bonds:
            if bond['atoms'][0] == atom_idx:
                neighbors.append(bond['atoms'][1])
            elif bond['atoms'][1] == atom_idx:
                neighbors.append(bond['atoms'][0])
        return neighbors
        
    def _analyze_ring(self, cycle: List[int]) -> Dict:
        """
        Analyze ring characteristics.
        
        Args:
            cycle: List of atom indices forming the ring
            
        Returns:
            Ring information dictionary
        """
        atoms = []
        is_aromatic = True
        total_electrons = 0
        
        for i in cycle:
            atomic_num = self.atoms[i][0]
            symbol = ELEMENTS[atomic_num][0]
            atoms.append(symbol)
            
            # Calculate valence electrons (simplified)
            if symbol in ['C', 'Si']:
                total_electrons += 4
            elif symbol in ['N', 'P']:
                total_electrons += 5
            elif symbol in ['O', 'S']:
                total_electrons += 6
            else:
                is_aromatic = False
        
        # Check Hückel's rule (4n+2)
        if is_aromatic and total_electrons % 4 != 2:
            is_aromatic = False
            
        return {
            'size': len(cycle),
            'atoms': atoms,
            'is_aromatic': is_aromatic,
            'indices': cycle
        }
        
    def identify_functional_groups(self) -> List[str]:
        """
        Identify functional groups.
        
        Returns:
            List of identified functional group names
        """
        groups = []
        
        # Define common functional group patterns
        patterns = {
            'alcohol': self._check_alcohol,
            'carboxyl': self._check_carboxyl,
            'amine': self._check_amine,
            'carbonyl': self._check_carbonyl,
            'ether': self._check_ether
        }
        
        for name, checker in patterns.items():
            if checker():
                groups.append(name)
                
        self.functional_groups = groups
        return groups
        
    def _check_alcohol(self) -> bool:
        """Check for alcohol hydroxyl group."""
        for i, (atomic_num, *_) in enumerate(self.atoms):
            if ELEMENTS[atomic_num][0] == 'O':
                neighbors = self._get_neighbors(i)
                if len(neighbors) == 2:
                    neighbor_symbols = [ELEMENTS[self.atoms[j][0]][0] for j in neighbors]
                    if 'H' in neighbor_symbols and 'C' in neighbor_symbols:
                        return True
        return False
        
    def _check_carboxyl(self) -> bool:
        """Check for carboxyl group."""
        for i, (atomic_num, *_) in enumerate(self.atoms):
            if ELEMENTS[atomic_num][0] == 'C':
                neighbors = self._get_neighbors(i)
                if len(neighbors) == 3:
                    neighbor_symbols = [ELEMENTS[self.atoms[j][0]][0] for j in neighbors]
                    if neighbor_symbols.count('O') == 2:
                        return True
        return False
        
    def _check_amine(self) -> bool:
        """Check for amine group."""
        for i, (atomic_num, *_) in enumerate(self.atoms):
            if ELEMENTS[atomic_num][0] == 'N':
                neighbors = self._get_neighbors(i)
                if len(neighbors) <= 3:
                    neighbor_symbols = [ELEMENTS[self.atoms[j][0]][0] for j in neighbors]
                    if neighbor_symbols.count('H') >= 1:
                        return True
        return False
        
    def _check_carbonyl(self) -> bool:
        """Check for carbonyl group."""
        for i, (atomic_num, *_) in enumerate(self.atoms):
            if ELEMENTS[atomic_num][0] == 'C':
                neighbors = self._get_neighbors(i)
                if len(neighbors) >= 2:
                    neighbor_symbols = [ELEMENTS[self.atoms[j][0]][0] for j in neighbors]
                    bonds = [b for b in self.bonds if i in b['atoms']]
                    for neighbor, bond in zip(neighbors, bonds):
                        if ELEMENTS[self.atoms[neighbor][0]][0] == 'O' and bond['type'] == BondType.DOUBLE:
                            return True
        return False
        
    def _check_ether(self) -> bool:
        """Check for ether group."""
        for i, (atomic_num, *_) in enumerate(self.atoms):
            if ELEMENTS[atomic_num][0] == 'O':
                neighbors = self._get_neighbors(i)
                if len(neighbors) == 2:
                    neighbor_symbols = [ELEMENTS[self.atoms[j][0]][0] for j in neighbors]
                    if neighbor_symbols.count('C') == 2:
                        return True
        return False

    def analyze_topology(self) -> Dict:
        """
        Analyze molecular topology.
        
        Returns:
            Dictionary containing connectivity, branching and ring connectivity information
        """
        return {
            'connectivity': self._analyze_connectivity(),
            'branching': self._analyze_branching(),
            'ring_connectivity': self._analyze_ring_connectivity()
        }
        
    def _analyze_connectivity(self) -> Dict:
        """Analyze atomic connectivity."""
        connectivity = defaultdict(int)
        for bond in self.bonds:
            atom1, atom2 = bond['atoms']
            symbol1 = ELEMENTS[self.atoms[atom1][0]][0]
            symbol2 = ELEMENTS[self.atoms[atom2][0]][0]
            connectivity[symbol1] += 1
            connectivity[symbol2] += 1
        return dict(connectivity)
        
    def _analyze_branching(self) -> Dict:
        """Analyze branching patterns."""
        branching = {
            'terminal': 0,    # Terminal atoms
            'linear': 0,      # Linear chain atoms
            'branched': 0     # Branch points
        }
        
        for i in range(len(self.atoms)):
            neighbors = self._get_neighbors(i)
            if len(neighbors) == 1:
                branching['terminal'] += 1
            elif len(neighbors) == 2:
                branching['linear'] += 1
            elif len(neighbors) > 2:
                branching['branched'] += 1
                
        return branching
        
    def _analyze_ring_connectivity(self) -> Dict:
        """Analyze ring system connectivity."""
        if not self.rings:
            return {'fused_rings': 0, 'spiro_rings': 0, 'isolated_rings': 0}
            
        fused = 0    # Fused rings
        spiro = 0    # Spiro rings
        isolated = 0 # Isolated rings
        
        for i, ring1 in enumerate(self.rings):
            shares_atoms = False
            for j, ring2 in enumerate(self.rings):
                if i != j:
                    common_atoms = set(ring1['indices']) & set(ring2['indices'])
                    if len(common_atoms) == 2:
                        fused += 1
                    elif len(common_atoms) == 1:
                        spiro += 1
                    elif len(common_atoms) == 0:
                        isolated += 1
                        
        return {
            'fused_rings': fused // 2,  # Divide by 2 to avoid double counting
            'spiro_rings': spiro // 2,
            'isolated_rings': max(0, len(self.rings) - (fused + spiro) // 2)
        }

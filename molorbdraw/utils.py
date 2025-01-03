"""
Utility functions for molecular visualization and analysis.
"""
from typing import Tuple

from .constants import ELEMENTS

def get_element_info(atomic_number: int) -> Tuple[str, float, Tuple[float, float, float]]:
    """
    Get element information from atomic number.
    
    Args:
        atomic_number: The atomic number of the element
        
    Returns:
        Tuple containing:
            - str: Element symbol
            - float: Atomic radius
            - Tuple[float, float, float]: RGB color values
            
    Example:
        >>> symbol, radius, color = get_element_info(6)
        >>> print(symbol)  # 'C'
        >>> print(radius)  # 0.76
        >>> print(color)   # (0.5, 0.5, 0.5)
    """
    return ELEMENTS.get(atomic_number, ('X', 0.5, (0.8, 0.8, 0.8)))

def calculate_bond_midpoint(pos1: Tuple[float, float, float], 
                          pos2: Tuple[float, float, float]) -> Tuple[float, float, float]:
    """
    Calculate the midpoint between two atomic positions.
    
    Args:
        pos1: (x, y, z) coordinates of first atom
        pos2: (x, y, z) coordinates of second atom
        
    Returns:
        Tuple containing the (x, y, z) coordinates of the midpoint
    """
    return tuple((p1 + p2) / 2 for p1, p2 in zip(pos1, pos2))

def validate_cube_file(filename: str) -> bool:
    """
    Validate if a file is in Gaussian cube format.
    
    Args:
        filename: Path to the cube file
        
    Returns:
        bool: True if file appears to be valid cube format, False otherwise
        
    Note:
        This is a basic validation that checks file structure but not content quality.
    """
    try:
        with open(filename, 'r') as f:
            # Skip first two comment lines
            f.readline()
            f.readline()
            
            # Read number of atoms and origin
            parts = f.readline().split()
            if len(parts) < 4:
                return False
                
            try:
                n_atoms = abs(int(float(parts[0])))
            except ValueError:
                return False
                
            # Check for three grid vectors
            for _ in range(3):
                parts = f.readline().split()
                if len(parts) < 4:
                    return False
                    
            # Check for atom lines
            for _ in range(n_atoms):
                parts = f.readline().split()
                if len(parts) < 5:
                    return False
                    
            return True
            
    except Exception:
        return False

def format_scientific(value: float, precision: int = 6) -> str:
    """
    Format a number in scientific notation.
    
    Args:
        value: Number to format
        precision: Number of decimal places (default: 6)
        
    Returns:
        str: Formatted string in scientific notation
        
    Example:
        >>> format_scientific(0.00123456)
        '1.234560E-03'
    """
    return f"{value:.{precision}E}"

def estimate_memory_usage(nx: int, ny: int, nz: int) -> float:
    """
    Estimate memory usage for a grid of given dimensions.
    
    Args:
        nx: Number of points in x direction
        ny: Number of points in y direction
        nz: Number of points in z direction
        
    Returns:
        float: Estimated memory usage in megabytes
        
    Note:
        Assumes 8 bytes per grid point (double precision)
    """
    bytes_per_point = 8  # double precision float
    total_bytes = nx * ny * nz * bytes_per_point
    return total_bytes / (1024 * 1024)  # Convert to MB

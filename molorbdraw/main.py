"""
Main entry point for molecular orbital visualization.
"""
import os
import sys
from typing import Optional

from .viewer import MoleculeViewer
from .utils import validate_cube_file, estimate_memory_usage

def view_cube_file(filename: str) -> None:
    """
    Load and display a cube file visualization.
    
    Args:
        filename: Path to the cube file to visualize
    """
    if not os.path.exists(filename):
        print(f"Error: File not found: {filename}")
        sys.exit(1)
        
    if not validate_cube_file(filename):
        print(f"Error: Invalid cube file format: {filename}")
        sys.exit(1)
        
    viewer = MoleculeViewer()
    
    try:
        # Set up visualization environment
        viewer.setup_visualization()
        
        # Set up controls (must be done before loading data)
        viewer.setup_controls()
        
        # Load data and create visualization
        print(f"Reading file: {filename}")
        viewer.load_cube_file(filename)
        viewer.create_visualization()
        
        # Adjust camera position
        viewer.renderer.ResetCamera()
        camera = viewer.renderer.GetActiveCamera()
        camera.Elevation(45)
        camera.Azimuth(45)
        
        # Start interaction
        viewer.render_window.Render()
        viewer.interactor.Initialize()
        viewer.interactor.Start()
        
    except Exception as e:
        print(f"Error: {str(e)}")
        raise

def main() -> None:
    """Main entry point."""
    if len(sys.argv) != 2:
        print("Usage: python -m molorbdraw <cube file path>")
        sys.exit(1)
    
    try:
        cube_file = sys.argv[1]
        print(f"Opening file: {cube_file}")
        view_cube_file(cube_file)
    except Exception as e:
        print(f"Error occurred: {str(e)}")
        sys.exit(1)

if __name__ == "__main__":
    main()

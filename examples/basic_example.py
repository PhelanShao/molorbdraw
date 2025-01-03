"""
Basic example demonstrating MolOrbDraw usage.
"""
import os
from molorbdraw import MoleculeViewer

def main():
    # Get the path to the example cube file
    current_dir = os.path.dirname(os.path.abspath(__file__))
    cube_file = os.path.join(current_dir, "DD.cub")
    
    if not os.path.exists(cube_file):
        print(f"Error: Example file not found: {cube_file}")
        return
    
    print("Creating molecular orbital visualization...")
    
    # Create viewer instance
    viewer = MoleculeViewer()
    
    # Set up visualization environment
    viewer.setup_visualization()
    
    # Set up interactive controls
    viewer.setup_controls()
    
    # Load and display the cube file
    print(f"Loading cube file: {cube_file}")
    viewer.load_cube_file(cube_file)
    viewer.create_visualization()
    
    # Adjust initial camera position
    viewer.renderer.ResetCamera()
    camera = viewer.renderer.GetActiveCamera()
    camera.Elevation(45)
    camera.Azimuth(45)
    
    print("\nVisualization Controls:")
    print("- Left Mouse: Rotate")
    print("- Middle Mouse: Pan")
    print("- Right Mouse: Zoom")
    print("- 'A': Toggle atoms")
    print("- 'B': Toggle bonds")
    print("- 'I': Toggle isosurface")
    print("- '1/2/3': Change color schemes")
    print("\nUse control panel sliders to adjust:")
    print("- Isosurface value")
    print("- Opacity")
    print("- Bond thickness")
    
    # Start interaction
    viewer.render_window.Render()
    viewer.interactor.Initialize()
    viewer.interactor.Start()

if __name__ == "__main__":
    main()

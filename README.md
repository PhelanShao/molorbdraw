# MolOrbDraw

MolOrbDraw is a powerful visualization tool for molecular orbitals from Gaussian cube files. It provides an interactive 3D visualization interface with advanced molecular structure analysis capabilities.

## Features

- **Interactive 3D Visualization**: View molecular orbitals with real-time adjustable isosurfaces
- **Molecular Structure Analysis**: 
  - Automatic bond detection and visualization
  - Ring system identification
  - Functional group recognition
  - Topology analysis
- **Customizable Display**:
  - Adjustable isosurface values
  - Multiple color schemes
  - Controllable opacity
  - Adjustable bond thickness
- **Comprehensive Chemical Information**:
  - Support for extended periodic table
  - Accurate bond length ranges
  - Multiple bond type visualization (single, double, triple, aromatic)
- **User-Friendly Interface**:
  - Interactive control panel
  - Keyboard shortcuts for common operations
  - Real-time molecular information display

## Installation

1. Ensure you have Python 3.7 or newer installed
2. Install the package using pip:
```bash
pip install molorbdraw
```

## Dependencies

- Python >= 3.7
- NumPy >= 1.19.0
- VTK >= 9.0.0

## Usage

### Command Line Interface

```bash
# View a cube file
molorbdraw path/to/your/file.cube
```
![Description](https://github.com/PhelanShao/molorbdraw/raw/main/examples/example.png)

### Python API

```python
from molorbdraw import MoleculeViewer

# Create viewer instance
viewer = MoleculeViewer()

# Load and display cube file
viewer.load_cube_file("path/to/your/file.cube")
viewer.setup_visualization()
viewer.setup_controls()
viewer.create_visualization()

# Start interaction
viewer.render_window.Render()
viewer.interactor.Start()
```

## Controls

### Keyboard Shortcuts

- `A`: Toggle atom visibility
- `B`: Toggle bond visibility
- `I`: Toggle isosurface visibility
- `1`: Red-Blue color scheme
- `2`: Green-Purple color scheme
- `3`: Yellow-Cyan color scheme

### Mouse Controls

- Left Mouse Button: Rotate camera
- Middle Mouse Button: Pan camera
- Right Mouse Button: Zoom camera
- Mouse Wheel: Zoom in/out

### Control Panel

The control panel provides sliders for:
- Isosurface value adjustment
- Opacity control
- Bond thickness adjustment

## File Format Support

Currently supports Gaussian cube files (.cube) with the following specifications:
- Standard cube file format
- Contains atomic coordinates and volumetric data
- Supports both molecular orbitals and electron density data

## Contributing

Contributions are welcome! Please feel free to submit a Pull Request.

## License

This project is licensed under the MIT License - see the LICENSE file for details.

## Acknowledgments

This project builds upon the work of various molecular visualization tools and chemical informatics libraries. Special thanks to:
- VTK (Visualization Toolkit) for 3D rendering capabilities
- The computational chemistry community for cube file format specifications

## Citation

If you use MolOrbDraw in your research, please cite:

```bibtex
@software{molorbdraw2023,
  author = {MolOrbDraw Contributors},
  title = {MolOrbDraw: A Molecular Orbital Visualization Tool},
  year = {2023},
  version = {0.1.0}
}

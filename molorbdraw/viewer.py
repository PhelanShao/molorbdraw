"""
Molecular visualization module using VTK.
"""
import vtk
import numpy as np
from typing import List, Tuple, Optional

from .constants import BOHR_TO_ANGSTROM, BondType
from .analyzer import MoleculeAnalyzer
from .utils import get_element_info

class MoleculeViewer:
    """3D molecular visualization class using VTK."""
    
    def __init__(self):
        """Initialize the viewer."""
        # VTK objects
        self.renderer = None
        self.render_window = None
        self.interactor = None
        self.panel_renderer = None
        
        # Molecule data
        self.data = None
        self.atoms = None
        self.bonds = None
        self.analyzer = None
        self.grid = None
        
        # Visualization objects
        self.isosurface_actor = None
        self.atom_actors = []
        self.bond_actors = []
        self.info_text = None
        self.contour = None
        self.mapper = None
        
        # Current display settings
        self.iso_value = 0.02
        self.opacity = 0.7
        self.bond_radius = 0.1

    def load_cube_file(self, filename: str) -> None:
        """
        Load molecular data from a cube file.
        
        Args:
            filename: Path to the cube file
        """
        with open(filename, 'r') as f:
            lines = f.readlines()
        
        current_line = 2  # Skip comments
        
        # Read number of atoms and origin
        parts = lines[current_line].split()
        n_atoms = abs(int(float(parts[0])))
        origin = np.array([float(x) for x in parts[1:4]])
        current_line += 1
        
        # Read grid information
        grid_info = []
        for i in range(3):
            parts = lines[current_line + i].split()
            n = int(parts[0])
            d = float(parts[i + 1])
            grid_info.append((n, d))
        current_line += 3
        
        nx, dx = grid_info[0]
        ny, dy = grid_info[1]
        nz, dz = grid_info[2]
        
        # Read atom information
        self.atoms = []
        for i in range(n_atoms):
            parts = lines[current_line + i].split()
            try:
                atom_num = int(float(parts[0]))
                pos = [float(x) * BOHR_TO_ANGSTROM for x in parts[2:5]]
                self.atoms.append((atom_num, *pos))
            except ValueError:
                print(f"Warning: Error parsing atom line {i + 1}")
                continue
        
        current_line += n_atoms
        
        # Create molecule analyzer
        self.analyzer = MoleculeAnalyzer(self.atoms)
        
        # Read volume data
        tmp_data = []
        for line in lines[current_line:]:
            try:
                tmp_data.extend([float(x) for x in line.split()])
            except ValueError:
                continue
                
        # Convert data and grid
        self.data = np.array(tmp_data).reshape((nx, ny, nz))
        self.grid = (
            (origin[0] + np.arange(nx) * dx) * BOHR_TO_ANGSTROM,
            (origin[1] + np.arange(ny) * dy) * BOHR_TO_ANGSTROM,
            (origin[2] + np.arange(nz) * dz) * BOHR_TO_ANGSTROM
        )

    def setup_visualization(self) -> None:
        """Set up VTK visualization environment."""
        # Create main window
        self.render_window = vtk.vtkRenderWindow()
        self.render_window.SetSize(1200, 800)
        
        # Create renderer and set viewport
        self.renderer = vtk.vtkRenderer()
        self.renderer.SetBackground(1, 1, 1)  # White background
        self.renderer.SetViewport(0.2, 0, 1.0, 1.0)  # Reserve left space for control panel
        
        # Create control panel renderer
        self.panel_renderer = vtk.vtkRenderer()
        self.panel_renderer.SetBackground(1.0, 0.9, 0.9)  # Light red background
        self.panel_renderer.SetViewport(0, 0, 0.2, 1.0)
        
        self.render_window.AddRenderer(self.panel_renderer)
        self.render_window.AddRenderer(self.renderer)
        
        self.interactor = vtk.vtkRenderWindowInteractor()
        self.interactor.SetRenderWindow(self.render_window)
        
        style = vtk.vtkInteractorStyleTrackballCamera()
        self.interactor.SetInteractorStyle(style)

    def create_visualization(self) -> None:
        """Create the complete visualization."""
        # Create VTK image data
        image_data = vtk.vtkImageData()
        image_data.SetDimensions(self.data.shape)
        image_data.SetOrigin(self.grid[0][0], self.grid[1][0], self.grid[2][0])
        image_data.SetSpacing(
            (self.grid[0][-1] - self.grid[0][0]) / (len(self.grid[0]) - 1),
            (self.grid[1][-1] - self.grid[1][0]) / (len(self.grid[1]) - 1),
            (self.grid[2][-1] - self.grid[2][0]) / (len(self.grid[2]) - 1)
        )
        
        # Add data
        scalar_array = vtk.vtkDoubleArray()
        scalar_array.SetNumberOfComponents(1)
        for value in self.data.flatten(order='F'):
            scalar_array.InsertNextValue(value)
        image_data.GetPointData().SetScalars(scalar_array)
        
        # Create isosurface
        self.contour = vtk.vtkContourFilter()
        self.contour.SetInputData(image_data)
        self.contour.SetValue(0, self.iso_value)
        self.contour.SetValue(1, -self.iso_value)
        
        # Create mapper and actor
        self.mapper = vtk.vtkPolyDataMapper()
        self.mapper.SetInputConnection(self.contour.GetOutputPort())
        self.mapper.ScalarVisibilityOn()
        
        self.isosurface_actor = vtk.vtkActor()
        self.isosurface_actor.SetMapper(self.mapper)
        self.isosurface_actor.GetProperty().SetOpacity(self.opacity)
        
        self.renderer.AddActor(self.isosurface_actor)
        
        # Add atoms and bonds
        self.add_atoms_and_bonds()
        
        # Update information display
        self.update_info_display()

    def add_atoms_and_bonds(self) -> None:
        """Add atoms and chemical bonds visualization."""
        if not self.atoms:
            return
            
        # Analyze molecular structure
        analysis = self.analyzer.analyze_structure()
        self.bonds = analysis['bonds']
            
        # Add atoms
        for atomic_num, x, y, z in self.atoms:
            symbol, radius, color = get_element_info(atomic_num)
            
            # Create sphere for atom
            sphere = vtk.vtkSphereSource()
            sphere.SetCenter(x, y, z)
            sphere.SetRadius(radius * 0.4)
            sphere.SetPhiResolution(20)
            sphere.SetThetaResolution(20)
            
            mapper = vtk.vtkPolyDataMapper()
            mapper.SetInputConnection(sphere.GetOutputPort())
            
            actor = vtk.vtkActor()
            actor.SetMapper(mapper)
            actor.GetProperty().SetColor(color)
            
            # Add label
            label = vtk.vtkVectorText()
            label.SetText(symbol)
            
            label_mapper = vtk.vtkPolyDataMapper()
            label_mapper.SetInputConnection(label.GetOutputPort())
            
            label_actor = vtk.vtkFollower()
            label_actor.SetMapper(label_mapper)
            label_actor.SetScale(0.3, 0.3, 0.3)
            label_actor.SetPosition(x + radius * 0.5, y + radius * 0.5, z + radius * 0.5)
            label_actor.SetCamera(self.renderer.GetActiveCamera())
            
            self.atom_actors.extend([actor, label_actor])
            self.renderer.AddActor(actor)
            self.renderer.AddActor(label_actor)
        
        # Add bonds
        for bond in self.bonds:
            i, j = bond['atoms']
            pos1 = np.array(self.atoms[i][1:])
            pos2 = np.array(self.atoms[j][1:])
            self.add_bond(pos1, pos2, bond['type'])

    def add_bond(self, pos1: np.ndarray, pos2: np.ndarray, bond_type: BondType) -> None:
        """
        Add a chemical bond visualization.
        
        Args:
            pos1: Position of first atom
            pos2: Position of second atom
            bond_type: Type of chemical bond
        """
        if bond_type == BondType.SINGLE:
            actors = self.create_single_bond(pos1, pos2)
        elif bond_type == BondType.DOUBLE:
            actors = self.create_double_bond(pos1, pos2)
        elif bond_type == BondType.TRIPLE:
            actors = self.create_triple_bond(pos1, pos2)
        elif bond_type == BondType.AROMATIC:
            actors = self.create_aromatic_bond(pos1, pos2)
        else:
            actors = self.create_single_bond(pos1, pos2)
        
        for actor in actors:
            self.bond_actors.append(actor)
            self.renderer.AddActor(actor)

    def create_single_bond(self, pos1: np.ndarray, pos2: np.ndarray, 
                         radius: Optional[float] = None) -> List[vtk.vtkActor]:
        """Create visualization for a single bond."""
        if radius is None:
            radius = self.bond_radius
            
        line = vtk.vtkLineSource()
        line.SetPoint1(*pos1)
        line.SetPoint2(*pos2)
        
        tube = vtk.vtkTubeFilter()
        tube.SetInputConnection(line.GetOutputPort())
        tube.SetRadius(radius)
        tube.SetNumberOfSides(20)
        
        mapper = vtk.vtkPolyDataMapper()
        mapper.SetInputConnection(tube.GetOutputPort())
        
        actor = vtk.vtkActor()
        actor.SetMapper(mapper)
        actor.GetProperty().SetColor(0.8, 0.8, 0.8)
        actor.tube = tube  # Save reference for later updates
        
        return [actor]

    def create_double_bond(self, pos1: np.ndarray, pos2: np.ndarray) -> List[vtk.vtkActor]:
        """Create visualization for a double bond."""
        # Calculate bond direction vector
        direction = pos2 - pos1
        direction = direction / np.linalg.norm(direction)
        
        # Calculate perpendicular vector
        if abs(direction[0]) < abs(direction[1]):
            perpendicular = np.cross(direction, [1, 0, 0])
        else:
            perpendicular = np.cross(direction, [0, 1, 0])
        perpendicular = perpendicular / np.linalg.norm(perpendicular)
        
        # Create two parallel bonds
        offset = perpendicular * self.bond_radius * 2
        pos1_1 = pos1 + offset
        pos2_1 = pos2 + offset
        pos1_2 = pos1 - offset
        pos2_2 = pos2 - offset
        
        return (self.create_single_bond(pos1_1, pos2_1, self.bond_radius*0.7) +
                self.create_single_bond(pos1_2, pos2_2, self.bond_radius*0.7))

    def create_triple_bond(self, pos1: np.ndarray, pos2: np.ndarray) -> List[vtk.vtkActor]:
        """Create visualization for a triple bond."""
        # Middle bond
        actors = self.create_single_bond(pos1, pos2)
        
        # Calculate perpendicular offset
        direction = pos2 - pos1
        direction = direction / np.linalg.norm(direction)
        
        if abs(direction[0]) < abs(direction[1]):
            perpendicular = np.cross(direction, [1, 0, 0])
        else:
            perpendicular = np.cross(direction, [0, 1, 0])
        perpendicular = perpendicular / np.linalg.norm(perpendicular)
        
        # Add two additional bonds
        offset = perpendicular * self.bond_radius * 2.5
        pos1_1 = pos1 + offset
        pos2_1 = pos2 + offset
        pos1_2 = pos1 - offset
        pos2_2 = pos2 - offset
        
        actors.extend(self.create_single_bond(pos1_1, pos2_1, self.bond_radius*0.7))
        actors.extend(self.create_single_bond(pos1_2, pos2_2, self.bond_radius*0.7))
        
        return actors

    def create_aromatic_bond(self, pos1: np.ndarray, pos2: np.ndarray) -> List[vtk.vtkActor]:
        """Create visualization for an aromatic bond."""
        # Solid line part
        actors = self.create_single_bond(pos1, pos2)
        
        # Calculate dashed line position
        direction = pos2 - pos1
        direction = direction / np.linalg.norm(direction)
        
        if abs(direction[0]) < abs(direction[1]):
            perpendicular = np.cross(direction, [1, 0, 0])
        else:
            perpendicular = np.cross(direction, [0, 1, 0])
        perpendicular = perpendicular / np.linalg.norm(perpendicular)
        
        # Add dashed line
        offset = perpendicular * self.bond_radius * 2
        pos1_1 = pos1 + offset
        pos2_1 = pos2 + offset
        
        # Create dashed effect
        num_segments = 10
        segment_vector = (pos2_1 - pos1_1) / num_segments
        
        for i in range(0, num_segments, 2):
            start = pos1_1 + segment_vector * i
            end = pos1_1 + segment_vector * (i + 1)
            actors.extend(self.create_single_bond(start, end, self.bond_radius*0.7))
        
        return actors

    def setup_controls(self) -> None:
        """Set up control panel interface."""
        # Title
        self.create_control_text("Parameter Controls", 0.02, 0.95)
        
        # Color scheme
        self.create_control_text("Color Scheme:", 0.02, 0.90)
        self.create_control_text("1: Red-Blue", 0.03, 0.87)
        self.create_control_text("2: Green-Purple", 0.03, 0.84)
        self.create_control_text("3: Yellow-Cyan", 0.03, 0.81)
        
        # Isosurface control
        self.iso_slider = self.create_slider("Isosurface Value", (0.001, 0.1), 0.02, (0.02, 0.75))
        self.iso_slider.AddObserver("InteractionEvent", self.update_isosurface)
        
        # Opacity control
        self.opacity_slider = self.create_slider("Opacity", (0.0, 1.0), 0.7, (0.02, 0.65))
        self.opacity_slider.AddObserver("InteractionEvent", self.update_opacity)
        
        # Bond thickness control
        self.bond_slider = self.create_slider("Bond Thickness", (0.05, 0.3), 0.1, (0.02, 0.55))
        self.bond_slider.AddObserver("InteractionEvent", self.update_bond_thickness)
        
        # Display controls
        self.create_control_text("Display Controls:", 0.02, 0.45)
        self.create_control_text("A: Show/Hide Atoms", 0.03, 0.42)
        self.create_control_text("B: Show/Hide Bonds", 0.03, 0.39)
        self.create_control_text("I: Show/Hide Isosurface", 0.03, 0.36)
        
        # Molecule information display area
        self.info_text = self.create_control_text("", 0.02, 0.25)
        
        # Add keyboard event handler
        self.interactor.AddObserver("KeyPressEvent", self.handle_keypress)

    def create_slider(self, title: str, value_range: Tuple[float, float], 
                     initial_value: float, position: Tuple[float, float]) -> vtk.vtkSliderWidget:
        """Create a slider control widget."""
        slider = vtk.vtkSliderWidget()
        slider_rep = vtk.vtkSliderRepresentation2D()
        
        slider_rep.SetMinimumValue(value_range[0])
        slider_rep.SetMaximumValue(value_range[1])
        slider_rep.SetValue(initial_value)
        slider_rep.SetTitleText(title)
        
        slider_rep.GetPoint1Coordinate().SetCoordinateSystemToNormalizedDisplay()
        slider_rep.GetPoint1Coordinate().SetValue(position[0], position[1])
        slider_rep.GetPoint2Coordinate().SetCoordinateSystemToNormalizedDisplay()
        slider_rep.GetPoint2Coordinate().SetValue(position[0] + 0.15, position[1])
        
        slider_rep.GetSliderProperty().SetColor(0.1, 0.1, 0.1)
        slider_rep.GetTitleProperty().SetColor(0, 0, 0)
        slider_rep.GetLabelProperty().SetColor(0, 0, 0)
        slider_rep.GetSelectedProperty().SetColor(0.2, 0.2, 0.8)
        
        slider.SetInteractor(self.interactor)
        slider.SetRepresentation(slider_rep)
        slider.EnabledOn()
        
        return slider

    def create_control_text(self, text: str, x: float, y: float) -> vtk.vtkTextActor:
        """Create control panel text."""
        text_actor = vtk.vtkTextActor()
        text_actor.SetInput(text)
        text_prop = text_actor.GetTextProperty()
        text_prop.SetColor(0, 0, 0)
        text_prop.SetFontSize(12)
        text_prop.SetJustificationToLeft()
        
        text_actor.GetPositionCoordinate().SetCoordinateSystemToNormalizedDisplay()
        text_actor.GetPositionCoordinate().SetValue(x, y)
        
        self.panel_renderer.AddActor2D(text_actor)
        return text_actor

    def update_info_display(self) -> None:
        """Update molecule information display."""
        if self.analyzer:
            analysis = self.analyzer.analyze_structure()
            info = f"Molecule Info:\n"
            info += f"Composition: {analysis['composition']}\n"
            info += f"Total Bonds: {len(analysis['bonds'])}\n"
            info += f"Rings: {len(analysis['rings'])}\n"
            if analysis['rings']:
                info += f"Aromatic: {sum(1 for r in analysis['rings'] if r['is_aromatic'])}\n"
            info += f"Groups: {', '.join(analysis['functional_groups'])}"
            
            self.info_text.SetInput(info)

    def update_isosurface(self, obj: vtk.vtkObject, event: str) -> None:
        """Update isosurface value."""
        value = obj.GetRepresentation().GetValue()
        if hasattr(self, 'contour'):
            self.contour.SetValue(0, value)
            self.contour.SetValue(1, -value)
            self.render_window.Render()

    def update_opacity(self, obj: vtk.vtkObject, event: str) -> None:
        """Update opacity value."""
        value = obj.GetRepresentation().GetValue()
        if self.isosurface_actor:
            self.isosurface_actor.GetProperty().SetOpacity(value)
            self.render_window.Render()

    def update_bond_thickness(self, obj: vtk.vtkObject, event: str) -> None:
        """Update bond thickness."""
        value = obj.GetRepresentation().GetValue()
        for actor in self.bond_actors:
            if hasattr(actor, 'tube'):
                actor.tube.SetRadius(value)
        self.render_window.Render()

    def handle_keypress(self, obj: vtk.vtkObject, event: str) -> None:
        """Handle keyboard events."""
        key = obj.GetKeySym().upper()
        
        if key in ['1', '2', '3']:
            colors = {
                '1': ((1,0,0), (0,0,1)),   # Red-Blue
                '2': ((0,1,0), (1,0,1)),   # Green-Purple
                '3': ((1,1,0), (0,1,1))    # Yellow-Cyan
            }
            self.update_isosurface_colors(*colors[key])
        elif key == 'A':
            self._toggle_actors(self.atom_actors)
        elif key == 'B':
            self._toggle_actors(self.bond_actors)
        elif key == 'I':
            if self.isosurface_actor:
                self.isosurface_actor.SetVisibility(
                    not self.isosurface_actor.GetVisibility()
                )
        
        self.render_window.Render()

    def _toggle_actors(self, actors: List[vtk.vtkActor]) -> None:
        """Toggle visibility of a group of actors."""
        if actors:
            visible = not actors[0].GetVisibility()
            for actor in actors:
                actor.SetVisibility(visible)

    def update_isosurface_colors(self, pos_color: Tuple[float, float, float], 
                                neg_color: Tuple[float, float, float]) -> None:
        """Update isosurface colors."""
        if hasattr(self, 'mapper'):
            lut = vtk.vtkLookupTable()
            lut.SetNumberOfColors(2)
            lut.Build()
            lut.SetTableValue(0, *neg_color, 1.0)
            lut.SetTableValue(1, *pos_color, 1.0)
            self.mapper.SetLookupTable(lut)
            self.mapper.SetScalarRange(-self.iso_value, self.iso_value)
            self.render_window.Render()

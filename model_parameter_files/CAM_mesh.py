#!/usr/bin/env python3
"""Generate a tri or quad mesh of a subduction zone vertical profile using Gmsh, making
use of the built-in geometry engine.

Points have been projected from longitude/latitude into a local
transverse Mercator projection. PyLith uses the Proj.4 library
for geographic projections. The proj parameters are:

+proj=tmerc +datum=WGS84 +lon_0=142.0 +lat_0=38.0 +k=0.9996

so that the local origin is at a longitude of 142.0 degrees (WGS84)
and a latitude of 38.0 degrees (WGS84).

Run `generate_gmsh.py --help` to see the command line options.
"""
import gmsh
import numpy as np
from pylith.meshio.gmsh_utils import (VertexGroup, MaterialGroup, GenerateMesh, group_exclude)

class App(GenerateMesh):
    """
    Application for generating the mesh.
    """
    # The elastic thickness of the CAM Plate is 13.7 +/- 1.4km (Hunter & Watts, 2016).
    # The convergence rate of the CAM plate is ~7-8 cm/yr
    # X_WEST = -15.0e+3
    X_WEST = -110.0e+3
    # X_EAST = 20.0e+3
    X_EAST = 60.0e+3
    Y_BOT = -15.0e+3
    Y_TOP = 0.0e+3

    # Top of slab from Slab 2.0 (Hayes et al., 2012, https://doi.org/10.1029/2011JB008524) 
    # OUTER_RISE_FAULT_SURFACES_X = 115e3 - np.array([25, 20.5, 19, 17, 14, 11.8, 9.7, 8.5, 7, 5.5, 4, 3, 2, 1]) * 1e3
    OUTER_RISE_FAULT_SURFACES_X = 0.0 - np.array([78, 71, 67, 62, 58.5, 55.5, 49, 46, 40.5, 39, 36, \
                                                  32, 29, 27, 25.5, 22, 19.5, 15, 12.5, 9.5, 7, 3.5]) * 1e3
    OUTER_RISE_FAULT_SURFACES_Y = np.zeros(len(OUTER_RISE_FAULT_SURFACES_X))
    surface_points = np.array([OUTER_RISE_FAULT_SURFACES_X, OUTER_RISE_FAULT_SURFACES_Y]).T
    fault_dips = np.ones(len(OUTER_RISE_FAULT_SURFACES_X)) * np.deg2rad(60)
    OUTER_RISE_FAULT_BURIED_EDGES_Y = np.ones(len(OUTER_RISE_FAULT_SURFACES_X)) * Y_BOT / 2
    OUTER_RISE_FAULT_BURIED_EDGES_X = surface_points[:, 0] - OUTER_RISE_FAULT_BURIED_EDGES_Y / np.tan(fault_dips)

    def __init__(self):
        """Constructor.
        """
        # Set the cell choices available through command line options
        # with the default cell type `tri` matching the PyLith parameter files.
        self.cell_choices = {
            "default": "tri",
            "choices": ["tri"],
            }
        self.filename = "500m_fault_mesh.msh"

    def create_geometry(self):
        """Create geometry.
        """
        BOTLEFT_PT = gmsh.model.geo.add_point(self.X_WEST, self.Y_BOT, 0.0)
        BOTRIGHT_PT = gmsh.model.geo.add_point(self.X_EAST, self.Y_BOT, 0.0)
        TOPRIGHT_PT = gmsh.model.geo.add_point(self.X_EAST, self.Y_TOP, 0.0)
        TOPLEFT_PT = gmsh.model.geo.add_point(self.X_WEST, self.Y_TOP, 0.0)

        # Create Domain Boundary Curves
        self.c_leftbound = gmsh.model.geo.add_line(TOPLEFT_PT, BOTLEFT_PT)
        self.c_botbound = gmsh.model.geo.add_line(BOTLEFT_PT, BOTRIGHT_PT)
        self.c_rightbound = gmsh.model.geo.add_line(BOTRIGHT_PT, TOPRIGHT_PT)
        
        # Create a fault in the outer rise
        # Perhaps I need to loop through all the points in outer_rise_fault_surface_x, and then
        # I would need to add the corresponding buried edge, create the curve for the outer rise,
        # and then split the slab top to the left and right of the outer rise fault. This split
        # has to be done in the correct order, I will be adding the faults based on increasing 
        # distance from the trench, and so the first time I split the curve I will need to take
        # the seaward side of the split, and then split that curve with the next outer rise fault,
        # propagating this split forward until all outer rise faults have been added.

        self.c_outer_rise_faults = np.zeros(len(self.OUTER_RISE_FAULT_SURFACES_X), dtype=int)
        self.outer_rise_fault_surface_points = np.zeros(len(self.OUTER_RISE_FAULT_SURFACES_X), dtype=int)
        self.outer_rise_fault_buried_edges = np.zeros(len(self.OUTER_RISE_FAULT_BURIED_EDGES_X), dtype=int)

        for i in range(len(self.OUTER_RISE_FAULT_SURFACES_X)):
            self.outer_rise_fault_buried_edges[i] = gmsh.model.geo.add_point(self.OUTER_RISE_FAULT_BURIED_EDGES_X[i], \
                                                                             self.OUTER_RISE_FAULT_BURIED_EDGES_Y[i], \
                                                                             0.0)

            self.outer_rise_fault_surface_points[i] = gmsh.model.geo.add_point(self.OUTER_RISE_FAULT_SURFACES_X[i], \
                                                                               self.OUTER_RISE_FAULT_SURFACES_Y[i], \
                                                                               0.0)

            self.c_outer_rise_faults[i] = gmsh.model.geo.add_polyline([self.outer_rise_fault_surface_points[i], self.outer_rise_fault_buried_edges[i]])
        
        # Slice the Slab top with the first outer rise fault to start the 'cascading curve splitting' 
        # for every other outer rise fault seaward of the trench.
        self.subducting_top_bound_points = TOPRIGHT_PT
        self.subducting_top_bound_points = np.insert(self.subducting_top_bound_points, 0, self.outer_rise_fault_surface_points)
        self.subducting_top_bound_points = np.insert(self.subducting_top_bound_points, 0, TOPLEFT_PT)

        self.c_topbound = gmsh.model.geo.add_polyline(self.subducting_top_bound_points)

        self.all_curves = np.zeros(len(self.outer_rise_fault_surface_points) + 1, dtype=int)
        curves = gmsh.model.geo.split_curve(self.c_topbound, [self.outer_rise_fault_surface_points[0]])
        self.all_curves[0] = curves[0]
        self.all_curves[1] = curves[1]
        for i in range(1, len(self.outer_rise_fault_surface_points)):
            curves = gmsh.model.geo.split_curve(self.all_curves[i], [self.outer_rise_fault_surface_points[i]])
            self.all_curves[i] = curves[0]
            self.all_curves[i + 1] = curves[1]

        # Create surfaces from bounding curves
        # In order to automate this stage of the outer rise implementation I think
        # the best way to do it is to create subarrays the make out the slab porition
        # to the left and right of each outer rise fault + the outer rise fault and then
        # append this to an array that gets fed into add_curve_loop.
        loop_array = np.array([self.c_leftbound,
                               self.c_botbound,
                               self.c_rightbound])

        for j in range(len(self.outer_rise_fault_surface_points)):
            loop_array = np.append(loop_array, -self.all_curves[j])
            loop_array = np.append(loop_array, self.c_outer_rise_faults[j])
            loop_array = np.append(loop_array, -self.c_outer_rise_faults[j])
        loop_array = np.append(loop_array, -self.all_curves[-1])

        loop = gmsh.model.geo.add_curve_loop(loop_array)
        self.s_slab = gmsh.model.geo.add_plane_surface([loop])
        print(self.all_curves)
        gmsh.model.geo.synchronize()


    def mark(self):
        """Mark geometry for materials, boundary conditions, faults, etc.

        This method is abstract in the base class and must be implemented.
        """
        # Create materials matching surfaces.
        materials = (
            MaterialGroup(tag=1, entities=[self.s_slab]),
        )
        for material in materials:
            material.create_physical_group()

        top_boundary_entities = self.all_curves

        # Create physical groups for the boundaries and the fault.
        vertex_groups = (
            VertexGroup(name="bndry_top", tag=10, dim=1, entities=top_boundary_entities),
            VertexGroup(name="bndry_west", tag=11, dim=1, entities=[self.c_leftbound]),
            VertexGroup(name="bndry_east", tag=12, dim=1, entities=[self.c_rightbound]),
            VertexGroup(name="bndry_bot", tag=14, dim=1, entities=[self.c_botbound]),
            VertexGroup(name="bndry_top_fluid", tag=15, dim=1, entities=top_boundary_entities),

        )
        for i in range(len(self.outer_rise_fault_surface_points)):
            vertex_groups = np.append(vertex_groups, VertexGroup(name="outer_rise_fault_" + str(i), \
                                                                 tag=201 + i, dim=1, entities=[self.c_outer_rise_faults[i]]))

            vertex_groups = np.append(vertex_groups, VertexGroup(name="edge_outer_rise_fault_" + str(i), tag=301 + i, dim=0, \
                                                                 entities=[int(self.outer_rise_fault_buried_edges[i])]))
        for group in vertex_groups:
            group.create_physical_group()

    def generate_mesh(self, cell):
        """Generate the mesh.
        """
        # Set discretization size with geometric progression from distance to the fault.
        # We turn off the default sizing methods.
        gmsh.option.set_number("Mesh.MeshSizeFromPoints", 0)
        gmsh.option.set_number("Mesh.MeshSizeFromCurvature", 0)
        gmsh.option.set_number("Mesh.MeshSizeExtendFromBoundary", 0)

        # First, we setup a field `field_distance` with the distance from the fault.
        fault_distance = gmsh.model.mesh.field.add("Distance")
        gmsh.model.mesh.field.setNumber(fault_distance, "Sampling", 200)

        mesh_refinement_array = np.empty(0)
        for i in range(len(self.c_outer_rise_faults)):
            mesh_refinement_array = np.append(mesh_refinement_array, self.c_outer_rise_faults[i])
        for i in range(len(self.all_curves)):
            mesh_refinement_array = np.append(mesh_refinement_array, self.all_curves[i])
        gmsh.model.mesh.field.setNumbers(fault_distance, "CurvesList", mesh_refinement_array)

        # Second, we setup a field `field_size`, which is the mathematical expression
        # for the cell size as a function of the cell size on the fault, the distance from
        # the fault (as given by `field_size`, and the bias factor.
        # The `GenerateMesh` class includes a special function `get_math_progression` 
        # for creating the string with the mathematical function.

        field_size = gmsh.model.mesh.field.add("MathEval")
        math_exp = GenerateMesh.get_math_progression(fault_distance, min_dx=0.5e+3, bias=1.075)
        gmsh.model.mesh.field.setString(field_size, "F", math_exp)

        ## Finally, we use the field `field_size` for the cell size of the mesh.
        gmsh.model.mesh.field.setAsBackgroundMesh(field_size)

        if cell == "quad":
            gmsh.option.setNumber("Mesh.Algorithm", 8)
            gmsh.model.mesh.generate(2)
            gmsh.model.mesh.recombine()
        else:
            gmsh.model.mesh.generate(2)
        gmsh.model.mesh.optimize("Laplace2D")

if __name__ == "__main__":
    App().main()

# End of file

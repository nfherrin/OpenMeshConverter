# OpenMeshConverter
An open source utility to convert various publicly available mesh formats

This converter takes a version 4 [Gmsh](https://gmsh.info/) file (tested with version 4.1) and converts it to the THOR mesh format.
There are plans to extend this converter to intake other versions of Gmsh and exodus and even perhaps add other output formats in addition to the current \ac{THOR} mesh output.

To compile OpenMeshConverter, navigate to the source folder and then make the converter by typing:
```
  >> make
```
A successful compilation of OpenMeshConverter will conclude with the line:
```
  >> mv ./OpenMeshConverter.exe ../
```

To run OpenMeshConverter, simply invoke the OpenMeshConverter binary and follow it immediately with the Gmsh input file (where `<gmsh_file>` is the name of the Gmsh file):
```
  >> <path_to_OpenMeshConverter>/OpenMeshConverter.exe <gmsh_file>
```
The output file will be titled `<gmsh_file>_out.thrm`. This output will set all boundary conditions to vacuum.
If the user desires to set boundary conditions to reflective or incoming flux boundary conditions, then boundary conditions can be specified on the command line when invoking the Mesh Converter by using the `-bc` indicator.
If the `-bc` indicator is called, then the next six entries will be assumed to be the boundary conditions (integer values) on each of the six primary directions.
The order for the boundary conditions specified in this manner are as follows:
```
 -x +x -y +y -z +z
```
0 is the integer value for vacuum boundary conditions, 1 is the integer value for reflective boundary conditions, and 2 is the integer value for incident flux boundary conditions.
i.e. The following use of OpenMeshConverter will convert the Gmsh file and assign reflective boundary conditions to the $-x$, $-y$, and $+y$ boundary faces, and all other boundary conditions will be set to vacuum.
```
  >> <path_to_Mesh_Converter>/Mesh_Converter.exe <gmsh_file> -bc 1 0 1 1 0 0
```
It should be noted that if reflective boundary conditions are specified, then the reflective boundaries must all reside on flat boundary surfaces.
If the user tries to assign reflective boundary conditions to a direction with a non-flat boundary, then the converter utility will throw an error and terminate.
This check is ignored if all boundary conditions are reflective.

---
## Examples
---

The folder `examples` contains an example of a Gmsh geometry file, `example.geo`, that has used Gmsh to generate a tet mesh file, `example.msh`.
The mesh is them converted to a THOR mesh using the utility, `example.thrm`.
These meshes are for an octant subsection of a sphere in a sphere that takes advantage of spherical symmetry and assigns reflective BCs to each of the negative (flat) directions.
The Gmsh generated mesh can be converted to a THOR mesh by the user if they so wish to compare to the example in the folder.
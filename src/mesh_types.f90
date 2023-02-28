!OpenMeshConverter is licensed under the MIT License.
!-------------------------------------------------------------------------------
!>    This module contains variables and functions used to store and manipulate
!!    data common to both the input mesh and output mesh.
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE mesh_types
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: vertex_type, element_type_3d

    !vertex pointer necessary to make array of pointers
    TYPE :: vert_ptr
      !actual pointer to the vertex type
      TYPE(vertex_type), POINTER :: p
    ENDTYPE

    !the vertex type
    TYPE :: vertex_type
      !location for x,y,z
      REAL(8) :: x=0
      REAL(8) :: y=0
      REAL(8) :: z=0
      !vertex id/index
      INTEGER :: id=0
      CONTAINS
        !compute distance to another point
        PROCEDURE :: distance
    ENDTYPE

    !the base element type
    TYPE, ABSTRACT :: base_element_type
      !element region
      INTEGER :: reg=0
      !element id/index
      INTEGER :: id=0
    ENDTYPE

    !the specific 3d element type (a tetrahedron)
    TYPE, EXTENDS(base_element_type) :: element_type_3d
      !4 corners of the tet
      TYPE(vert_ptr) :: corner(4)
      !Adjacent tet id for faces 1 to 4
      INTEGER :: adj_id(4)=0
      !Adjacent tet face for faces 1 to 4
      INTEGER :: adj_face(4)=0
    ENDTYPE
  CONTAINS

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !computes the distance to another vertex
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    REAL(8) FUNCTION distance(a,b)
      CLASS(vertex_type), INTENT(IN) :: a
      TYPE(vertex_type), INTENT(IN) :: b

      distance=SQRT((a%x-b%x)**2+(a%y-b%y)**2+(a%z-b%z)**2)
    ENDFUNCTION distance

  END MODULE mesh_types

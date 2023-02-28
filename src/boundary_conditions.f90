!OpenMeshConverter is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to create boundary
!! conditions and adjacency info for a given set of elements
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE boundary_conditions
    USE globals
    IMPLICIT NONE
    PRIVATE
    PUBLIC :: adjacency_calc

    INTEGER, ALLOCATABLE :: tbound_cond(:,:)
    INTEGER, ALLOCATABLE :: bc_side(:)
CONTAINS

  SUBROUTINE adjacency_calc()
    INTEGER :: i,og_face(3),adj_idx,j,ii,jj
    LOGICAL :: adj_found=.FALSE.

    DO i=1,num_tets
      CALL orderverts(tet(i))
    ENDDO

    ALLOCATE(tbound_cond(num_tets*4,3))
    tbound_cond=0
    !loop over all tets
    adj_idx=0
    num_bcf=0
    prog=0
    WRITE(*,'(A)',ADVANCE='NO')'Progress:'
    !loop over all tets
    DO i=1,num_tets
      IF(MOD(i,CEILING(num_tets*1.0/(max_prog-1.0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
      !loop over all faces
      DO j=1,4
        !we're not exiting the face yet (we do once we find the adjacency)
        adj_found=.FALSE.
        !make sure that the adjacency hasn't already been found
        IF(tet(i)%adj_id(j) .EQ. 0)THEN
          !loop over all other tets
          DO ii=1,num_tets
            !make sure they're not the same element
            IF(ii .NE. i)THEN
              !loop over all faces for the other tet
              DO jj=1,4
                !if the sides match, then assign the adjacencies
                IF(check_face(tet(i),j,tet(ii),jj))THEN
                  tet(i)%adj_id(j)=ii
                  tet(ii)%adj_id(jj)=i
                  tet(i)%adj_face(j)=jj
                  tet(ii)%adj_face(jj)=j
                  adj_found=.TRUE.
                ENDIF
                IF(adj_found)EXIT
              ENDDO
            ENDIF
            IF(adj_found)EXIT
          ENDDO
          !if we never found an adjacency, then the face is a bc face
          IF(.NOT. adj_found)THEN
            num_bcf=num_bcf+1
            tbound_cond(num_bcf,1)=i
            tbound_cond(num_bcf,2)=j
            !assign bc to 0 for now (this will change later)
            tbound_cond(num_bcf,3)=0
          ENDIF
        ENDIF
      ENDDO
    ENDDO
    ALLOCATE(bc_data(num_bcf,3),bc_side(num_bcf))
    bc_data=0
    bc_side=0
    DO i=1,num_bcf
      bc_data(i,:)=tbound_cond(i,:)
    ENDDO

    bc_side=0
    !determine which sides are flat
    IF(MINVAL(side_bc) .EQ. MAXVAL(side_bc))THEN
      bc_data(:,3)=side_bc(1)
    ELSE
      CALL det_side_flatness()
    ENDIF
    DO i=prog,max_prog
      WRITE(*,'(A)',ADVANCE='NO')'*'
    ENDDO
    WRITE(*,*)
    DEALLOCATE(tbound_cond,bc_side)
  ENDSUBROUTINE adjacency_calc

  LOGICAL FUNCTION check_face(tet1,face1,tet2,face2)
    TYPE(element_type_3d),INTENT(IN) :: tet1,tet2
    INTEGER,INTENT(IN) :: face1,face2
    INTEGER :: ids1(3)=0,ids2(3)=0,j=0,i=0

    check_face=.FALSE.

    !assign face ids for face1
    i=0
    DO j=1,4
      IF(j .NE. face1)THEN
        i=i+1
        ids1(i)=tet1%corner(j)%p%id
      ENDIF
    ENDDO
    !assign face ids for face2
    i=0
    DO j=1,4
      IF(j .NE. face2)THEN
        i=i+1
        ids2(i)=tet2%corner(j)%p%id
      ENDIF
    ENDDO

    IF(ids1(1) .EQ. ids2(1) .AND. ids1(2) .EQ. ids2(2) .AND. ids1(3) .EQ. ids2(3))check_face=.TRUE.
  ENDFUNCTION check_face

  SUBROUTINE det_side_flatness()
    INTEGER :: i,j,el_id,face_idx
    REAL(8) :: face_point(3,3),ext_point(3),norm_vec(3),lambda,offset

    !sides are assumed to be flat, and if they are not this is set to false.
    side_flat=.TRUE.
    DO i=1,num_bcf
      el_id=bc_data(i,1)
      face_idx=0
      !assign extruded and face boundary points
      DO j=1,4
        IF(bc_data(i,2) .EQ. j)THEN
          !if it equals the face index the it's the extruded point
          ext_point(:)=(/tet(el_id)%corner(j)%p%x,tet(el_id)%corner(j)%p%y,&
              tet(el_id)%corner(j)%p%z/)
        ELSE
          !if doesn't equal the face index the it's one of the face points
          face_idx=face_idx+1
          face_point(face_idx,:)=(/tet(el_id)%corner(j)%p%x, &
              tet(el_id)%corner(j)%p%y,tet(el_id)%corner(j)%p%z/)
        ENDIF
      ENDDO
      !get the outward going unit normal vector for the tet for this face
      norm_vec=cross(face_point(2,:)-face_point(1,:), face_point(3,:)-face_point(1,:))
      offset=face_point(1,1)*norm_vec(1)+face_point(1,2)*norm_vec(2)+face_point(1,3)*norm_vec(3)
      lambda=(offset-norm_vec(1)*ext_point(1)-norm_vec(2)*ext_point(2)-norm_vec(3)*ext_point(3)) &
          /(norm_vec(1)**2+norm_vec(2)**2+norm_vec(3)**2)
      norm_vec=norm_vec*lambda
      norm_vec=norm_vec/(SQRT(norm_vec(1)**2+norm_vec(2)**2+norm_vec(3)**2))

      !figure out which side of the problem this bc faces
      IF(ABS(MAXVAL(norm_vec)) .GT. ABS(MINVAL(norm_vec)))THEN
        !face on a positive side, assign the side for that face
        IF(MAXLOC(norm_vec,1) .EQ. 1)THEN
          bc_side(i)=2
        ELSEIF(MAXLOC(norm_vec,1) .EQ. 2)THEN
          bc_side(i)=4
        ELSE
          bc_side(i)=6
        ENDIF
        IF(ABS(MAXVAL(norm_vec)-1.0) .LE. 1.0E-14)THEN
          !face is flat, do nothing
        ElSE
          !face is not flat, it only takes 1 for the side to not be flat
          side_flat(bc_side(i))=.FALSE.
        ENDIF
      ELSE
        !face on a negative side, assign the side for that face
        IF(MINLOC(norm_vec,1) .EQ. 1)THEN
          bc_side(i)=1
        ELSEIF(MINLOC(norm_vec,1) .EQ. 2)THEN
          bc_side(i)=3
        ELSE
          bc_side(i)=5
        ENDIF
        IF(ABS(MINVAL(norm_vec)+1.0) .LE. 1.0E-14)THEN
          !face is flat, do nothing
        ElSE
          !face is not flat, it only takes 1 for the side to not be flat
          side_flat(bc_side(i))=.FALSE.
        ENDIF
      ENDIF
    ENDDO
    DO i=1,6
      IF(.NOT. side_flat(i) .AND. side_bc(i) .EQ. 1)THEN
        WRITE(*,*)'ERROR: side ',i,' is reflective but not flat'
        STOP
      ENDIF
    ENDDO
    DO i=1,num_bcf
      bc_data(i,3)=side_bc(bc_side(i))
    ENDDO
  ENDSUBROUTINE det_side_flatness

  FUNCTION cross(a, b)
    REAL(8) :: cross(3)
    REAL(8), INTENT(IN) :: a(3), b(3)

    cross(1) = a(2) * b(3) - a(3) * b(2)
    cross(2) = a(3) * b(1) - a(1) * b(3)
    cross(3) = a(1) * b(2) - a(2) * b(1)
  END FUNCTION cross

  SUBROUTINE orderverts(this_tet)
    TYPE(element_type_3d), INTENT(INOUT) :: this_tet
    TYPE(vertex_type), POINTER :: temp_vert
    INTEGER :: i,changes

    !bubble sort algorithm, pretty cheap for only 4 points
    DO
      changes=0
      DO i=1,3
        IF(this_tet%corner(i)%p%id .GE. this_tet%corner(i+1)%p%id)THEN
          temp_vert => vertex(this_tet%corner(i)%p%id)
          this_tet%corner(i)%p => vertex(this_tet%corner(i+1)%p%id)
          this_tet%corner(i+1)%p => vertex(temp_vert%id)
          changes=changes+1
        ENDIF
      ENDDO
      IF(changes .EQ. 0)EXIT
    ENDDO
  ENDSUBROUTINE orderverts
END MODULE boundary_conditions

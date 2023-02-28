!OpenMeshConverter is licensed under the MIT License.
!-------------------------------------------------------------------------------
!> This module contains the functionality necessary to output a file in the
!! Thor_mesh format (.thrm)
!> @author Nicholas F. Herring
!-------------------------------------------------------------------------------
MODULE output_thrm
  USE globals
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: output_thrm_file, calcvols
CONTAINS

  SUBROUTINE output_thrm_file()
    INTEGER :: i,j
    WRITE(*,'(A)',ADVANCE='NO')'Progress:'

    OPEN(UNIT=30,FILE=TRIM(ADJUSTL(mesh_outfile)),ACTION='WRITE',STATUS='REPLACE')

    !print out base data
    WRITE(30,'(I0)')num_verts
    WRITE(30,'(I0)')num_tets
    WRITE(30,'(I0)')1
    WRITE(30,'(I0)')1
    !print out vertices
    DO i=1,num_verts
      WRITE(30,'(I0,3ES24.16)')i,vertex(i)%x,vertex(i)%y,vertex(i)%z
    ENDDO
    !print out tet regions and source. Assume each region has its own source. User can change this later
    DO i=1,num_tets
      WRITE(30,'(I0,A,I0,A,I0)')i,' ',tet(i)%reg,' ',tet(i)%reg
    ENDDO
    !print out tet composition
    DO i=1,num_tets
      WRITE(30,'(I0,A,I0,A,I0,A,I0,A,I0)')i,' ',tet(i)%corner(1)%p%id,' ',tet(i)%corner(2)%p%id,' ', &
          tet(i)%corner(3)%p%id,' ',tet(i)%corner(4)%p%id
    ENDDO
    !print out boundary conditions
    WRITE(30,'(I0)')num_bcf
    DO i=1,num_bcf
      WRITE(30,'(I0,A,I0,A,I0)')bc_data(i,1),' ',bc_data(i,2)-1,' ',bc_data(i,3)
    ENDDO
    !print out adjacency list
    WRITE(30,'(I0)')num_tets*4
    prog=0
    DO i=1,num_tets
      IF(MOD(i,CEILING(num_tets/(max_prog-1.0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
      DO j=1,4
        WRITE(30,'(I0,A,I0,A,I0,A,I0)')tet(i)%id,' ',j-1,' ',tet(i)%adj_id(j),' ' &
            ,MAX(tet(i)%adj_face(j)-1,0)
      ENDDO
    ENDDO

    CLOSE(30)
    DO i=prog,max_prog
      WRITE(*,'(A)',ADVANCE='NO')'*'
    ENDDO
    WRITE(*,*)
  ENDSUBROUTINE output_thrm_file

  SUBROUTINE calcvols()
    REAL(8) :: totalvol1
    TYPE(vertex_type) :: a,b,c,d
    REAL(8), ALLOCATABLE :: tetvol(:),regvol(:)
    INTEGER, ALLOCATABLE :: tets_in_reg(:)
    INTEGER :: i,minreg,maxreg

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'

    minreg=MINVAL(tet(:)%reg)
    maxreg=MAXVAL(tet(:)%reg)
    ALLOCATE(tetvol(num_tets),regvol(minreg:maxreg),tets_in_reg(minreg:maxreg))
    tets_in_reg=0
    tetvol=0
    totalvol1=0
    regvol=0
    prog=0
    !compute tet volumes and add to both total volumes and region volumes
    DO i=1,num_tets
      a=tet(i)%corner(1)%p
      b=tet(i)%corner(2)%p
      c=tet(i)%corner(3)%p
      d=tet(i)%corner(4)%p
      tetvol(i)=ABS((-c%y*d%x+b%y*(-c%x+d%x)+b%x*(c%y-d%y)+c%x*d%y)*(a%z-d%z)+(a%x-d%x) &
          *(-c%z*d%y+b%z*(-c%y+d%y)+b%y*(c%z-d%z)+c%y*d%z)+(a%y-d%y)*(b%z*(c%x-d%x) &
          +c%z*d%x-c%x*d%z+b%x*(-c%z+d%z)))/6.0D0
      regvol(tet(i)%reg)=regvol(tet(i)%reg)+tetvol(i)
      tets_in_reg(tet(i)%reg)=tets_in_reg(tet(i)%reg)+1
      totalvol1=totalvol1+tetvol(i)
      IF(MOD(i,CEILING(num_tets*1.0/(max_prog-1.0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
    ENDDO
    DO i=prog,max_prog
      WRITE(*,'(A)',ADVANCE='NO')'*'
    ENDDO
    WRITE(*,*)

    DO i=minreg,maxreg
      WRITE(*,'(A,I0,A,I0)')'Region ',i,' tets: ',tets_in_reg(i)
      WRITE(*,'(A,I0,A,ES24.16)')'Region ',i,' volume: ',regvol(i)
    ENDDO
    WRITE(*,'(A,I0)')'Total number of tets: ',SUM(tets_in_reg)
    WRITE(*,'(A,ES24.16)')'Total system volume: ',totalvol1
  ENDSUBROUTINE calcvols
END MODULE output_thrm

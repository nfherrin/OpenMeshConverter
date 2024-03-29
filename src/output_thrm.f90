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
    INTEGER :: i
    WRITE(*,'(A)',ADVANCE='NO')'Progress:'

    OPEN(UNIT=30,FILE=TRIM(ADJUSTL(mesh_outfile)),ACTION='WRITE',STATUS='REPLACE')

    !print out base data
    WRITE(30,'(I0)')num_verts
    WRITE(30,'(I0)')num_tets
    WRITE(30,'(I0)')1
    WRITE(30,'(I0)')1
    !print out vertices
    DO i=1,num_verts
      WRITE(30,'(I0,3ES24.16)')i,vertex(i,:)
    ENDDO
    !print out tet regions and source. Assume each region has its own source. User can change this later
    DO i=1,num_tets
      WRITE(30,'(I0,A,I0,A,I0)')i,' ',el_tag(i),' ',el_tag(i)
    ENDDO
    !print out tet composition
    DO i=1,num_tets
      WRITE(30,'(I0,A,I0,A,I0,A,I0,A,I0)')i,' ',element(i,1),' ',element(i,2),' ',element(i,3),' ' &
          ,element(i,4)
    ENDDO
    !print out boundary conditions
    WRITE(30,'(I0)')num_bcf
    DO i=1,num_bcf
      WRITE(30,'(I0,A,I0,A,I0)')bc_data(i,1),' ',bc_data(i,2),' ',bc_data(i,3)
    ENDDO
    !print out adjacency list
    WRITE(30,'(I0)')num_tets*4
    prog=0
    DO i=1,num_tets*4
      IF(MOD(i,CEILING(num_tets*4.0/(max_prog-1.0))) .EQ. 0)THEN
        WRITE(*,'(A)',ADVANCE='NO')'*'
        prog=prog+1
      ENDIF
      WRITE(30,'(I0,A,I0,A,I0,A,I0)')adj_list(i,1),' ',adj_list(i,2),' ',adj_list(i,3),' ',adj_list(i,4)
    ENDDO

    CLOSE(30)
    DO i=prog,max_prog
      WRITE(*,'(A)',ADVANCE='NO')'*'
    ENDDO
    WRITE(*,*)
  ENDSUBROUTINE output_thrm_file

  SUBROUTINE calcvols()
    REAL(8) :: totalvol1
    REAL(8) :: a(3),b(3),c(3),d(3)
    REAL(8), ALLOCATABLE :: tetvol(:),regvol(:)
    INTEGER, ALLOCATABLE :: tets_in_reg(:)
    INTEGER :: i,minreg,maxreg

    WRITE(*,'(A)',ADVANCE='NO')'Progress:'

    minreg=MINVAL(el_tag(:))
    maxreg=MAXVAL(el_tag(:))
    ALLOCATE(tetvol(num_tets),regvol(minreg:maxreg),tets_in_reg(minreg:maxreg))
    tets_in_reg=0
    tetvol=0
    totalvol1=0
    regvol=0
    prog=0
    !compute tet volumes and add to both total volumes and region volumes
    DO i=1,num_tets
      a(:)=vertex(element(i,1),:)
      b(:)=vertex(element(i,2),:)
      c(:)=vertex(element(i,3),:)
      d(:)=vertex(element(i,4),:)
      tetvol(i)=ABS((-c(2)*d(1)+b(2)*(-c(1)+d(1))+b(1)*(c(2)-d(2))+c(1)*d(2))*(a(3)-d(3))+(a(1)-d(1)) &
          *(-c(3)*d(2)+b(3)*(-c(2)+d(2))+b(2)*(c(3)-d(3))+c(2)*d(3))+(a(2)-d(2))*(b(3)*(c(1)-d(1)) &
          +c(3)*d(1)-c(1)*d(3)+b(1)*(-c(3)+d(3))))/6
      regvol(el_tag(i))=regvol(el_tag(i))+tetvol(i)
      tets_in_reg(el_tag(i))=tets_in_reg(el_tag(i))+1
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
      IF(tets_in_reg(i) .GT. 0)THEN
        WRITE(*,'(A,I0,A,I0)')'Region ',i,' tets: ',tets_in_reg(i)
        WRITE(*,'(A,I0,A,ES24.16)')'Region ',i,' volume: ',regvol(i)
      ENDIF
    ENDDO
    WRITE(*,'(A,I0)')'Total number of tets: ',SUM(tets_in_reg)
    WRITE(*,'(A,ES24.16)')'Total system volume: ',totalvol1
  ENDSUBROUTINE calcvols
END MODULE output_thrm

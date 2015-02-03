! ====================
!  README Gmsh2SU3D ver 0.02
! ====================
!
! Original version by Ceanwang@gmail.com, Jan 2012
!
! Adapted from gmsh2dolfyn.f90 developed by dolfyn team.
! For dolfyn, please visit http://www.dolfyn.net/index_en.html
!
! Support Gmsh 2.5.0
!
! Purpose
! -------
! This Fortran95 program translates a mesh file from Gmsh (.msh) format
! to Su2 format.
!
! Input and Output
! ----------------
! Input : A Gmsh .msh file (version 2.0, ascii format).
! Output: SU2 files.
!
! Running the Program
!--------------------
! First compile it using a Fortran95 compiler (eg g95 or gfortran).
! Run it from the command line.
! The program prompts for the name of the input file.
!
! Bug reports
! -----------
! Please report bugs to ceanwang@gmail.com
!
! Important note about the Gmsh msh format and Physical Groups.
! -------------------------------------------------------------
! In order to define boundary conditions, the Gmsh geometry-builder allows a
! group of faces to be assigned a common 'physical group' label.  The mesh
! inherits this label, and the label is used in the .su2 file.
!
! When saving the mesh, the default is to save only mesh elements with a
! physical group label.  This means that some mesh elements will be missing,
! unless every mesh element belongs to a physical group.
!
! For example in the adapted gmsh tutorial t2.geo enter:
!
! Physical Volume ("Fluid") = {119,120};
! Physical Surface("Inlet") = {111};
! Physical Surface("Outlet") = {132};
!
! Mesh 3D and save it as t2.msh
!
!========================================================================
!========================================================================
!========================================================================
SUBROUTINE UPPERCASE(STR)

      IMPLICIT NONE

      CHARACTER(LEN=*), INTENT(IN OUT) :: STR
      INTEGER :: I, DEL

      DEL = IACHAR('a') - IACHAR('A')

      DO I = 1, LEN_TRIM(STR)
         IF (LGE(STR(I:I),'a') .AND. LLE(STR(I:I),'z')) THEN
            STR(I:I) = ACHAR(IACHAR(STR(I:I)) - DEL)
         END IF
      END DO

      RETURN

END SUBROUTINE UPPERCASE

integer function lens(string)

   character(len=*) string

   do i=len(string),0,-1
     if( string(i:i) .ne. ' ') goto 10
   end do
   i = 0
10 continue

   lens = i

end function lens
!======================================================================
subroutine openfile(iunit,casename,extension,reqform,status,idebug)

   character(len=*)  casename
   character(len=*)  extension
   character(len=*)  reqform
   character(len=*)  status
   character(len=48) filename
   character(len=11) form

   logical exists

   filename = casename(1:lens(casename))//extension(1:lens(extension))
   length   = lens(filename)

   if( idebug > 2 )write(*,*) 'Opening ',filename(1:length)

   if( status(1:3) == 'OLD' )then
     inquire(file=filename(1:length),exist=exists,form=form)
     if( .not. exists )then
       write(*,*) '*** Error: File ',filename(1:length),' does not exist'
       stop
     endif
   endif

   open(iunit,file=filename(1:length),form=reqform,status=status)

   if( idebug >= 2 ) write(*,*) 'File ',filename(1:length),' opened'

end subroutine openfile
!==============================================================================
program gmsh2SU2

   implicit none

   integer, parameter :: IOinp = 13, IOcel = 14   ! I/O file numbers
   integer, parameter :: IOdbg = 63, IOcfg = 12
   integer, parameter :: IOgmsh= 24                                      !Gmsh mesh file
   
   integer :: Ninlet = 0                                     
   integer :: Noutlet = 0                                     
   integer :: Nsurface = 0                                     
   integer :: isur = 0                                     

   integer :: debug = 0

   integer, parameter :: version = 0530

   character(len=128) :: line

   integer, parameter :: MaxNames = 100
   integer, parameter :: MaxNperBnd = 500
   integer, parameter :: MaxNodes = 90000
   
   character(len=64), dimension(MaxNames) :: Names
   character(len=64), dimension(MaxNames) :: Regions
   integer, dimension(MaxNames)           :: ICTID    = -1
   integer, dimension(MaxNames)           :: Partition=  1

   logical, dimension(MaxNames)           :: Fluid    = .false.
   logical, dimension(MaxNames)           :: Boundary = .false.

   character(len=64)  :: casename = 'su2'
   character(len=72)  :: c_input1, c_input2, c_input3

   integer i, j, k, ie, icel, ibnd, iloop,ii
   integer tbnd(maxnames)
   integer nbnd
   !
   ! nodes/vertices
   !
   integer n_nodes,inode
   real  node(MaxNodes,3)
   integer  mytags(MaxNames)
   integer  tv0(MaxNperBnd,MaxNames)
   integer  tv1(MaxNperBnd,MaxNames)
   integer  tv2(MaxNperBnd,MaxNames)
   integer  tv3(MaxNperBnd,MaxNames)
   integer  tv4(MaxNperBnd,MaxNames)

   !
   ! there are 19 gmsh element types: \
   !
   !    1 : 2-node line
   !    2 : 3-node triangle (face)
   !    3 : 4-node quadrangle (face)
   !    4 : 4-node tetrahedron
   !    5 : 8-node hexahedron (eg cube)
   !    6 : 6-node triangular-prism
   !    7 : 5-node pyramid
   !
   !  8-14: 'second-order' elements.  Ref Gmsh manual.
   !   15 : 1-node point
   ! 16-19: more second-order FEM elements
   !
   ! the nodes/vertices for each element are read into the
   ! v array.
   !
   ! each element can have several tags.
   ! the first tag gives the physical_group number.
   ! all elements on the same boundary have the same physical_group number.
   !
   integer, parameter :: element_type(19) = &
                      (/ 2,3,4,4,8,6,5,3,6,9,10,27,18,14,1,8,20,15,13 /)
			
   integer :: n_elements, ielement, ielement_type, n_tags, n_names, lens
   integer :: tags(64), v(27)
   integer :: bmarknew,bmarkold
   integer :: i3, i4q, i4, i5, i6, i8
   integer :: i2
   integer :: ivs = 0
   


   if( size(v) /= maxval(element_type) )then
     stop'bug: error in dimensions of array v'
   endif

   !
   ! read the gmsh filename, then open the .msh file
   !
   write(*,*) 'Gmsh2SU2: Converts a Gmsh mesh file to SU2 format.'
   write(*,*) '(Input must be in Gmsh version 2.0 ascii format.'
   write(*,*) ' Output is in SU2 format.)'
   write(*,*) ' '
   write(*,*) 'Input Gmsh filename, excluding the .msh suffix'
   read(*,'(A)') casename

   write(*,*) 'Opening the Gmsh file'
   call openfile(IOgmsh,casename,'.msh','FORMATTED','OLD',debug)

   !
   ! read the Gmsh file header
   !
   write(*,*)'Reading MeshFormat'
   read(IOgmsh,*) c_input1
   call check_input_character(c_input1,'$MeshFormat')

   read(IOgmsh,*) c_input1,c_input2,c_input3
   if( c_input1 == '2.2' )then
     ivs = 22
   else if( c_input1 == '2.1' )then
     ivs = 21
   else if( c_input1 == '2' )then
     ivs = 20
   else
     write(*,*) '*** WARNING: unknown Gmsh version'
     write(*,*) '*** Unexpected results might happen'
     ivs = 21
   endif

   if( ivs == 20 )then
     call check_input_character(c_input1,'2')
   else if( ivs == 21 )then
     call check_input_character(c_input1,'2.1')
   else if( ivs == 22 )then
     call check_input_character(c_input1,'2.2')
   else
     write(*,*) '*** Version found ',c_input1
   endif

   call check_input_character(c_input2,'0')
   call check_input_character(c_input3,'8')

   write(*,*) 'MeshFormat: ', c_input1(1:lens(c_input1)),' ', &
                              c_input2(1:lens(c_input2)),' ', &
                              c_input3(1:lens(c_input3))

   read(IOgmsh,*) c_input1
   call check_input_character(c_input1,'$EndMeshFormat')

   !
   ! read the Gmsh PhysicalNames
   !
   write(*,*)'Reading PhysicalNames'
   read(IOgmsh,*) c_input1
   call check_input_character(c_input1,'$PhysicalNames')

   read(IOgmsh,*) n_names
   if( n_names <= 0 )then
     write(*,*) 'error: number of names must be a positive number'
     stop
   endif

   if( ivs == 20 )then
     do i=1,n_names
      !read(IOgmsh,*) k,j,c_input1
       read(IOgmsh,*)   j,c_input1
       write(*,*) 'Name ',j,'-> ', c_input1
       Names(j) = c_input1
     end do
   else if( ivs == 21 )then
     do i=1,n_names
       read(IOgmsh,*) k,j,c_input1
       write(*,*) 'Name ',j,'-> ', c_input1
       Names(j) = c_input1
     end do
   else
     do i=1,n_names
       read(IOgmsh,*) k,j,c_input1
       write(*,*) 'Name ',j,'-> ', c_input1
       Names(j) = c_input1
     end do
   endif


   do i=1,n_names
       CALL UPPERCASE(Names(i))
   enddo

   read(IOgmsh,*) c_input1
   call check_input_character(c_input1,'$EndPhysicalNames')
   !
   ! read the nodes from the .msh file and write them
   ! to the .vrt file.
   !
   write(*,*)'Reading Nodes'
   read(IOgmsh,*) c_input1
   call check_input_character(c_input1,'$Nodes')
   !nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn
   read(IOgmsh,*) n_nodes
   if( n_nodes <= 0 )then
     write(*,*) 'error: number of nodes must be a positive number'
     stop
   endif

   if( n_nodes > MaxNodes )then
       write(*,*) 'error: The Gmsh file contains ',n_nodes,' nodes.'
       write(*,*) 'Gmsh2Su2 is hard-wired for a maximum of ',MaxNodes,&
        	  'nodes.  The dimension of this array needs to be increased.'
     stop
   endif

   !
   ! open the su2 .vrt.su2 file
   !
   !write(*,*) 'Creating the su2 .vrt.su2 file'
   !call openfile(IOvrt,casename,'.vrt.su2','FORMATTED','UNKNOWN',debug)

   nodes: do iloop=1,n_nodes
     read(IOgmsh,*) inode,(node(iloop,i), i=1,3)
     !write(IOvrt,'(3g16.9,6x,i9)') (node(iloop,i),i=1,3),inode-1
   enddo nodes

   write(*,*) 'Nodes written ',n_nodes
   !
   ! close the su2 .vrt.su2 file
   !
   !close(IOvrt)

   read(IOgmsh,*) c_input1
   call check_input_character(c_input1,'$EndNodes')
   !eeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeeee
   !
   ! read the elements from the .msh file and write them
   ! to the .cel and .bnd files.
   !
   read(IOgmsh,*) c_input1
   call check_input_character(c_input1,'$Elements')

   read(IOgmsh,*) n_elements
   if( n_elements <= 0 )then
     write(*,*) 'error: number of elements must be a positive number'
     stop
   endif

   write(*,*) 'Total Gmsh elements to be read in:',n_elements
   !
   ! open the su2 .cel files
   !
   write(*,*) 'Creating the .su2 files'
   call openfile(IOcel,casename,'.su2','FORMATTED','UNKNOWN',debug)
   write(IOcel,101) 3
101 format('NDIME= ',i1)
 

   !
   ! note in Gmsh fluid cells and boudaries can be mixed
   ! we just keep track on them both
   ! remind default region is not assigned
   !
   icel = 0
   ibnd = 0
   nbnd=0
   tbnd(nbnd)=0

   i2   = 0
   i3   = 0
   i4q  = 0
   i4   = 0
   i5   = 0
   i6   = 0
   i8   = 0
   
   bmarkOld=0
   do ie=1,n_elements

     read(IOgmsh,*) ielement, ielement_type, n_tags
     if( ivs <= 21 )then
       if( n_tags /= 3 ) write(*,*) 'tag error n_tags /= 3:',ielement,n_tags
     else
       if( n_tags /= 2 ) write(*,*) 'tag error n_tags /= 2:',ielement,n_tags     
     endif
     call check_element_type(ielement_type,element_type)
     call check_n_tags(n_tags,tags)
     backspace(IOgmsh)

     !
     ! we need to circumvent backspace but
     ! advance='no' requires fixed format
     ! just keep it for now.
     !

     !
     ! now we know what to to expect to find on the line
     !
     read(IOgmsh,*) ielement, ielement_type, &
     		    n_tags, (tags(i),i=1,n_tags),&
     		   (v(i),i=1,element_type(ielement_type))
           
     do ii=1,element_type(ielement_type)      
         v(ii)=v(ii)-1
     end do
           
     bmarkNew=tags(1)
     if( 4 <= ielement_type .and. ielement_type <= 7 )then
       if (icel==0) then
           write(IOcel,121) n_elements-ibnd
121 format('NELEM= ',i10)
       endif
       icel = icel + 1

       if( .not. Fluid(tags(1)) ) Fluid(tags(1)) = .true.
       if( Boundary(tags(1)) )then
         write(*,*) 'Inconsistent data: Physical names ids overlap 1'
       endif

     1 format(i8,8(1x,i8),2(1x,i4))
       select case(ielement_type)

         case(4) ! 4-node tet

           write(IOcel,*) 10,v(1),v(2),v(3),v(4),icel-1
           i4 = i4 + 1

         case(5) ! 8-node hex

           write(IOcel,*) 12,v(1),v(2),v(3),v(4),  v(5),v(6),v(7),v(8),icel-1
     	     	!tags(1),tags(3)
           i8 = i8 + 1

         case(6) ! 6-node prism or wedge

           write(IOcel,*) 13,v(1),v(2),v(3), v(4),v(5),v(6),icel-1
           i6 = i6 + 1

         case(7) ! 5-node pyramid

           write(IOcel,*) 14,v(1),v(2),v(3),v(4),  v(5),icel-1
           i5 = i5 + 1

         case default
	
	   write(*,*)'internal error 1'

       end select

     elseif( ielement_type == 2 .or. ielement_type == 3 )then

       if (bmarkNew/=bmarkOld) then
          bmarkOld=bmarkNew
          nbnd=nbnd+1
          tbnd(nbnd)=0
       endif
       
       ibnd = ibnd + 1
       tbnd(nbnd) = tbnd(nbnd) + 1
       
       
       if( .not. Boundary(tags(1)) ) Boundary(tags(1)) = .true.
       if( Fluid(tags(1)) )then
         write(*,*) 'Inconsistent data: Physical names ids overlap 2'
       endif

       select case(ielement_type)

         case(2) ! 3-node tri

           !write(IObnd,*) 5, v(1),v(2),v(3)
           mytags(nbnd)=tags(1)
           tv0(tbnd(nbnd),nbnd)=5
           tv1(tbnd(nbnd),nbnd)=v(1)
           tv2(tbnd(nbnd),nbnd)=v(2)
           tv3(tbnd(nbnd),nbnd)=v(3)
           i3 = i3 + 1

         case(3) ! 4-node quad

           !write(IObnd,1) 8, v(1),v(2),v(3),v(4)
           mytags(nbnd)=tags(1)
           tv0(tbnd(nbnd),nbnd)=8
           tv1(tbnd(nbnd),nbnd)=v(1)
           tv2(tbnd(nbnd),nbnd)=v(2)
           tv3(tbnd(nbnd),nbnd)=v(3)
           tv4(tbnd(nbnd),nbnd)=v(4)
           i4q = i4q + 1

         case default
	
	   write(*,*)'internal error 2'

       end select

     elseif( ielement_type == 1 )then
         ! 2- nodes line

           !write(IObnd,*) 3, v(1),v(2)
           mytags(nbnd)=tags(1)
           tv0(tbnd(nbnd),nbnd)=3           
           tv1(tbnd(nbnd),nbnd)=v(1)
           tv2(tbnd(nbnd),nbnd)=v(2)
           i2 = i2 + 1

     else
     
       write(*,*)'internal error 3'

     endif

   end do

   read(IOgmsh,*) c_input1
   call check_input_character(c_input1,'$EndElements')
   !------------------------------------------------------------
   ! write out points
   !------------------------------------------------------------
   write(IOcel,91) n_nodes
91 format('NPOIN= ',i10)

   do iloop=1,n_nodes
     write(IOcel,'(3g16.9,6x,i9)') (node(iloop,i),i=1,3),iloop-1
   enddo 
   
   write(IOcel,111) n_names-1
111 format('NMARK= ',i10)
   do j=1,n_names-1
       write(IOcel,141) mytags(j)
       write(IOcel,151) tbnd(j)
       do i=1,tbnd(j)
            if  (tv0(i,j)==5) then
                write(IOcel,*) tv0(i,j), tv1(i,j),tv2(i,j),tv3(i,j)
            endif
201 format(1x, i10,3f15.6)                
            if  (tv0(i,j)==8) then
                write(IOcel,*) tv0(i,j), tv1(i,j),tv2(i,j),tv3(i,j),tv4(i,j)
            endif
211 format(1x, i10,4f15.6)

       enddo
   end do
141 format('MARKER_TAG= ',i3)
151 format('MARKER_ELEMS= ', i10)
   close(IOcel)


   if( i3  > 0 ) write(*,*) 'Triangle boundaries: ',i3
   if( i4q > 0 ) write(*,*) 'Quad boundaries:     ',i4q
   if( i4  > 0 ) write(*,*) 'Tetrahedral cells:   ',i4
   if( i5  > 0 ) write(*,*) 'Pyramid cells:       ',i5
   if( i6  > 0 ) write(*,*) 'Prism cells:         ',i6
   if( i8  > 0 ) write(*,*) 'Hexahedral cells:    ',i8

   !------------------------------------------------------------
   ! finally write out the boundary names
   !------------------------------------------------------------
   do i=1,n_Names
       if (Names(i)=='INLET') then 
          Ninlet=Ninlet+1
       endif   
       if (Names(i)=='OUTLET') then 
           Noutlet=Noutlet+1
       endif
   enddo
   Nsurface=N_names-Ninlet-Noutlet-1
   
   write(*,*) 'Writing the .inp file'
   call openfile(IOinp,casename,'_inp.txt','FORMATTED','UNKNOWN',debug)

   write(IOinp,'('''')') 
   write(IOinp,'(''% -------------------- BOUNDARY CONDITION DEFINITION --------------------------%'')') 
   write(IOinp,'(''%'')') 
   write(IOinp,'(''% Euler wall boundary marker(s) (NONE = no marker)'')') 
   !write(IOinp,'(''MARKER_EULER='')')  !( 3,4,5,6,7,8,9,10 )
   if (N_names>0) then
    write(IOinp,'( "MARKER_EULER=(" )',ADVANCE = "NO") 
    do i=1,N_names
       if (Names(i) .ne. 'INLET' .and. Names(i) .ne. 'OUTLET' .and. Names(i) .ne. 'FLUID') then
          write(IOinp,'( i4 )',ADVANCE = "NO") i
          isur=isur+1
          if (isur<Nsurface) then
              write(IOinp,'( "," )',ADVANCE = "NO") 
          endif
       endif
    end do  
    write(IOinp,'( ")" )') 
   else
           write(IOinp,'(''MARKER_EULER = NONE'')') 
   endif
   write(IOinp,'(''%'')') 
   
   write(IOinp,'(''% Inlet boundary marker(s) (NONE = no marker) '')') 
   write(IOinp,'(''% Format: ( inlet marker, total temperature, total pressure, flow_direction_x, '')') 
   write(IOinp,'(''%           flow_direction_y, flow_direction_z, ... ) where flow_direction is'')') 
   write(IOinp,'(''%           a unit vector.'')') 
   if (Ninlet>0) then
    do i=1,MaxNames
       if (Names(i)=='INLET') then
          write(IOinp,'( "MARKER_INLET=(",i4,",288.6, 102010.0, 1.0, 0.0, 0.0)" )') i
       endif
    end do  
   else
           write(IOinp,'(''MARKER_INLET= NONE'')') 
   endif
   write(IOinp,'(''%'')') 
  !1001 format(1x,'MARKER_INLET=(',i10,'288.6, 102010.0, 1.0, 0.0, 0.0)')
   
   write(IOinp,'(''% Outlet boundary marker(s) (NONE = no marker)'')') 
   write(IOinp,'(''% Format: ( outlet marker, back pressure (static), ... )'')') 
   if (Noutlet>0) then
     do i=1,MaxNames
       if (Names(i)=='OUTLET') then
           write(IOinp,'( "MARKER_OUTLET=(",i4,",101300.0)" )') i
       endif    
     end do
   else
           write(IOinp,'(''MARKER_OUTLET= NONE'')') 
   endif
   write(IOinp,'(''%'')') 
   
   write(IOinp,'(''% Marker(s) of the surface to be plotted or designed'')') 
   write(IOinp,'( ''MARKER_PLOTTING=(4,5)'' )') ! ( 5,4 )
   write(IOinp,'(''%'')') 
   
   write(IOinp,'(''% Marker(s) of the surface where the functional (Cd, Cl, etc.) will be evaluated'')') 
   write(IOinp,'( ''MARKER_MONITORING=(4,5)'' )') ! ( 5, 4 )
   write(IOinp,'('''')') 

  !do i=1,MaxNames
  !   if( Boundary(i) )then
  !     if( i <= 9 )then
  !       write(IOinp,'(''rname,'',i1,'','',A64)') i,Names(i)
  !     elseif( i <= 99 )then
  !       write(IOinp,'(''rname,'',i2,'','',A64)') i,Names(i)
  !     elseif( i <= 999 )then
  !       write(IOinp,'(''rname,'',i3,'','',A64)') i,Names(i)
  !     elseif( i <= 9999 )then
  !       write(IOinp,'(''rname,'',i4,'','',A64)') i,Names(i)
  !     elseif( i <= 99999 )then
  !       write(IOinp,'(''rname,'',i5,'','',A64)') i,Names(i)
  !     else
  !       write(IOinp,'(''rname,'',i6,'','',A64)') i,Names(i)
  !     endif
  !   endif
  ! end do

   close(IOinp)

   write(*,*) 'Done gmsh2su2'

   contains
   !------------------------------------------------------------------------------------
   subroutine check_input_character(c1,c2)

     implicit none

     character (len=*) :: c1, c2

     if( c1(1:len(c2)) /= c2 )then
       write(*,*)  'error reading Gmsh input file: ',&
                   'the following two characters should be the ',&
		   'same but differ ',c1(1:len(c2)),c2
       stop
     endif

   end subroutine
   subroutine check_element_type(ielement_type,element_type)

     implicit none
     integer ielement_type
     integer element_type(:)

     if( ielement_type < 0 )then
       write(*,*) 'error reading Gmsh file: element type must be positive'
       write(*,*) 'element type = ',ielement_type
       stop
     endif

     if( ielement_type > size(element_type) )then
       write(*,*) 'error reading Gmsh file: unrecognised element type'
       write(*,*) 'element type ',ielement_type
       write(*,*) 'max recognised element type ',size(element_type)
       stop
     endif

   end subroutine
   subroutine check_n_tags(ntags,itags)

     implicit none

     integer ntags
     integer itags(:)

     if( ntags > size(itags) )then
       write(*,*) 'error: The Gmsh file contains ',ntags,' tags per element'
       write(*,*) 'Gmsh2Su2 is hard-wired for a maximum of ',size(itags),&
        	  'tags.  The dimension of this array needs to be increased.'
       stop
     endif

   end subroutine
end


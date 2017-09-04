program su2_periodic

   !-----------------------------------------------------------------------------!
   !         Algorith to create periodic data structure (SU2 format)             ! 
   !-----------------------------------------------------------------------------!
   ! Author: Prof. Antonio Ghidoni                                               !
   ! email:  antonio.ghidoni@unibs.it                                            !
   ! University of Brescia, Department of Mechanical and Industrial Engineering  !
   !-----------------------------------------------------------------------------!

   implicit none

   character(30)           :: rmesh, wmesh
   integer(4), allocatable :: ne(:), np(:), nmark(:), nb(:,:)
   integer(4), allocatable :: neg(:), npg(:), nmarkg(:), nbg(:,:)
   integer(4)              :: nf, nd, nperio, id2np(15)

   ! Grid arrays
   real(8), allocatable       :: xp(:,:,:)
   real(8), allocatable       :: xpg(:,:,:)
   integer(4), allocatable    :: s2p(:,:,:), topo(:,:), b2p(:,:,:,:), npp2(:), &
                                 bp2bp(:,:,:), bperio(:,:), npp1(:), p2p(:,:), p3p4(:,:)
   character(20), allocatable :: cbc(:,:)

   ! Temporary grid arrays
   integer(4), allocatable    :: s2ptmp(:,:,:), b2ptmp(:,:,:,:), po2pn(:)
   real(8), allocatable       :: xptmp(:,:,:)

   ! Periodic bcs
   character(20)   :: perio(2,20)
   integer(4)      :: dperio(20)
   real(8)         :: lperio(9,20), ax_perio(3,20), lnper(20)
   logical         :: axial

   integer(4), parameter :: mtb = 20, mb = 100000, me = 8000000, mp = 8000000, &
                            md  = 3

   real(8), parameter :: tol0 = 1.e-8 

   real(8) :: scaling

   id2np = 0
   ! 2D
   id2np(3)  = 2
   id2np(5)  = 3
   id2np(9)  = 4
   ! 3D
   id2np(10) = 4
   id2np(12) = 8
   id2np(13) = 6
   id2np(14) = 5

   write(*,*)
   write(*,*)'===================================================='
   write(*,*)'READING cfg file'
   write(*,*)
   call read_cfg
   write(*,*)'===================================================='

   write(*,*)
   write(*,*)'===================================================='
   write(*,*)'READING mesh file: ',rmesh
   write(*,*)
   call read_mesh
   call write_tec
   write(*,*)'===================================================='

   write(*,*)
   write(*,*)'===================================================='
   write(*,*)'BUILD periodic data structure'
   write(*,*)
   call build_perio
   write(*,*)'===================================================='

   write(*,*)
   write(*,*)'===================================================='
   write(*,*)'WRITING mesh file: ',wmesh
   call write_mesh
   call write_tec_perio
   write(*,*)'===================================================='

contains

   subroutine read_cfg

      character(200)  :: str200
      integer(4)      :: i, l, idf, m(12), mm
      logical         :: rmsh, wmsh
 

      scaling = 1.0
      nperio = 0
      rmsh = .false.
      wmsh = .false.
      do l = 1, 10000
         do i = 1, 200
            str200(i:i) = ' '
         enddo

         read (*,'(a200)',end=9) str200

         if ( str200(1:16) == 'MARKER_PERIODIC=' ) then
            nperio = nperio +1
            mm  = 0
            m   = 0
            do i = 1, 200
               if ( str200(i:i) == '(' ) then
                  mm = mm +1
                  m(mm)  = i+1
               elseif ( str200(i:i) == ',' ) then
                  mm = mm +1
                  m(mm)  = i-1
               elseif ( str200(i:i) == ')' ) then
                  mm = mm +1
                  m(mm)  = i-1
               endif
            enddo
            read (str200(m(1):m(2)),*)perio(1,nperio)
            read (str200(m(2)+2:m(3)),*) perio(2,nperio)
            read (str200(m(3)+2:m(4)),*)lperio(1,nperio) 
            read (str200(m(4)+2:m(5)),*)lperio(2,nperio) 
            read (str200(m(5)+2:m(6)),*)lperio(3,nperio) 
            read (str200(m(6)+2:m(7)),*)lperio(4,nperio) 
            read (str200(m(7)+2:m(8)),*)lperio(5,nperio) 
            read (str200(m(8)+2:m(9)),*)lperio(6,nperio) 
            read (str200(m(9)+2:m(10)),*)lperio(7,nperio) 
            read (str200(m(10)+2:m(11)),*)lperio(8,nperio) 
            read (str200(m(11)+2:m(12)),*)lperio(9,nperio) 

            write(*,*)'Periodic condition:',nperio
            write(*,'(3a25)')'FLAG=',perio(1,nperio), perio(2,nperio)
            write(*,*)'Rotation_center_(x,y,z)'
            write(*,'(3f15.8)')lperio(1:3,nperio)
            write(*,*)'Rotation_angle_(x,y,z)-axis'
            write(*,'(3f15.8)')lperio(4:6,nperio)
            write(*,*)'Translation_(x,y,z)'
            write(*,'(3f15.8)')lperio(7:9,nperio)
         elseif ( str200(1:14) == 'MESH_FILENAME=' ) then
            read (str200(15:34),*)rmesh
            rmsh = .true.
         elseif ( str200(1:18) == 'MESH_OUT_FILENAME=' ) then
            read (str200(19:38),*)wmesh
            wmsh = .true.
         elseif ( str200(1:16) == 'MSH_PERIO_SCALE=' ) then
            read (str200(17:38),*)scaling
            wmsh = .true.
         endif
      enddo

   9  continue

      if ( .not. rmsh ) then
         write(*,*)'ERROR! INPUT mesh file not specified!'
         stop
      endif
      if ( .not. wmsh ) then
         write(*,*)'WARNING! OUTPUT mesh file not specified!'
         write(*,*)'Default name will be used: mesh_out.su2'
         wmesh = 'mesh_out.su2'
         stop
      endif

      if ( lperio(4,1) == 0. .and. lperio(5,1) == 0. .and. lperio(6,1) == 0. ) then
         axial = .true.
      else
         axial = .false.
      endif

      do i = 1, nperio
         if ( lperio(4,i) /= 0. ) then
            dperio(i)     = 1
            lnper(i)     = lperio(4,i)*acos(-1.0)/180.
            ax_perio(1,i) = 1.
            ax_perio(2,i) = 0.
            ax_perio(3,i) = 0.
         elseif ( lperio(5,i) /= 0. ) then
            dperio(i)     = 2
            lnper(i)     = lperio(5,i)*acos(-1.0)/180.
            ax_perio(1,i) = 0.
            ax_perio(2,i) = 1.
            ax_perio(3,i) = 0.
         elseif ( lperio(6,i) /= 0. ) then
            dperio(i)     = 3
            lnper(i)     = lperio(6,i)*acos(-1.0)/180.
            ax_perio(1,i) = 0.
            ax_perio(2,i) = 0.
            ax_perio(3,i) = 1.
         endif
      enddo

   end subroutine read_cfg

   subroutine read_mesh

      integer(4)     :: iunit, i, j, k, io, last, m, npri, npyr, ntet, nhex, ntri, nqua
      character(6)   :: c6, control
      character(200) :: str200
      real(8)        :: ct(3,4) 

      
      iunit = 15
      open(unit=iunit,file=rmesh,form="formatted")

      ! nf: number of zones
      control = 'NZONE='
      call tstlin2(iunit,control,io)
      if ( io == 0 ) then
         read(iunit,'(a6,i5)') c6, nf
      else
         nf = 1
      endif


      allocate (ne(nf), np(nf), nmark(nf))

      allocate (xp(md,mp,nf))
      allocate (p2p(mp,nf))
      allocate (p3p4(mp,nf))
      allocate (s2p(8,me,nf))
      allocate (topo(me,nf))

      allocate (nb(mtb,nf))
      allocate (b2p(5,mb,mtb,nf))
      allocate (cbc(mtb,nf))

      allocate (xptmp(md,mp,nf))
      allocate (s2ptmp(8,me,nf))
      allocate (b2ptmp(5,mb,mtb,nf))
      allocate ( po2pn(mp) )

      do i = 1, nf

         rewind(iunit)
         control = 'NDIME='
         call tstlin2(iunit,control,io)
         if ( io == 0 ) then
            read(iunit,'(a6,i5)') c6, nd
         else
            write(*,*)'Error! Cannot read mesh dimension'
            stop
         endif

         ! ne(:): number of elements for each zone
         rewind(iunit)
         control = 'NELEM='
         call tstlin2(iunit,control,io)
         if ( io == 0 ) then
            read(iunit,'(a6,i15)') c6, ne(i)
         else
            write(*,*)'Error! Cannot read mesh elements'
            stop
         endif

         ! topo(:,:): element ID
         ! s2p(:,:,:): connectivity of each element
         do j = 1, ne(i)
            read(iunit,*)topo(j,i), s2ptmp(1:id2np(topo(j,i)),j,i)
         enddo

         ! np(:): number of points for each zone
         rewind(iunit)
         control = 'NPOIN='
         call tstlin2(iunit,control,io)
         if ( io == 0 ) then
            read(iunit,'(a6,i15)') c6, np(i)
         else
            write(*,*)'Error! Cannot read points coordinates'
            stop
         endif

         ! xp(:,:,:): points coordinates for each zone
         do j = 1, np(i)
            read(iunit,*)xptmp(1:nd,j,i)
         enddo
         p2p(:,i) = -1
         p3p4(:,i) = -1

         rewind(iunit)
         control = 'NMARK='
         call tstlin2(iunit,control,io)

         ! nmark(:): number of bcs for each zone
         if ( io == 0 ) then
            read(iunit,'(a6,i15)') c6, nmark(i)
         else
            write(*,*)'Error! Cannot read number of bcs'
            stop
         endif
        
         ! cbc(:,:): ID of the bcs for each zone
         ! nb(:,:): number of boundary elements for each bcs 
         ! b2p(:,:,:,:): connectivity of every boundary element for each bcs 
         do j = 1, nmark(i)
            read (iunit,'(a200)') str200
            read (str200(12:31),*)cbc(j,i)
            read (iunit,'(a200)') str200
            read (str200(14:28),*)nb(j,i)
            do k = 1, nb(j,i)
               read(iunit,*)b2ptmp(1,k,j,i), b2ptmp(2:id2np(b2ptmp(1,k,j,i))+1,k,j,i)
            enddo
         enddo
      enddo

      close(iunit)

      ! Change numbering

      do i = 1, nf

         po2pn = -1

         do j = 1, nmark(i)
            if ( cbc(j,i) == perio(1,i) ) then
               last = np(i)-1
               do k = 1, nb(j,i)
                  do m = 2, id2np(b2ptmp(1,k,j,i))+1
                     if ( po2pn(b2ptmp(m,k,j,i)+1) == -1 ) then
                        po2pn(b2ptmp(m,k,j,i)+1) = last
                        last = last -1
                     endif 
                  enddo
               enddo
            endif
         enddo

         last = 0
         do j = 1, ne(i)
            do k = 1, id2np(topo(j,i)) 
               if ( po2pn(s2ptmp(k,j,i)+1) == -1 ) then
                  po2pn(s2ptmp(k,j,i)+1) = last
                  last = last +1
               endif
            enddo 
         enddo

         do j = 1, ne(i)
            do k = 1, id2np(topo(j,i)) 
               s2p(k,j,i) = po2pn(s2ptmp(k,j,i)+1)
               xp(:,s2p(k,j,i)+1,i) = xptmp(:,s2ptmp(k,j,i)+1,i)
            enddo 
         enddo

         do j = 1, nmark(i)
            do k = 1, nb(j,i)
               b2p(1,k,j,i) = b2ptmp(1,k,j,i)
               do m = 2, id2np(b2ptmp(1,k,j,i))+1
                  b2p(m,k,j,i) = po2pn(b2ptmp(m,k,j,i)+1)
               enddo
            enddo
         enddo

         npri = 0
         ntet = 0
         nhex = 0
         npyr = 0
         ntri = 0
         nqua = 0
         do j = 1, ne(i)
            if ( topo(j,i) == 5 ) ntri = ntri +1
            if ( topo(j,i) == 9 ) nqua = nqua +1
            if ( topo(j,i) == 10 ) ntet = ntet +1
            if ( topo(j,i) == 12 ) nhex = nhex +1
            if ( topo(j,i) == 13 ) npri = npri +1
            if ( topo(j,i) == 14 ) npyr = npyr +1
         enddo

        
         write(*,*)'Domain=',i
         if ( nd == 2 ) then
            write(*,'(a12,i15)')'Triangles=', ntri
            write(*,'(a12,i15)')'Quads=', nqua
         else
            write(*,'(a12,i15)')'Tetrahedra=', ntet
            write(*,'(a12,i15)')'Prysm=', npri
            write(*,'(a12,i15)')'Pyramid=', npyr
            write(*,'(a12,i15)')'Hexahedra=', nhex
         endif
         write(*,*)
         write(*,*)'Boundary conditions'
         do j = 1, nmark(i)
            write(*,*)cbc(j,i),'=',nb(j,i)
         enddo 
      enddo


   end subroutine read_mesh

   subroutine build_perio

      integer(4) :: i, j, k, ii, iii, e1, e2, m, n, kkk, w, npri, ntet, nhex, npyr, ntri, nqua
      integer(4) :: jj, jjj, kk, coun
      logical    :: found, exist
      real(8)    :: xa(3), xb(3), r1, r2, dx(3), da, db, dc, alp1, alp2, &
                    delta, delteta, oo(3), xap(3), tol


      allocate( bperio(2,nf) )
      bperio = 0

      ! bperio(:,:): numbering of two periodic bc for every zone
      do j = 1, nf
         do k = 1, nmark(j)
            if ( cbc(k,j) == perio(1,j) ) then
               bperio(1,j) = k
            elseif ( cbc(k,j) == perio(2,j) ) then
               bperio(2,j) = k
            endif
         enddo
      enddo       

      allocate( bp2bp(4,mb+1,nf) )
      allocate( npp1(nf) )

      bp2bp    = 0
      npp1    = 0

      ! bp2bp(1,:,:): list of points associated to side e1
      ! npp1(:): number of periodic points on the edge e1
      do j = 1, nf
         e1 = bperio(1,j)
               
         do i = 1, nb(e1,j)
            do k = 2, id2np(b2p(1,i,e1,j))+1
               found = .false.
               do ii = 1, npp1(j)
                  if ( bp2bp(1,ii,j) == b2p(k,i,e1,j) ) then
                     found = .true.
                     exit
                  endif
               enddo
               if ( .not. found ) then
                  npp1(j) = npp1(j) +1
                  bp2bp(1,npp1(j),j) = b2p(k,i,e1,j)
               endif 
            enddo
         enddo

!         do i = 1, npp1(j)
!            write(555,*)xp(1:nd,bp2bp(1,i,j)+1,j)
!         enddo
      enddo

      ! bp2bp(2,:,:): list of points associated to side e2
      ! p2p(:,:): for each periodic point of e1, the number of the corresponding point
      !           on e2 is given

      xa = 0.
      xb = 0. 


      do j = 1, nf

         tol = tol0

         m = 0
         n = 0
         do kkk = 1, 3
            if ( dperio(j) == kkk ) cycle
            if ( m == 0 ) then
               m = kkk
            else
               n = kkk
            endif
         enddo

         e1 = bperio(1,j)
         e2 = bperio(2,j) 

         do i = 1, npp1(j)
            xa(1:nd) = xp(:,bp2bp(1,i,j)+1,j)

            if ( axial ) then
               xap = xa +lperio(7:9,j)
            else
               oo = lperio(1:nd,j)
               oo(dperio(j)) = xa(dperio(j))
               dx    = xa -oo
               r1    = sqrt(dx(1)**2 +dx(2)**2 +dx(3)**2)
               alp1  = atan2(dx(n),dx(m))
   
               xap(1:nd) = xa
               xap(m)    = lperio(m,j) +r1*cos(alp1 +lnper(j))
               xap(n)    = lperio(n,j) +r1*sin(alp1 +lnper(j))
            endif

     100    found = .false.
            do k = 1, nb(e2,j)
               do ii = 2, id2np(b2p(1,k,e2,j))+1
                  xb(1:nd) = xp(:,b2p(ii,k,e2,j)+1,j) 

                  da = sqrt((xap(1)-xb(1))**2+(xap(2)-xb(2))**2+(xap(3)-xb(3))**2)

                  if ( da < tol ) then
                     write(*,*)i, da
                     bp2bp(2,i,j) = b2p(ii,k,e2,j)
                     p2p(bp2bp(1,i,j)+1,j) = bp2bp(2,i,j)
                     found = .true.
                     exit
                  endif
               enddo
               if ( found ) exit
            enddo
            if ( .not. found .and. tol <= 1.e-3) then
               tol = tol*10.
               goto 100
            elseif ( .not. found) then
               write(*,*)'ERROR! Cannot find periodic point'
               stop
            endif
         enddo
!         do i = 1, npp1(j)
!            write(777,*)xp(1:nd,bp2bp(2,i,j)+1,j)
!         enddo
      enddo

      allocate ( neg(nf), npg(nf), nmarkg(nf), nbg(mtb,nf) )
      allocate( npp2(nf) )
      allocate (xpg(md,mp,nf))

      npp2 = 0

      xpg = 0.0
      do i = 1, nf
         do j = 1, np(i)
            xpg(:,j,i) = xp(:,j,i)
         enddo
      enddo

      ! nmarkg(:): number of bcs for each zone with halos
      do i = 1, nf
         nmarkg(i) = nmark(i) +2
      enddo

      do i = 1, nf
         do j = nmark(i)+1, nmarkg(i)
            cbc(j,i)="SEND_RECEIVE"
         enddo
      enddo
     
      ! neg(:): number of elements for each zone (with halos)
      ! npg(:): number of points for each zone (with halos)
      ! bp2bp(3,:,:): list of points (inside domain) for the elements sharing the side e1
      ! bp2bp(4,:,:): list of points for the elements sharing the side e1
      do i = 1, nf
         neg(i) = ne(i)
         npg(i) = np(i)
         do j = 1, ne(i)
            found = .false.
            do ii = 1, id2np(topo(j,i))
               if ( p2p(s2p(ii,j,i)+1,i) >= 0 ) then
                  neg(i) = neg(i) +1
                  topo(neg(i),i)     = topo(j,i)
                  found = .true.
                  exit
               endif
            enddo
            if ( found ) then
               do ii = 1, id2np(topo(j,i))
                  if ( p2p(s2p(ii,j,i)+1,i) >= 0 ) then
                     s2p(ii,neg(i),i) = p2p(s2p(ii,j,i)+1,i)
                  else
                     exist = .false.
                     do iii = 1, npp2(i)
                        if ( bp2bp(3,iii,i) == s2p(ii,j,i) ) then
                           exist = .true.
                           exit
                        endif
                     enddo
                     if ( .not. exist ) then 
                        npp2(i)            = npp2(i) +1
                        npg(i)             = npg(i) +1
                        s2p(ii,neg(i),i)   = npg(i) -1
                        bp2bp(3,npp2(i),i) = s2p(ii,j,i) 
                        bp2bp(4,npp2(i),i) = s2p(ii,neg(i),i) 
                        if ( axial ) then
                           xpg(1:nd,npg(i),i) = xp(1:nd,s2p(ii,j,i)+1,i) +lperio(7:7+nd-1,i)
                        else
                           m = 0
                           n = 0
                           do kkk = 1, 3
                              if ( dperio(i) == kkk ) cycle
                              if ( m == 0 ) then
                                 m = kkk
                              else
                                 n = kkk
                              endif
                           enddo

                           oo = lperio(1:nd,i)
                           oo(dperio(i)) = xp(dperio(i),s2p(ii,j,i)+1,i) 
                           dx    = xp(1:nd,s2p(ii,j,i)+1,i) -oo 
                           r1    = sqrt(dx(1)**2 +dx(2)**2 +dx(3)**2)
                           alp1  = atan2(dx(n),dx(m))

                           xpg(1:nd,npg(i),i) = xp(1:nd,s2p(ii,j,i)+1,i)
                           xpg(m,npg(i),i)    = lperio(m,i) +r1*cos(alp1 +lnper(i))
                           xpg(n,npg(i),i)    = lperio(n,i) +r1*sin(alp1 +lnper(i))
                        endif
                     else
                        s2p(ii,neg(i),i) = bp2bp(4,iii,i)
                     endif
                  endif
               enddo
            endif
         enddo
         npri = 0
         ntet = 0
         nhex = 0
         npyr = 0
         ntri = 0
         nqua = 0
         do j = 1, neg(i)
            if ( topo(j,i) == 5 ) ntri = ntri +1
            if ( topo(j,i) == 9 ) nqua = nqua +1
            if ( topo(j,i) == 10 ) ntet = ntet +1
            if ( topo(j,i) == 12 ) nhex = nhex +1
            if ( topo(j,i) == 13 ) npri = npri +1
            if ( topo(j,i) == 14 ) npyr = npyr +1
         enddo


         write(*,*)'Domain=',i
         if ( nd == 2 ) then
            write(*,'(a12,i15)')'Triangles=', ntri
            write(*,'(a12,i15)')'Quads=', nqua
         else
            write(*,'(a12,i15)')'Tetrahedra=', ntet
            write(*,'(a12,i15)')'Prysm=', npri
            write(*,'(a12,i15)')'Pyramid=', npyr
            write(*,'(a12,i15)')'Hexahedra=', nhex
         endif
      enddo

      do i = 1, nf
         do j = 1, npp2(i)
            p3p4(bp2bp(3,j,i)+1,i) = bp2bp(4,j,i)
!            write(888,*)xpg(1:nd,bp2bp(3,j,i)+1,i)
!            write(999,*)xpg(1:nd,bp2bp(4,j,i)+1,i)
         enddo
      enddo

      nbg = nb
      do i = 1, nf
         do j = 1, nmark(i)
            if ( j == bperio(1,i) .or. j == bperio(2,i) ) cycle
!
!Bisogna verififcare quali boundary intersecano il periodico
!            do w = 1, nwall
!               if ( cbc(j,i) == cwall(w,i) ) goto 100   
!            enddo

            do k = 1, nb(j,i)
               iii = 0
               jjj = 0
               do ii = 2, id2np(b2p(1,k,j,i))+1
                  if ( p2p(b2p(ii,k,j,i)+1,i) >= 0 ) jjj = jjj +1
                  if ( p3p4(b2p(ii,k,j,i)+1,i) >= 0 ) iii = iii +1
               enddo
               if ( (iii+jjj) == id2np(b2p(1,k,j,i)) .and. jjj > 0) then
                  nbg(j,i) = nbg(j,i) +1
                  do ii = 2, id2np(b2p(1,k,j,i))+1
                     if ( p2p(b2p(ii,k,j,i)+1,i) >= 0 ) then
                        b2p(ii,nbg(j,i),j,i) = p2p(b2p(ii,k,j,i)+1,i)
                     elseif ( p3p4(b2p(ii,k,j,i)+1,i) >= 0 ) then 
                        b2p(ii,nbg(j,i),j,i) = p3p4(b2p(ii,k,j,i)+1,i)
                     endif
                  enddo
                  b2p(1,nbg(j,i),j,i) = b2p(1,k,j,i)
               endif
            enddo
!       100  continue
         enddo
         write(*,*)
         write(*,*)'Boundary conditions'
         do j = 1, nmark(i)
            write(*,*)cbc(j,i),'=',nbg(j,i)
         enddo


!!!!         do j = 1, nmark(i)
!!!!            do jj = 1, nbg(j,i)
!!!!               do k = 1, neg(i)
!!!!                  coun = 0
!!!!                  found = .false.
!!!!                  do jjj = 2, id2np(b2p(1,jj,j,i))+1
!!!!                     do kk = 1, id2np(topo(k,i))
!!!!                        if ( b2p(jjj,jj,j,i) == s2p(kk,k,i) ) coun = coun +1
!!!!                     enddo
!!!!                     if ( coun == id2np(b2p(1,jj,j,i)) ) then
!!!!                        found = .true.
!!!!                        exit
!!!!                     endif
!!!!                  enddo
!!!!                  if ( found ) exit
!!!!               enddo
!!!!               if ( .not. found ) then
!!!!                  write(*,*)'ERROR! Face without corresponding element'
!!!!                  write(*,*)b2p(1:id2np(b2p(1,jj,j,i))+1,jj,j,i)
!!!!                  write(*,*)xpg(1:nd,b2p(2,jj,j,i),i)
!!!!                  write(*,*)xpg(1:nd,b2p(3,jj,j,i),i)
!!!!                  write(*,*)xpg(1:nd,b2p(4,jj,j,i),i)
!!!!                  stop
!!!!               endif
!!!!            enddo
!!!!         enddo
      enddo

   end subroutine build_perio

   subroutine write_tec

      integer(4) :: i, j, idf


      idf = 10
      open(unit=idf,file='periodic_original.dat')

      if ( nd == 2 ) then
         do i = 1, nf
            write (idf,*) 'VARIABLES = "X",  "Y"'
            write (idf,*) 'ZONE N=',np(i),', E=',ne(i),', F=FEPOINT, ET=QUADRILATERAL'
            do j = 1, np(i)
               write (idf,'(2f20.10)') xp(1:nd,j,i)
            enddo
            write(idf,*)
            do j = 1, ne(i)
               if ( topo(j,i) == 5 ) then
                  write (idf,'(4i15)') s2p(1:3,j,i)+1, s2p(1,j,i)+1
               else
                  write (idf,'(4i15)') s2p(1:4,j,i)+1
               endif
            enddo
         enddo
      else
         do i = 1, nf
            write (idf,*) 'VARIABLES = "X",  "Y",  "Z"'
            write (idf,*) 'ZONE N=',np(i),', E=',ne(i),', & 
                           F=FEPOINT, ET=BRICK'
            do j = 1, np(i)
               write (idf,'(2f20.10)') xp(:,j,i)
            enddo
            write(idf,*)
            do j = 1, ne(i)
               if ( topo(j,i) == 10 ) then
                  ! TET
                  write (idf,'(8i10)') s2p(1:3,j,i)+1, s2p(3,j,i)+1, s2p(4,j,i)+1, s2p(4,j,i)+1, s2p(4,j,i)+1, s2p(4,j,i)+1
               elseif ( topo(j,i) == 12 ) then
                  ! HEX
                  write (idf,'(8i10)') s2p(1:8,j,i)+1
               elseif ( topo(j,i) == 13 ) then
                  ! PRISM
                  write (idf,'(8i10)') s2p(1:3,j,i)+1, s2p(3,j,i)+1, s2p(4:6,j,i)+1, s2p(6,j,i)+1
               elseif ( topo(j,i) == 14 ) then
                  ! PYRAMID
                  write (idf,'(8i10)') s2p(1:4,j,i)+1, s2p(5,j,i)+1, s2p(5,j,i)+1, s2p(5,j,i)+1, s2p(5,j,i)+1
               endif
            enddo
         enddo
      endif

      close(idf)

   end subroutine write_tec

   subroutine write_tec_perio

      integer(4) :: i, j, idf


      idf = 10
      open(unit=idf,file='periodic_halo.dat')

      if ( nd == 2 ) then
         do i = 1, nf
            write (idf,*) 'VARIABLES = "X",  "Y"'
            write (idf,*) 'ZONE N=',npg(i),', E=',neg(i),', F=FEPOINT, ET=QUADRILATERAL'
            do j = 1, npg(i)
               write (idf,'(2f20.10)') xpg(1:nd,j,i)
            enddo
            write(idf,*)
            do j = 1, neg(i)
               if ( topo(j,i) == 5 ) then
                  write (idf,'(4i15)') s2p(1:3,j,i)+1, s2p(1,j,i)+1
               else
                  write (idf,'(4i15)') s2p(1:4,j,i)+1
               endif
            enddo
         enddo
      else
         do i = 1, nf
            write (idf,*) 'VARIABLES = "X",  "Y",  "Z"'
            write (idf,*) 'ZONE NODES=',npg(i),', ELEMENTS=',neg(i),', & 
                           DATAPACKING=POINT, ZONETYPE=FEBRICK'
            do j = 1, npg(i)
               write (idf,'(3f20.10)') xpg(:,j,i)
            enddo
            write(idf,*)
            do j = 1, neg(i)
               if ( topo(j,i) == 10 ) then
                  ! TET
                  write (idf,'(8i15)') s2p(1:3,j,i)+1, s2p(3,j,i)+1, s2p(4,j,i)+1, s2p(4,j,i)+1, s2p(4,j,i)+1, s2p(4,j,i)+1
               elseif ( topo(j,i) == 12 ) then
                  ! HEX
                  write (idf,'(8i15)') s2p(1:8,j,i)+1
               elseif ( topo(j,i) == 13 ) then
                  ! PRISM
                  write (idf,'(8i15)') s2p(1:2,j,i)+1, s2p(2,j,i)+1,s2p(3,j,i)+1, s2p(4:5,j,i)+1, s2p(5,j,i)+1,s2p(6,j,i)+1
               elseif ( topo(j,i) == 14 ) then
                  ! PYRAMID
                  write (idf,'(8i15)') s2p(1:4,j,i)+1, s2p(5,j,i)+1, s2p(5,j,i)+1, s2p(5,j,i)+1, s2p(5,j,i)+1
               endif
            enddo
         enddo
      endif

      close(idf)

   end subroutine write_tec_perio

   subroutine write_mesh

      integer(4)   :: idf, if, i, j, ip, iel
      character(2) :: sz

      write(*,*)'Scaling=',scaling
      idf = 25
      open(unit=idf,file=wmesh,form="formatted")

      write(idf,'(a6,i3)') 'NZONE=', nf
      
      do if = 1, nf

      write(idf,*)
      write(idf,'(a6,i3)') 'IZONE=', if

      write(idf,*) '%'
      write(idf,*) '% Problem dimension'
      write(idf,*) '%'
      write(idf,'(a6,i3)') 'NDIME=', nd

      write(idf,*) '%'
      write(idf,*) '% Inner element connectivity'
      write(idf,*) '%'

      write(idf,'(a6,i10)') 'NELEM=', neg(if)

      ! Write the interior element connectivity...

      iel = 0
      do i = 1, neg(if)
         sz = int2ch2( id2np(topo(i,if)) )
         write (idf,'(1i3,'//sz//'i10,i10)') topo(i,if), ( s2p(j,i,if), j=1,id2np(topo(i,if)) ), iel
         iel = iel +1
      enddo

      write(idf,*) '%'
      write(idf,*) '% Node coordinates'
      write(idf,*) '%'
      write(idf,'(a6,2i10)') 'NPOIN=', npg(if), npg(if)-npp1(if)-npp2(if) 

      do ip = 1, npg(if)
         sz = int2ch2( nd )
         write (idf,'('//sz//'f25.15,i10)') scaling*xpg(1:nd,ip,if), ip-1
      enddo

      write(idf,'(a6,i3)') 'NMARK=', nmarkg(if)

      do j = 1, nmark(if)
         write(idf,'(a11,a10)') 'MARKER_TAG=', cbc(j,if)
         write(idf,'(a13,i10)') 'MARKER_ELEMS=', nbg(j,if) 
         do i = 1, nbg(j,if)
            sz = int2ch2( id2np(b2p(1,i,j,if)))
            write(idf,'(1i3,'//sz//'i10)') b2p(1,i,j,if), b2p(2:id2np(b2p(1,i,j,if))+1,i,j,if)
         enddo
      enddo

      j = nmark(if)+1
      write(idf,'(a11,a20)') 'MARKER_TAG=', cbc(j,if)
      write(idf,'(a13,i10)') 'MARKER_ELEMS=', npp1(if) +npp2(if) 
      write(idf,'(a10)') 'SEND_TO= 1'
      do i = 1, npp1(if)
         write(idf,*)1, bp2bp(2,i,if), 1
      enddo
      do i = 1, npp2(if)
         write(idf,*)1, bp2bp(3,i,if), 2
      enddo

      j = nmarkg(if)
      write(idf,'(a11,a20)') 'MARKER_TAG=', cbc(j,if)
      write(idf,'(a13,i10)') 'MARKER_ELEMS=',  npp1(if) +npp2(if) 
      write(idf,'(a10)') 'SEND_TO=-1'
      do i = 1, npp1(if)
         write(idf,*)1, bp2bp(1,i,if), 1
      enddo
      do i = 1, npp2(if)
         write(idf,*)1, bp2bp(4,i,if), 2
      enddo

      write(idf,'(a12)') 'NPERIODIC= 3'
      write(idf,'(a15,i4)') 'PERIODIC_INDEX=',0
      write(idf,'(3f15.9)') lperio(1:3,if)
      write(idf,'(3f15.9)') lperio(1:3,if)
      write(idf,'(3f15.9)') lperio(1:3,if)
      write(idf,'(a15,i4)') 'PERIODIC_INDEX=',1
      write(idf,'(3f15.9)') lperio(1:3,if)
      write(idf,'(3f15.9)') lperio(4:6,if)*acos(-1.0)/180.
      write(idf,'(3f15.9)') scaling*lperio(7:9,if)
      write(idf,'(a15,i4)') 'PERIODIC_INDEX=',2
      write(idf,'(3f15.9)')-lperio(1:3,if)
      write(idf,'(3f15.9)')-lperio(4:6,if)*acos(-1.0)/180.
      write(idf,'(3f15.9)')-scaling*lperio(7:9,if)

      enddo

      close(idf)

   end subroutine write_mesh

   subroutine tstlin2 ( iunit, string, io )

      integer(4)   :: iunit, io
      character(6) :: string

      character(80) :: text

 2002 read(iunit,'(A80)',IOSTAT=io)text

      if ( io == 0 ) then
         if ( text(1:6) /= string ) goto 2002

         backspace (iunit)
      endif

   end subroutine tstlin2

   function int2ch2 ( n ) result ( c2 )

      integer(4),   intent(in) :: n
      character(2)             :: c2

      integer(4) :: m, n0, n1


      m       =  n
      n1      =  m/10
      n0      =  m -10*n1

      m       =  ichar('0')
      c2      =  char(m+n1)//char(m+n0)

   end function int2ch2

   function voltet( ct ) result ( vol )

      ! tetrahedron volume

      real(8) :: ct(3,4), v(3,4), vol

      integer(4) :: j, i


      do j=2,4
         do i=1,3
            v(i,j) = ct(i,j) - ct(i,1)
         enddo
      enddo

      vol = v(1,4)*(v(2,2)*v(3,3)-v(3,2)*v(2,3))+ &
               v(2,4)*(v(3,2)*v(1,3)-v(1,2)*v(3,3))+ &
               v(3,4)*(v(1,2)*v(2,3)-v(2,2)*v(1,3))
      vol =(1./6.)*vol

   end function voltet

end program su2_periodic

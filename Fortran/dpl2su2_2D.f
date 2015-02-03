      PROGRAM dpl2su2

C    Code provided by Dr. Aldo Bonfiglioli
C    Associate professor of Fluid Flow Machinery
C    Scuola di Ingegneria
C    Universita' della Basilicata

C
C    This routine reads a .dpl unstructured datafile
C    and creates a SU2 mesh file
C
C    dplot is a:
C     write (*,*) 'STRUCTURED/UNSTRUCTURED CONTOUR PLOTTING PACKAGE'
C     write (*,*) '    -Written by  DARREN DE ZEEUW                '
C     write (*,*) '    -Last Modified 02/93                        '
C     write (*,*) ' '
C
C    the dplot format is described in format 110-130
C
C    Unstructured dplot file are created using the
C    Delaundo mesh generator
C    http://www.cerfacs.fr/~muller/delaundo.html 
C
      IMPLICIT NONE
C
      INTEGER NELEMAX,NODEMAX,NFACMAX,NDIM,NOFVERT
      PARAMETER (NELEMAX=525000,NODEMAX=265000,
     &           NFACMAX=5000,NDIM=2,NOFVERT=NDIM+1)
C
      INTEGER nCells,nNodes,nCorners,nBoundaryFaces,NFAC,nBodies,
     +        nBodyFaces,iBodyType
      INTEGER node1,node2,ntriangle
      INTEGER i,j,k,npo
      INTEGER NELEM,NPOIN,NBFAC,NHOLE
C
      INTEGER Corner(3,NELEMAX)
C
      DOUBLE PRECISION XY(2,NODEMAX)
      CHARACTER*80 filename,title
C
      EQUIVALENCE (nCells,NELEM)
      EQUIVALENCE (nNodes,NPOIN)
      EQUIVALENCE (NFAC,NBFAC)
C
      WRITE (6,*)
      WRITE (6,*) '    Give the input filename in dpl format :'
      WRITE (6,*)
      READ (5,'(A80)') filename
C
      WRITE(6,150)
    5 OPEN (1,FILE=filename,STATUS='OLD',FORM='formatted')
      OPEN(UNIT=12,FILE="MESH.su2")
C
C    UNSTRUCTURED DATA FILE
C
      READ (1,'(A80)') title
      WRITE (6,'(A80)') title
c
c..Needs "unstr" or "UNSTR" as first five characters
c
      READ (1,*) nCells
      IF (nCells.GT.NELEMAX) THEN
         WRITE(6,*)'Actual nof cells ',nCells,' > Max allowed = ',
     &   NELEMAX 
         STOP 'Increase NELEMAX within dpl2su2' 
      ENDIF
      DO 19 i = 1,nCells
          READ (1,*) nCorners, (Corner(j,i),j=1,nCorners)
   19 CONTINUE
      WRITE (6,*) nCells,'    Cells read !'
      WRITE(12,FMT=201)NDIM
      WRITE(12,FMT=202)NELEM
      do I= 1,NELEM
         WRITE(12,FMT=80)5,((Corner(j,I)-1),j=1,NOFVERT),I-1
      enddo
      WRITE (6,*) nCells,'    Cells have been written !'
C
      READ (1,*) nNodes
      IF (nNodes.GT.NODEMAX) THEN
         WRITE(6,*)'Actual nof nodes ',nNodes,' > Max allowed = ',
     &   NODEMAX 
         STOP 'Increase NODEMAX within dpl2su2'
      ENDIF
      READ (1,*)
c
c..read gridpoints
c
      DO 21 i = 1,nNodes
          READ (1,*) XY(1,i),XY(2,i)
   21 CONTINUE
      WRITE (6,*) nNodes,'    Nodes read !'
C
      WRITE(12,FMT=203) NPOIN
      do 125 I=1,NPOIN
          WRITE(12,FMT=*)(XY(j,I),j=1,NDIM),I-1
  125 CONTINUE
      WRITE (6,*) nNodes,'    Nodes written !'
C
      READ (1,*) nBodies
      Write(6,*)nBodies,' have been found in the dpl file'
      WRITE(12,FMT=204)nBodies
      DO 25 j = 1,nBodies
          READ (1,*) nBodyFaces,iBodyType
          WRITE (6,*) 'There are ',nBodyFaces,' edges on bndry # ',j,
     &' type is ',iBodyType
         write(12,FMT=205)iBodyType
         write(12,FMT=206)nBodyFaces
         DO 25 i = 1,nBodyFaces
              READ (1,*) node1,node2,ntriangle
              write(12,FMT=95)3,node1-1,node2-1
   25 CONTINUE
c
c  REM: dplot files created by Delaundo have nBoundaryFaces=0
c
c
      READ (1,*) nBoundaryFaces
      IF (nBoundaryFaces.NE.0) THEN
          WRITE (*,*) '    nBoundaryFaces is NOT 0 but ... ',
     +      nBoundaryFaces
          CALL EXIT(1)

      ENDIF
      CLOSE (1)
C
C
      WRITE(6,*)
      WRITE(6,*)'REMEMBER:'
      WRITE(6,*)
      WRITE(6,*)'Give a meaningful name to MARKER_TAG in the mesh file'
      WRITE(6,*)'It needs to have a corresponding boundary condition'
      WRITE(6,*)'in the configuration file!!!!'
C
C
C
   80 FORMAT(I1,1X,5(I6,1X))
   95 FORMAT(I1,1X,4(I6,1X))
 110  format ('                                                       '/
     +        'UNSTRUCTURED DATA FILE                                 '/
     +        '   write (1,"(a22)") "unstructured grid data"          '/
     +        'c..Needs "unstr" or "UNSTR" as first five characters   '/
     +        '   write (1,*) nCells                                  '/
     +        '   do i=1,nCells                                       '/
     +        '      write (1,*) nCorners,(Corner(j,i),j=1,nCorners)  '/
     +        '   end do                                              ')
 120  format ('   write (1,*) nNodes                                  '/
     +        '   write (1,*) U1(inf),U2(inf),U3(inf),U4(inf)         '/
     +        'c..Freestream state values                             '/
     +        '   do i=1,nNodes                                       '/
     +        '      write (1,*) X(i),Y(i),U1(i),U2(i),U3(i),U4(i)    '/
     +        '   end do                                              ')
 130  format ('   write (1,*) nBodies                                 '/
     +        '   do j=1,nBodies                                      '/
     +        '     write (1,*) nBodyFaces(j)                         '/
     +        '     do i=1,nBodyFaces(j)                              '/
     +        '       write (1,*) node1(i),node2(i)                   '/
     +        '     end do                                            '/
     +        '   end do                                              '/
     +        'c..Two nodes for face on body numbered as above        '/
     +        '   write (1,*) nBoundaryFaces                          '/
     +        '   do i=1,nBoundaryFaces                               '/
     +        '     write (1,*) node1(i),node2(i)                     '/
     +        '   end do                                              '/
     +        'c..Two nodes for face on boundary numbered as above    '/
     +        '                                                       ')
  111 FORMAT('.',$)
  135 FORMAT(5X,'Writing the coordinates ... ')
  145 FORMAT(5X,'Writing the variables ... ')
  150 FORMAT(5X,'I am opening the Dplot file ... ',/)
  175 FORMAT(' done !',/)
  180 FORMAT(5X,'Su2 file WRITTEN',/)
  201 FORMAT('NDIME= ',I1)
  202 FORMAT('NELEM= ',I8)
  203 FORMAT('NPOIN= ',I8)
  204 FORMAT('NMARK= ',I8)
  205 FORMAT('MARKER_TAG= ',I2)
  206 FORMAT('MARKER_ELEMS= ',I8)
      END

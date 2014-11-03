C -------------------------------------------------------------
C CHECKDESIGN cheks design.dll or design.m
C Developed by A.Rossi, C.Planas and G.Fiorentini
C
C State-space format:   y(t) = c(t)z(t) + H(t)x(t)   + G(t)u(t)
C                       x(t) = a(t)     + F(t)x(t-1) + R(t)u(t)
C
C y(t) (ny x 1)          ny  = # of endogenous series
C z(t) (nz x 1)          nz  = # of exogenous series
C x(t) (nx x 1)          nx  = # of continous states
C u(t) (nu x 1)          nu  = # of shocks
C c(t) (ny x nz x ns1)   ns1 = # of states for c(t)
C H(t) (ny x nx x ns2)   ns2 = # of states for H(t)
C G(t) (ny x nu x ns3)   ns3 = # of states for G(t)
C a(t) (nx x ns4)        ns4 = # of states for a(t)
C F(t) (nx x nx x ns5)   ns5 = # of states for F(t)
C R(t) (nx x nu x ns6)   ns6 = # of states for R(t)
C
C Copyright (C) 2010-2014 European Commission
C
C This file is part of Program DMM
C
C DMM is free software developed at the Joint Research Centre of the
C European Commission: you can redistribute it and/or modify it under
C the terms of the GNU General Public License as published by
C the Free Software Foundation, either version 3 of the License, or
C (at your option) any later version.
C
C DMM is distributed in the hope that it will be useful,
C but WITHOUT ANY WARRANTY; without even the implied warranty of
C MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C GNU General Public License for more details.
C
C You should have received a copy of the GNU General Public License
C along with DMM.  If not, see <http://www.gnu.org/licenses/>.
C -------------------------------------------------------------
	SUBROUTINE CHECKDESIGN(ny,nz,nx,nu,ns,nt,d,theta,PATH,NMLNAME)
#if defined(__CYGWIN32__) || defined(_WIN32)
#ifdef __INTEL_COMPILER
      USE dfwin
#endif
#else
#ifdef __GFORTRAN__
      USE gfortran
#endif
      USE ISO_C_BINDING
      USE ISO_C_UTILITIES
      USE DLFCN
#endif

C INPUT
	INTEGER ny,nz,nx,nu,ns(6),nt,d(2)
	DOUBLE PRECISION theta(nt)
	CHARACTER*200 NMLNAME,PATH,FILEOUT

C LOCALS
	INTEGER I,J,maxnz,IFAIL,ESTABLE
	DOUBLE PRECISION WORK(4*nx),WR(nx),WI(nx),VR(1),
	1 VI(1),W(nx)  !WRY(ny),WORK1(64*ny)
	DOUBLE PRECISION,ALLOCATABLE::c(:,:,:),H(:,:,:),
	1 G(:,:,:),a(:,:),F(:,:,:),R(:,:,:)  !,HRG(:,:),HRGRH(:,:)
	CHARACTER*3 CJ
C EXTERNAL SUBROUTINES
      EXTERNAL DGEEV
#ifdef __GFORTRAN__
      INTEGER tmp
      CHARACTER*13 fmt
#endif

	ALLOCATE(c(ny,max(nz,1),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6))) !,HRG(ny,nu),HRGRH(ny,ny))

#if defined(ORIGDLL) || defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
	  CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
#else
#endif

	maxnz = max(1,nz)
	FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//'.CHK'
 	OPEN(11,FILE = FILEOUT, ACCESS='SEQUENTIAL')

C write ny,nx,etc
      WRITE(11,1000) nt,ny,nx,nu,nz
1000  FORMAT(' nt = ',I3,'; ny = ',I3,'; nx = ',I3,'; nu = ',I3,
     #       '; nz = ',I3,';')

C write theta
	 WRITE(11,*) 'theta(1:nt) = ['
	 WRITE(11,*) ' '
#ifdef __GFORTRAN__
      WRITE(fmt, '(a,i4,a)') '(', nt, '(F20.10))'
      WRITE(11,fmt) theta(1:nt)
#else
	 WRITE(11,'(<nt>(F20.10))') theta(1:nt)
#endif
       WRITE(11,*) ']'

C write c(ny,max(1,nz),ns(1))
	DO J = 1,ns(1)
	 WRITE(11,*) ' '
	 IF (J.LE.9) THEN
	  WRITE(CJ,'(I1)') J
	 ELSEIF ((J.GE.10).AND.(J.LE.99)) THEN
	  WRITE(CJ,'(I2)') J
	 ELSEIF ((J.GE.100).AND.(J.LE.999)) THEN
	  WRITE(CJ,'(I3)') J
	 ENDIF
	 WRITE(11,*) 'c(1:ny,1:nz,'//TRIM(CJ)// ') = ['
	 WRITE(11,*) ' '
#ifdef __GFORTRAN__
      WRITE(fmt, '(a,i4,a)') '(', maxnz, '(F20.10))'
      WRITE(11,fmt) (c(I,1:maxnz,J),I=1,ny)
#else
	 WRITE(11,'(<maxnz>(F20.10))') (c(I,1:maxnz,J),I=1,ny)
#endif
       WRITE(11,*) ']'
	ENDDO

C write H(ny,nx,ns(2))
	DO J = 1,ns(2)
	 WRITE(11,*) ' '
	 IF (J.LE.9) THEN
	  WRITE(CJ,'(I1)') J
	 ELSEIF ((J.GE.10).AND.(J.LE.99)) THEN
	  WRITE(CJ,'(I2)') J
	 ELSEIF ((J.GE.100).AND.(J.LE.999)) THEN
	  WRITE(CJ,'(I3)') J
	 ENDIF
	 WRITE(11,*) 'H(1:ny,1:nx,'//TRIM(CJ)// ') = ['
	 WRITE(11,*) ' '
#ifdef __GFORTRAN__
      WRITE(fmt, '(a,i4,a)') '(', nx, '(F20.10))'
      WRITE(11,fmt) (H(I,1:nx,J),I=1,ny)
#else
	 WRITE(11,'(<nx>(F20.10))') (H(I,1:nx,J),I=1,ny)
#endif
       WRITE(11,*) ']'
	ENDDO

C write G(ny,nu,ns(3))
	DO J = 1,ns(3)
	 WRITE(11,*) ' '
	 IF (J.LE.9) THEN
	  WRITE(CJ,'(I1)') J
	 ELSEIF ((J.GE.10).AND.(J.LE.99)) THEN
	  WRITE(CJ,'(I2)') J
	 ELSEIF ((J.GE.100).AND.(J.LE.999)) THEN
	  WRITE(CJ,'(I3)') J
	 ENDIF
	 WRITE(11,*) 'G(1:ny,1:nu,'//TRIM(CJ)// ') = ['
	 WRITE(11,*) ' '
#ifdef __GFORTRAN__
      WRITE(fmt, '(a,i4,a)') '(', nu, '(F20.10))'
      WRITE(11,fmt) (G(I,1:nu,J),I=1,ny)
#else
	 WRITE(11,'(<nu>(F20.10))') (G(I,1:nu,J),I=1,ny)
#endif
       WRITE(11,*) ']'
	ENDDO

C write a(nx,ns(4))
	DO J = 1,ns(4)
	 WRITE(11,*) ' '
	 IF (J.LE.9) THEN
	  WRITE(CJ,'(I1)') J
	 ELSEIF ((J.GE.10).AND.(J.LE.99)) THEN
	  WRITE(CJ,'(I2)') J
	 ELSEIF ((J.GE.100).AND.(J.LE.999)) THEN
	  WRITE(CJ,'(I3)') J
	 ENDIF
	 WRITE(11,*) 'a(1:nx,'//TRIM(CJ)// ') = ['
	 WRITE(11,*) ' '
#ifdef __GFORTRAN__
      WRITE(fmt, '(a,i4,a)') '(', 1, '(F20.10))'
      WRITE(11,fmt) (a(I,J),I=1,nx)
#else
	 WRITE(11,'(<1>(F20.10))') (a(I,J),I=1,nx)
#endif
       WRITE(11,*) ']'
	ENDDO

C write F(nx,nx,ns(5))
	DO J = 1,ns(5)
	 WRITE(11,*) ' '
	 IF (J.LE.9) THEN
	  WRITE(CJ,'(I1)') J
	 ELSEIF ((J.GE.10).AND.(J.LE.99)) THEN
	  WRITE(CJ,'(I2)') J
	 ELSEIF ((J.GE.100).AND.(J.LE.999)) THEN
	  WRITE(CJ,'(I3)') J
	 ENDIF
	 WRITE(11,*) 'F(1:nx,1:nx,'//TRIM(CJ)// ') = ['
	 WRITE(11,*) ' '
#ifdef __GFORTRAN__
      WRITE(fmt, '(a,i4,a)') '(', nx, '(F20.10))'
      WRITE(11,fmt) (F(I,1:nx,J),I=1,nx)
#else
	 WRITE(11,'(<nx>(F20.10))') (F(I,1:nx,J),I=1,nx)
#endif
       WRITE(11,*) ']'
      ENDDO

C write R(nx,nu,ns(6))
	DO J = 1,ns(6)
	 WRITE(11,*) ' '
	 IF (J.LE.9) THEN
	  WRITE(CJ,'(I1)') J
	 ELSEIF ((J.GE.10).AND.(J.LE.99)) THEN
	  WRITE(CJ,'(I2)') J
	 ELSEIF ((J.GE.100).AND.(J.LE.999)) THEN
	  WRITE(CJ,'(I3)') J
	 ENDIF
	 WRITE(11,*) 'R(1:nx,1:nu,'//TRIM(CJ)// ') = ['
	 WRITE(11,*) ' '
#ifdef __GFORTRAN__
      WRITE(fmt, '(a,i4,a)') '(', nu, '(F20.10))'
      WRITE(11,fmt) (R(I,1:nu,J),I=1,nx)
#else
	 WRITE(11,'(<nu>(F20.10))') (R(I,1:nu,J),I=1,nx)
#endif
       WRITE(11,*) ']'
      ENDDO

C Check unstable eigenvalues of F
      DO J = 1,ns(5)
       IF (d(2).GT.0) THEN
	  IFAIL=-1
C	  CALL F02EBF('N',d(2),F(1:d(2),1:d(2),J),d(2),
C	1              WR(1:d(2)),WI(1:d(2)),VR,1,VI,1,WORK,4*nx,IFAIL)
        CALL DGEEV('N','N',d(2),F(1:d(2),1:d(2),J),d(2),
     1             WR(1:d(2)),WI(1:d(2)),VR,1,VI,1,WORK,4*nx,IFAIL)

	  ESTABLE = 0
	  DO I = 1,d(2)
	   W(I) = WR(I)**2+WI(I)**2
#ifdef __GFORTRAN__
       tmp = W(I).GE.1.D0
       ESTABLE = ESTABLE + ABS(tmp)
#else
	   ESTABLE = ESTABLE + ABS(W(I).GE.1.D0)
#endif
	  ENDDO
        IF (ESTABLE.NE.d(2)) THEN
	   WRITE(11,*) ' '
	   WRITE(11,*) 'WARNING: the number of unstable eigenvalues for '
	   WRITE(11,*) 'F(1:nx,1:nx,'//TRIM(CJ)// 'is not equal to d(2) '
     	   WRITE(11,*) 'or non-stationary states are not placed in the'
         WRITE(11,*) 'first d(2) positions.'
	  ENDIF
       ENDIF

C Check stable eigenvalues of F
       IF (nx-d(2).GT.0) THEN
	  IFAIL=-1
c	  CALL F02EBF('N',nx-d(2),F(d(2)+1:nx,d(2)+1:nx,J),
c	1              nx-d(2),WR,WI,VR,1,VI,1,WORK,4*nx,IFAIL)
        CALL DGEEV('N','N',nx-d(2),F(d(2)+1:nx,d(2)+1:nx,J),
     #              nx-d(2),WR,WI,VR,1,VI,1,WORK,4*nx,IFAIL)

	  ESTABLE = 0
	  DO I = 1,nx-d(2)
         W(I) = WR(I)**2+WI(I)**2
#ifdef __GFORTRAN__
         tmp = W(I).LT.1.D0
         ESTABLE = ESTABLE + ABS(tmp)
#else
	   ESTABLE = ESTABLE + ABS(W(I).LT.1.D0)
#endif
	  ENDDO
        IF (ESTABLE.NE.(nx-d(2))) THEN
	   WRITE(11,*) ' '
	   WRITE(11,*) 'WARNING: the number of stable eigenvalues for '
	   WRITE(11,*) 'F(1:nx,1:nx,'//TRIM(CJ)//'is not equal to nx-d(2)'
     	   WRITE(11,*) 'or non-stationary states are not placed in the '
         WRITE(11,*) 'first d(2) positions.'
	  ENDIF
	 ENDIF
      ENDDO

	CLOSE(11)
      DEALLOCATE(c,H,G,a,F,R)

	RETURN
      PAUSE
      END

C Check rank{(HR+G)*(HR+G)'} this check is wrong!!!
c	DO J   = 1,ns(2) !H
c	 DO JJ  = 1,ns(3) !G
c	  DO JJJ = 1,ns(6) !R

c	  DO I =1,ny
c	   DO K =1,nu
c	    HRG(I,K) = SUM(H(I,1:nx,J)*R(1:nx,K,JJJ))+G(I,K,JJ)
c	   ENDDO
c	  ENDDO

c	  DO I =1,ny
c	   DO K =1,ny
c	    HRGRH(I,K) = SUM(HRG(I,1:nu)*HRG(K,1:nu))
c	   ENDDO
c	  ENDDO

c	  IFAIL = -1
c	  CALL F02FAF('N','U',ny,HRGRH,ny,WRY(1:ny),WORK1,64*ny,IFAIL)
c	  SRANK = 0
c	  DO 10 I=1,ny
c10      IF (WRY(I).GT.1.D-12) SRANK=SRANK+1

c        IF (SRANK.LT.ny) THEN
c	   WRITE(11,*) ' '
c	   WRITE(11,*) 'WARNING: the rank of the system computed looking '
c	   WRITE(11,*) 'at rank{(HR+G)*transpose(HR+G)} is less than ny  '
c	  ENDIF

c	  ENDDO
c	 ENDDO
c	ENDDO


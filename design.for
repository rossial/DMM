C --------------------------------------------------------------------
C This file crates matlabdll.dll whih is ment to read the .m file
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
C ---------------------------------------------------------------------
#include "fintrf.h"
	SUBROUTINE DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:'design_' :: DESIGN

C INPUT
	INTEGER ny,nz,nx,nu,ns(6),nt
	DOUBLE PRECISION theta(nt)
C OUTPUT
	DOUBLE PRECISION c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
     &                 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),
     &                 R(nx,nu,ns(6))
C
C COMMON
C
      INTEGER ep
      COMMON /X/ ep
      DATA ep/0/
      CHARACTER*200 mfile,pathmfile
      COMMON /M/ mfile,pathmfile
      CHARACTER*1024 buffer,buffer1
      COMMON /ERRBUFFER/ buffer
C
C LOCALS TO INTERFACE MatLab
C
      INTEGER engOpen, engGetVariable, mxCreateDoubleMatrix
      INTEGER mxGetPr
      INTEGER engPutVariable, engEvalString, engClose
      INTEGER i, temp, status
	DOUBLE PRECISION nsd(6)
	INTEGER ny_ptr,nz_ptr,nx_ptr,nu_ptr,ns_ptr,nt_ptr,theta_ptr
	INTEGER C_ptr, H_ptr, G_ptr, A_ptr, F_ptr, R_ptr
      CHARACTER*1024 matlaberror

C Try to open MatLab (just the first time)
      IF (ep .eq.0 ) THEN
        ep = engOpen('/Applications/MATLAB_R2014b.app/bin/matlab ')
        IF (ep .eq. 0) THEN ! Can''t start Matlab engine
#ifdef __GFORTRAN__
		   WRITE(*,*) ' '
		   WRITE(*,*) ' Can''t start MATLAB engine'
		   WRITE(*,*) ' Program aborting'
#else
		   TYPE *, ' '
		   TYPE *, ' Can''t start MATLAB engine'
		   TYPE *, ' Program aborting'
		   PAUSE
#endif
		   STOP
        ENDIF
        IF (engEvalString(ep,'cd ' // pathmfile).ne. 0) then ! Can''t find or open the MatLab funtion
#ifdef __GFORTRAN__
		   WRITE(*,*) ' '
		   WRITE(*,*) ' Can''t find or open the MatLab function'
		   WRITE(*,*) ' Program aborting'
#else
		   TYPE *, ' '
		   TYPE *, ' Can''t find or open the MatLab function'
		   TYPE *, ' Program aborting'
		   PAUSE
#endif
		   STOP
        ENDIF
      ENDIF

C
C Check MatLab INPUT
C
	ny_ptr = mxCreateDoubleScalar(ny*1.0d0)
#ifdef __GFORTRAN__
      status = engPutVariable(ep, 'ny', ny_ptr)
#else
      status = engPutVariable(ep, 'ny'C, ny_ptr)
#endif
      IF (status .ne. 0) THEN ! Can''t read ny in the Matlab file
#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' Can''t read ny in the MATLAB file'
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' Can''t read ny in the MATLAB file'
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF
	nz_ptr = mxCreateDoubleScalar(nz*1.0d0)
#ifdef __GFORTRAN__
      status = engPutVariable(ep, 'nz', nz_ptr)
#else
      status = engPutVariable(ep, 'nz'C, nz_ptr)
#endif
      IF (status .ne. 0) THEN ! Can''t read nz in the Matlab file
#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' Can''t read nz in the MATLAB file'
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' Can''t read nz in the MATLAB file'
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF

	nx_ptr = mxCreateDoubleScalar(nx*1.0d0)
#ifdef __GFORTRAN__
      status = engPutVariable(ep, 'nx', nx_ptr)
#else
      status = engPutVariable(ep, 'nx'C, nx_ptr)
#endif
      IF (status .ne. 0) THEN ! Can''t read nx in the Matlab file
#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' Can''t read nx in the MATLAB file'
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' Can''t read nx in the MATLAB file'
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF

	nu_ptr = mxCreateDoubleScalar(nu*1.0d0)
#ifdef __GFORTRAN__
      status = engPutVariable(ep, 'nu', nu_ptr)
#else
      status = engPutVariable(ep, 'nu'C, nu_ptr)
#endif
      IF (status .ne. 0) THEN ! Can''t read nu in the MatLab file
#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' Can''t read nu in the MATLAB file'
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' Can''t read nu in the MATLAB file'
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF

      ns_ptr = mxCreateDoubleMatrix(1, 6, 0)
	DO i=1,6
	  nsd(i)=ns(i)*1.d0
	ENDDO
      CALL mxCopyReal8ToPtr(nsd, mxGetPr(ns_ptr), 6)
#ifdef __GFORTRAN__
      status = engPutVariable(ep, 'ns', ns_ptr)
#else
      status = engPutVariable(ep, 'ns'C, ns_ptr)
#endif
      IF (status .ne. 0) THEN ! Can''t read ns in the Matlab file
#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' Can''t read ns in the MATLAB file'
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' Can''t read ns in the MATLAB file'
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF

      theta_ptr = mxCreateDoubleMatrix(1, nt, 0)
      CALL mxCopyReal8ToPtr(theta, mxGetPr(theta_ptr), nt)
#ifdef __GFORTRAN__
      status = engPutVariable(ep, 'theta', theta_ptr)
#else
      status = engPutVariable(ep, 'theta'C, theta_ptr)
#endif
      IF (status .ne. 0) THEN ! Can''t read theta in the Matlab file
#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' Can''t read nt in the MATLAB file'
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' Can''t read nt in the MATLAB file'
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF

C
C     Evaluate the MatLab DESIGN Funtion with the input data from FORTRAN
C
      buffer = ''
      status = engOutputBuffer(ep, buffer1)
      IF (engEvalString(ep, 'clear success;'//
     &   '[C,H,G,A,F,R]='//TRIM(mfile)//'( ny,nz,nx,'//
     &   'nu,ns,theta);'//'success=1;') .ne. 0) then ! engEvalString failed
#if defined(MEX)
		 CALL GETERRSTR(matlaberror)
#else
#endif

#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' the MATLAB funtion can not be executed:'
		 WRITE(*,*) trim(matlaberror)
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' the MATLAB funtion can not be executed:'
		 TYPE *, trim(matlaberror)
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF
#ifdef __GFORTRAN__
      C_ptr = engGetVariable(ep, 'success')
#else
      C_ptr = engGetVariable(ep, 'success'C)
#endif
      IF (C_ptr .eq. 0) then ! engEvalString failed
          buffer=buffer1
#if defined(MEX)
		 CALL GETERRSTR(matlaberror)
#else
#endif

#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' the MATLAB funtion can not be executed:'
		 WRITE(*,*) trim(matlaberror)
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' the MATLAB funtion can not be executed:'
		 TYPE *, trim(matlaberror)
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF

C
C Get the MatLab DESIGN created matrics back to FORTRAN
C
#ifdef __GFORTRAN__
      C_ptr = engGetVariable(ep, 'C')
#else
      C_ptr = engGetVariable(ep, 'C'C)
#endif
      IF(C_ptr.NE.0) THEN
       CALL mxCopyPtrToReal8(mxGetPr(C_ptr), c, ny*max(1,nz)*ns(1))
      ELSE
#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' C could not be assigned during the call'
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' C could not be assigned during the call '
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF

#ifdef __GFORTRAN__
      H_ptr = engGetVariable(ep, 'H')
#else
      H_ptr = engGetVariable(ep, 'H'C)
#endif
      IF(H_ptr.NE.0) THEN
       CALL mxCopyPtrToReal8(mxGetPr(H_ptr), H, ny*nx*ns(2))
      ELSE
#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' H could not be assigned during the call'
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' H could not be assigned during the call '
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF

#ifdef __GFORTRAN__
      G_ptr = engGetVariable(ep, 'G')
#else
      G_ptr = engGetVariable(ep, 'G'C)
#endif
      IF(G_ptr.NE.0) THEN
       CALL mxCopyPtrToReal8(mxGetPr(G_ptr), G, ny*nu*ns(3))
      ELSE
#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' G could not be assigned during the call'
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' G could not be assigned during the call '
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF

#ifdef __GFORTRAN__
      A_ptr = engGetVariable(ep, 'A')
#else
      A_ptr = engGetVariable(ep, 'A'C)
#endif
      IF(A_ptr.NE.0) THEN
       CALL mxCopyPtrToReal8(mxGetPr(A_ptr), a, nx*ns(4))
      ELSE
#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' A could not be assigned during the call'
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' A could not be assigned during the call '
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF

#ifdef __GFORTRAN__
      F_ptr = engGetVariable(ep, 'F')
#else
      F_ptr = engGetVariable(ep, 'F'C)
#endif
      IF(F_ptr.NE.0) THEN
       CALL mxCopyPtrToReal8(mxGetPr(F_ptr), F, nx*nx*ns(5))
      ELSE
#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' F could not be assigned during the call'
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' F could not be assigned during the call '
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF

#ifdef __GFORTRAN__
      r_ptr = enggetvariable(ep, 'R')
#else
      r_ptr = enggetvariable(ep, 'R'c)
#endif
      IF(R_ptr.NE.0) THEN
       CALL mxCopyPtrToReal8(mxGetPr(R_ptr), R, nx*nu*ns(6))
      ELSE
#ifdef __GFORTRAN__
		 WRITE(*,*) ' '
		 WRITE(*,*) ' R could not be assigned during the call'
		 WRITE(*,*) ' Program aborting'
#else
		 TYPE *, ' '
		 TYPE *, ' R could not be assigned during the call '
		 TYPE *, ' Program aborting'
		 PAUSE
#endif
		 STOP
      ENDIF
C
C Free dynamic memory allocated by MXCREATE function
C
      CALL mxDestroyArray(nz_ptr)
      CALL mxDestroyArray(nx_ptr)
      CALL mxDestroyArray(nu_ptr)
      CALL mxDestroyArray(ns_ptr)
      CALL mxDestroyArray(theta_ptr)
      CALL mxDestroyArray(C_ptr)
      CALL mxDestroyArray(H_ptr)
      CALL mxDestroyArray(G_ptr)
      CALL mxDestroyArray(A_ptr)
      CALL mxDestroyArray(F_ptr)
      CALL mxDestroyArray(R_ptr)
      CALL mxDestroyArray(ny_ptr)

      RETURN
      END

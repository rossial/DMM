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


C Try to open MatLab (just the first time)
      IF (ep .eq.0 ) THEN
        ep = engOpen('matlab ')
        IF (ep .eq. 0) THEN
          ny = 0 ! ' Can''t start MatLab engine'
          RETURN
        ENDIF
        IF (engEvalString(ep,'cd ' // pathmfile).ne. 0) then
          ny = -7 ! ' Can''t find or open the MatLab funtion'
          RETURN
        ENDIF
      ENDIF

C
C Check MatLab INPUT
C
	ny_ptr = mxCreateDoubleScalar(ny*1.0d0)
      status = engPutVariable(ep, 'ny'C, ny_ptr)
      IF (status .ne. 0) THEN
         ny = -1 ! ' Can''t read ny in the MatLab file'
         RETURN
      ENDIF

	nz_ptr = mxCreateDoubleScalar(nz*1.0d0)
      status = engPutVariable(ep, 'nz'C, nz_ptr)
      IF (status .ne. 0) THEN
         ny = -2 ! ' Can''t read nz in the MatLab file'
         RETURN
      ENDIF

	nx_ptr = mxCreateDoubleScalar(nx*1.0d0)
      status = engPutVariable(ep, 'nx'C, nx_ptr)
      IF (status .ne. 0) THEN
         ny = -3 ! ' Can''t read nx in the MatLab file'
         RETURN
      ENDIF

	nu_ptr = mxCreateDoubleScalar(nu*1.0d0)
      status = engPutVariable(ep, 'nu'C, nu_ptr)
      IF (status .ne. 0) THEN
         ny = -4 ! ' Can''t read nu in the MatLab file'
         RETURN
      ENDIF

      ns_ptr = mxCreateDoubleMatrix(1, 6, 0)
	DO i=1,6
	  nsd(i)=ns(i)*1.d0
	ENDDO
      CALL mxCopyReal8ToPtr(nsd, mxGetPr(ns_ptr), 6)
      status = engPutVariable(ep, 'ns'C, ns_ptr)
      IF (status .ne. 0) THEN
         ny = -5 ! ' Can''t read ns in the MatLab file'
         RETURN
      ENDIF

      theta_ptr = mxCreateDoubleMatrix(1, nt, 0)
      CALL mxCopyReal8ToPtr(theta, mxGetPr(theta_ptr), nt)
      status = engPutVariable(ep, 'theta'C, theta_ptr)
      IF (status .ne. 0) THEN
         ny = -6 ! ' Can''t read theta in the MatLab file'
         RETURN
      ENDIF

C
C     Evaluate the MatLab DESIGN Funtion with the input data from FORTRAN
C
      buffer = ''
      status = engOutputBuffer(ep, buffer1)
      IF (engEvalString(ep, 'clear success;'//
     &   '[C,H,G,A,F,R]='//TRIM(mfile)//'( ny,nz,nx,'//
     &   'nu,ns,theta);'//'success=1;') .ne. 0) then
         ny = -8   ! engEvalString failed
         RETURN
      ENDIF
      C_ptr = engGetVariable(ep, 'success'C)
      IF (C_ptr .eq. 0) then
          buffer=buffer1
          ny=-8   ! engEvalString failed
          RETURN
      ENDIF

C
C Get the MatLab DESIGN created matrics back to FORTRAN
C
      C_ptr = engGetVariable(ep, 'C'C)
      IF(C_ptr.NE.0) THEN
       CALL mxCopyPtrToReal8(mxGetPr(C_ptr), c, ny*max(1,nz)*ns(1))
      ELSE
       ny = -101
      ENDIF

      H_ptr = engGetVariable(ep, 'H'C)
      IF(H_ptr.NE.0) THEN
       CALL mxCopyPtrToReal8(mxGetPr(H_ptr), H, ny*nx*ns(2))
      ELSE
       ny = -102
      ENDIF

      G_ptr = engGetVariable(ep, 'G'C)
      IF(G_ptr.NE.0) THEN
       CALL mxCopyPtrToReal8(mxGetPr(G_ptr), G, ny*nu*ns(3))
      ELSE
       ny = -103
      ENDIF

      A_ptr = engGetVariable(ep, 'A'C)
      IF(A_ptr.NE.0) THEN
       CALL mxCopyPtrToReal8(mxGetPr(A_ptr), a, nx*ns(4))
      ELSE
       ny = -104
      ENDIF

      F_ptr = engGetVariable(ep, 'F'C)
      IF(F_ptr.NE.0) THEN
       CALL mxCopyPtrToReal8(mxGetPr(F_ptr), F, nx*nx*ns(5))
      ELSE
       ny = -105
      ENDIF

      R_ptr = engGetVariable(ep, 'R'C)
      IF(R_ptr.NE.0) THEN
       CALL mxCopyPtrToReal8(mxGetPr(R_ptr), R, nx*nu*ns(6))
      ELSE
       ny = -106
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

C -----------------------------------------
C To make dynamic the name of the .m file
C -----------------------------------------
	SUBROUTINE SETFILEM(string1,string2)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:'setfilem_' :: SETFILEM
	CHARACTER*200 string1,string2,mfile,pathmfile
      COMMON /M/ mfile, pathmfile
      mfile     = string1
      pathmfile = string2
	RETURN
      END

C --------------------
C To get MatLab errors
C --------------------
      SUBROUTINE GETERRSTR(matlaberror)
!DEC$ ATTRIBUTES DLLEXPORT, ALIAS:'geterrstr_' :: GETERRSTR
      CHARACTER*1024 matlaberror
      CHARACTER*1024 buffer
      COMMON /ERRBUFFER/ buffer
	 matlaberror = buffer
      RETURN
      END

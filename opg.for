C ------------------------------------------------------------
C OPG coumputes the Var-Covar matrix of Parameters inverting
C the Observed Information matrix calcuted by outer product
C of gradient estimator - e.g. Davidson MacKinnon pp 265-66
C Developed by A.Rossi, C.Planas and G.Fiorentini
C
C OUTPUT: SE = Standard deviations;
C         IFAIL =  0 Hessian
C         IFAIL =  1 OPG
C         IFAIL = -1 Failure
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
C ------------------------------------------------------------
	SUBROUTINE OPG(nobs,d,ny,nz,nx,nu,nt,ns,yk,IYK,S,
	1 theta,thetaprior,HESS,SE,XS,AKMSE,INN,IFAIL)
#if defined(__CYGWIN32__) || defined(_WIN32)
#ifdef __INTEL_COMPILER
      USE dfwin
#endif
#else
      USE ISO_C_BINDING
      USE ISO_C_UTILITIES
      USE DLFCN
#endif

! Input
	INTEGER nobs,d(2),ny,nz,nx,nu,nt,ns(6),IYK(nobs,ny+1),
	1 S(nobs,6)
      DOUBLE PRECISION yk(nobs,ny+nz),theta(nt),thetaprior(nt,4),
     1 HESS(nt*(nt+1)/2)
! Output
      INTEGER IFAIL,IFAILSY
      DOUBLE PRECISION SE(nt),XS(nobs,nx),AKMSE(nobs,nx),INN(nobs,ny)
! Locals
      INTEGER I,J,K,IFREE(nt),NFREE
      DOUBLE PRECISION DRI,RMAX
      DOUBLE PRECISION, ALLOCATABLE:: THETAV(:),DLL(:),DLLM(:),XSM(:,:),
     1 XT(:,:),PT(:,:,:),Xdd(:,:),Pdd(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE:: c(:,:,:),H(:,:,:),
	1 G(:,:,:),a(:,:),F(:,:,:),R(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE:: GRAD(:,:),P(:),LTR(:),
	1 W(:),OP(:,:),MAT(:,:),VC(:,:),pro(:,:)

	SE(:) = 0.D0
      DRI   = 1.D-3
	NFREE = 0
      DO I = 1,nt
         IF ((theta(I).GT.thetaprior(I,3)).AND.
     $        (theta(I).LT.thetaprior(I,4))) THEN
            NFREE = NFREE + 1
            IFREE(NFREE) = I
         END IF
      END DO

C Using Hessian from E04UCF
	ALLOCATE (LTR(NFREE*(NFREE+1)/2),W(NFREE))
      DO 25 I=1,NFREE
	DO 25 J=1,I
      K = IFREE(I)*(IFREE(I)+1)/2-IFREE(I)+IFREE(J)
25    LTR(I*(I+1)/2-I+J) = HESS(K)
      IFAIL = 0
      CALL SYMINV(LTR,NFREE,LTR,W,J,IFAIL,RMAX)

	ALLOCATE(GRAD(nobs,NFREE),P(NFREE),OP(NFREE,NFREE),MAT(nobs*nx,NFREE),
     1 VC(NFREE,NFREE),pro(nobs*nx,NFREE))
      ALLOCATE(c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6)))
      ALLOCATE(THETAV(nt),DLL(nobs),DLLM(nobs),XSM(nobs,nx),
     1 XT(0:nobs,nx),PT(0:nobs,nx,nx),Xdd(max(d(1),1),nx),
     1 Pdd(max(d(1),1),nx,nx))

#if defined(ORIGDLL) || defined(MEX)
	  CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
#else
#endif
      CALL IKF(d,ny,nz,nx,nu,ns,S(1:max(d(1),1),1:6),
	1         yk(1:max(d(1),1),1:ny+nz),IYK(1:max(d(1),1),1:ny+1),
     2         c,H,G,a,F,R,Xdd,Pdd,DLL(1:max(d(1),1)))
	XT(d(1),1:nx)      = Xdd(max(d(1),1),1:nx)
	PT(d(1),1:nx,1:nx) = Pdd(max(d(1),1),1:nx,1:nx)
	CALL KF(nobs,d,ny,nz,nx,nu,ns,S,yk,IYK,c,H,G,a,F,R,XT,PT,DLL)
	CALL KS(nobs,d,ny,nz,nx,nu,ns,S,yk,IYK,c,H,G,a,F,R,XS,
	1        PT(1:nobs,1:nx,1:nx))
      CALL INNOV(nobs,d,ny,nz,nx,nu,ns,nt,S,yk,IYK,theta,pdll,INN)

	AKMSE(1:d(1),:) = 0.D0
	DO I = d(1)+1,nobs
      DO J = 1,nx
      AKMSE(I,J) = PT(I,J,J)
      END DO
      END DO

      THETAV(1:NFREE) = theta(IFREE(1:NFREE))
	P(1:NFREE) = THETAV(1:NFREE)*DRI
	DO I=1,NFREE
       IF(DABS(P(I)).LT.1.D-13) P(I)=1.D-13
       IF (((theta(IFREE(I))+P(I)).GT.THETAPRIOR(IFREE(I),4)).OR.
	1    ((theta(IFREE(I))+P(I)).LT.THETAPRIOR(IFREE(I),3))) THEN
        IF (theta(IFREE(I)).GT.0.D0) THEN
         P(I) = THETAPRIOR(IFREE(I),3) - theta(IFREE(I))
        ELSE
          P(I) = THETAPRIOR(IFREE(I),4) - theta(IFREE(I))
        ENDIF
       ENDIF
      END DO

	DO 1000 I=1,NFREE
	 THETAV(I) = THETAV(I) + P(I)
       theta(IFREE(I)) = THETAV(I)
#if defined(ORIGDLL) || defined(MEX)
	   CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
#else
#endif
	 CALL IKF(d,ny,nz,nx,nu,ns,S(1:max(d(1),1),1:6),
	1          yk(1:max(d(1),1),1:ny+nz),IYK(1:max(d(1),1),1:ny+1),
     2          c,H,G,a,F,R,Xdd,Pdd,DLLM(1:max(d(1),1)))
	 XT(d(1),1:nx)      = Xdd(max(d(1),1),1:nx)
	 PT(d(1),1:nx,1:nx) = Pdd(max(d(1),1),1:nx,1:nx)
	 CALL KF(nobs,d,ny,nz,nx,nu,ns,S,yk,IYK,c,H,G,a,F,R,XT,PT,DLLM)
	 CALL KS(nobs,d,ny,nz,nx,nu,ns,S,yk,IYK,c,H,G,a,F,R,XSM,
	1         PT(1:nobs,1:nx,1:nx))
	 THETAV(I) = THETAV(I) - P(I)
	 theta(IFREE(I)) = THETAV(I)
       GRAD(2:nobs,I) = (DLLM(2:nobs) - DLL(2:nobs))/P(I)

	 DO J = 1,nx ! MAT (nobs x nx) x nfree
	  MAT(1+(J-1)*nobs:J*nobs,I) = (XSM(:,J)-XS(:,J))/P(I)
       ENDDO
1000  CONTINUE

C Use OPG if Hessian is bad
	 IF (IFAIL.GT.0) THEN
        DO 300 I=1,NFREE
	  DO 300 J=1,I
300	  OP(I,J) = SUM(GRAD(2:nobs,I)*GRAD(2:nobs,J))

	  DO 150 I=1,NFREE
	  DO 150 J=1,I
150     LTR(i*(i+1)/2-i+j)=OP(I,J)
        IFAIL = 1
	  CALL SYMINV(LTR,NFREE,LTR,W,J,IFAILSY,RMAX)
        IF (IFAILSY.NE.0) THEN
         IFAIL = -1
         GOTO 1111
        ENDIF
       ENDIF

C ------------------------------------------------------
C Computes MAT(i,:)*VC*MAT(:,i) for each i=1,2,..,nobs
C ------------------------------------------------------
       DO 170 i=1,NFREE
	 VC(i,i)=LTR(i*(i+1)/2)
       SE(IFREE(i)) = DSQRT(LTR(i*(i+1)/2))
	 DO 170 j=1,i-1
       VC(i,j)=LTR(i*(i+1)/2-i+j)
170    VC(j,i)=VC(i,j)

	 DO 301 I = 1,nobs*nx ! pro = MAT         x  VC
       DO 301 J = 1,Nfree  !      (nobsxnx)xnt x  nt x nt
301    pro(i,j) = SUM(MAT(i,1:NFREE)*VC(1:NFREE,j))

C AKMSE: first nobs = var of nobs estimate of first state element
	 DO I = d(1)+1,nobs
	 DO J = 1,nx
       AKMSE(I,J) = AKMSE(I,J) + SUM(pro(nobs*(J-1)+I,1:NFREE)*
	1               MAT(nobs*(J-1)+I,1:NFREE))
	 ENDDO
       ENDDO

1111  DO I = 1,nobs
       DO J = 1,nx
        AKMSE(I,J) = DSQRT(AKMSE(I,J))
	 ENDDO
	ENDDO

      DEALLOCATE(OP,W,LTR,P,GRAD,MAT,VC,pro)
      DEALLOCATE(c,H,G,a,F,R,THETAV,DLL,DLLM,XSM,XT,PT,Xdd,Pdd)

	RETURN
      END

cELSE
c       CALL KIM2(nobs,d,ny,nz,nx,nu,ns,PRODUCT(ns),nv,np,INFOS,yk,
c     1           c,H,G,a,F,R,psi,1,XS,AKMSE,SSMOOTH,DLL)
c        DO I=1,NFREE
c	   SE(IFREE(I)) = DSQRT(LTR(I*(I+1)/2))
c        ENDDO
c      ENDIF

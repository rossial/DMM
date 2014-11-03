C --------------------------------------------------------------------
C SIMDATA simulates (S, x, yk) given (theta,psi)
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
C H(t) (ny x nx x ns2)   ns2 = # of states for S2(t)
C G(t) (ny x nu x ns3)   ns3 = # of states for S3(t)
C a(t) (nx x ns4)        ns4 = # of states for S4(t)
C F(t) (nx x nx x ns5)   ns5 = # of states for S5(t)
C R(t) (nx x nu x ns6)   ns6 = # of states for S6(t)
C
C Set latent variables: INFOS (9 x nv)
C by cols: S1,S2,...,Snv; with nv <=6
C by row: the 1st contains the # of matrices affected by Si
C         the 2nd-3rd etc point to c (1),H (2),G (3),a (4),F (5),R (6)
C         the 8-th row contains the # of states
C         the 9-th row spec the dynamics for Sj

C OUTPUT:
C
C  S     ~ p(S|psi)
C  STATE ~ p(x|S,theta)
C  yk    ~ p(y|S,x,theta)
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
C --------------------------------------------------------------------
	SUBROUTINE SIMDATA(nobs,d,ny,nz,nx,nu,ns,nstot,nt,nv,np,INFOS,
	1                   theta,psi,Z,STATE,yk)
#if defined(__CYGWIN32__) || defined(_WIN32)
#ifdef __INTEL_COMPILER
      USE dfwin
#endif
#else
      USE ISO_C_BINDING
      USE ISO_C_UTILITIES
      USE DLFCN
#endif

C INPUT
	INTEGER nobs,d(2),ny,nz,nx,nu,ns(6),nstot,nt,nv,np(3),INFOS(9,6)
	DOUBLE PRECISION theta(nt),psi(np(1))
C OUTPUT
	INTEGER Z(nobs)
	DOUBLE PRECISION STATE(nobs,nx),yk(nobs,ny+nz)

C LOCALS
	INTEGER ISEQ,I,J,K,it,IFAIL,IPIV(nx)
	INTEGER S(nobs,6),SEQ(nv)
	DOUBLE PRECISION P1(INFOS(8,1),INFOS(8,1)),
	1 P2(INFOS(8,2),INFOS(8,2)),P3(INFOS(8,3),INFOS(8,3)),
     2 P4(INFOS(8,4),INFOS(8,4)),P5(INFOS(8,5),INFOS(8,5)),
     3 P6(INFOS(8,6),INFOS(8,6)),PMAT(nstot,nstot),PE(nstot)

C EXTERNAL FUNCTIONS
      DOUBLE PRECISION genunf,gennor
C EXTERNAL SUBROUTINES
      EXTERNAL DGETRF,DGETRI,DESIGNZ,PPROD,ERGODIC,INT2SEQ,LYAP,SETGMN,
     1 GENMN

	DOUBLE PRECISION U,AUX
	DOUBLE PRECISION,ALLOCATABLE::R(:,:,:),c(:,:,:),H(:,:,:),
	1 G(:,:,:),a(:,:),F(:,:,:)
	DOUBLE PRECISION,ALLOCATABLE:: Xdd(:,:),Pdd(:,:,:),WORK(:),
	1 FP(:,:),WORK1(:),UP(:)


	ALLOCATE (R(nx,nu,ns(6)),c(ny,max(nz,1),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),
     3 Xdd(MAX(d(1),1),nx),Pdd(MAX(d(1),1),nx,nx),
     3 WORK((nx+2)*(nx+1)/2),FP(nx,nx),WORK1(64*nx),UP(nu) )

#if defined(ORIGDLL) || defined(MATLAB_MEX_FILE) || defined(OCTAVE_MEX_FILE)
	  CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
#else
#endif

C DRAW Z ~ Pr(S1 x ... x Snv|psi)
	S(:,:) = 1
      IF (nv.GT.0) THEN
	 CALL DESIGNZ(nv,np(1),psi,INFOS,P1,P2,P3,P4,P5,P6)
C PALL(i,j) = Pr[Z(t+1)=i|Z(t)=j], Z = S1 x S2 x ... x Snv
	 CALL PPROD(nv,nstot,INFOS,P1,P2,P3,P4,P5,P6,PMAT)
C ERGODIC solves PE: PE*(I-P') = 0
       CALL ERGODIC(nstot,PMAT,PE)
C	 U = G05CAF(U) ! Sampling from U(0,1)
       U = genunf(0.D0,1.D0)   ! Sampling from U(0,1)
	 ISEQ = 1
	 AUX  = PE(1)
	 DO 5 WHILE (AUX.LT.U)
	 ISEQ = ISEQ + 1
5	 AUX  = AUX  + PE(ISEQ)
	 Z(1) = ISEQ
	 CALL INT2SEQ(Z(1),nv,INFOS,SEQ,S(1,:))
	 DO it=2,nobs
C	  U = G05CAF(U) ! Sampling from U(0,1)
        U = genunf(0.D0,1.D0)   ! Sampling from U(0,1)
	  ISEQ = 1
	  AUX  = PMAT(1,Z(it-1))
	  DO 10 WHILE (AUX.LT.U)
	  ISEQ = ISEQ + 1
10	  AUX  = AUX  + PMAT(ISEQ,Z(it-1))
	  Z(it) = ISEQ
 	  CALL INT2SEQ(Z(it),nv,INFOS,SEQ,S(it,:))
	 ENDDO
	ELSE
	 S(:,:) = 1
	 Z(:)   = 1
	ENDIF

	yk(1:nobs,1:ny) = 0.D0
C DRAW x(1) ~ N[x(1|0),P(1|0)]
	STATE(1,1:nx) = 0.D0
	IF (d(2).LT.nx) THEN
	 Xdd(1,1:nx) = 0.D0
	 Pdd(1,1:nx,1:nx) = 0.D0
	 IF(SUM(ABS(a(d(2)+1:nx,S(1,4)))).NE.0.D0) THEN
	   FP(d(2)+1:nx,d(2)+1:nx) = -F(d(2)+1:nx,d(2)+1:nx,S(1,5))
	   DO 25 I = d(2)+1,nx
25	   FP(I,I) = 1.D0+FP(I,I)
         IFAIL = -1
C	   CALL F07ADF(nx-d(2),nx-d(2),FP(d(2)+1:nx,d(2)+1:nx),nx-d(2),
C	1               IPIV(d(2)+1:nx),IFAIL)
C	   CALL F07AJF(nx-d(2),FP(d(2)+1:nx,d(2)+1:nx),nx-d(2),
C	1               IPIV(d(2)+1:nx),WORK1,64*nx,IFAIL)
	   CALL DGETRF(nx-d(2),nx-d(2),FP(d(2)+1:nx,d(2)+1:nx),nx-d(2),
	1               IPIV(d(2)+1:nx),IFAIL)
	   CALL DGETRI(nx-d(2),FP(d(2)+1:nx,d(2)+1:nx),nx-d(2),
	1               IPIV(d(2)+1:nx),WORK1,64*nx,IFAIL)

	   DO 30 I = d(2)+1,nx
30	   Xdd(1,I) = SUM(FP(I,d(2)+1:nx)*a(d(2)+1:nx,S(1,4))) ! inv(I-F)*a
       ENDIF
	 CALL LYAP(nx-d(2),nu,1.D-3,F(d(2)+1:nx,d(2)+1:nx,S(1,5)),
	1           R(d(2)+1:nx,1:nu,S(1,6)),Pdd(1,d(2)+1:nx,d(2)+1:nx))
	 IFAIL = -1
C	 CALL G05EAF(Xdd(1,d(2)+1:nx),nx-d(2),Pdd(1,d(2)+1:nx,d(2)+1:nx),
C	1             nx-d(2),10.D-14,WORK,(nx+2)*(nx+1)/2,IFAIL)
C	 CALL G05EZF(STATE(1,d(2)+1:nx),nx-d(2),WORK,(nx+2)*(nx+1)/2,
C	1             IFAIL)
       CALL setgmn(Xdd(1,d(2)+1:nx),Pdd(1,d(2)+1:nx,d(2)+1:nx),nx-d(2),
     #             nx-d(2),WORK(1:(nx-d(2)+2)*(nx-d(2)+1)/2))
       CALL genmn(WORK(1:(nx-d(2)+2)*(nx-d(2)+1)/2),STATE(1,d(2)+1:nx),
     #            WORK1(1:nx-d(2)))

      ENDIF

C DRAW u(1)
	DO 35 J = 1,nu
C35    UP(J) = G05DDF(0.0D0,1.D0)
35    UP(J) = gennor(0.0D0,1.D0)


C COMPUTE y(1)
	DO 36 K = 1,ny
36	yk(1,K) = SUM(H(K,1:nx,S(1,2))*STATE(1,1:nx))
     #        + SUM(G(K,1:nu,S(1,3))*UP(1:nu))
     #        + SUM(c(K,1:nz,S(1,1))*yk(1,ny+1:ny+nz))

	DO 100 it = 2,nobs
C DRAW u ~ N(0,I)
	 DO 40 J = 1,nu
C40     UP(J) = G05DDF(0.0D0,1.D0)
40     UP(J) = gennor(0.0D0,1.D0)

C COMPUTE x
	 DO 50 K = 1,nx
50	 STATE(it,K) = a(K,S(it,4))
     #             + SUM(F(K,1:nx,S(it,5))*STATE(it-1,1:nx))
     #             + SUM(R(K,1:nu,S(it,6))*UP(1:nu))

C COMPUTE yk
	 DO 60 K = 1,ny
60	 yk(it,K) = SUM(H(K,1:nx,S(it,2))*STATE(it,1:nx))
     #          + SUM(G(K,1:nu,S(it,3))*UP(1:nu))
     #          + SUM(c(K,1:nz,S(it,1))*yk(it,ny+1:ny+nz))

100   CONTINUE

      DEALLOCATE (R,c,H,G,a,F,Xdd,Pdd,WORK,FP,WORK1,UP)

	RETURN
	END

C --------------------------------------------------------------------
C SIMSTATE2 (no missing values) IMPLEMENTS THE SIMULATION SMOOTHER
C in Durbin and Koopman (2002), "A simple and efficient simulation smoother
C for state space time series analysis". Biometrika, 89, 3, 603-15
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
C OUTPUT:
C
C   STATE ~ p(x|y,theta,Z) (nobs x nx)
C
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
	SUBROUTINE SIMSTATE2(nobs,d,ny,nz,nx,nu,ns,nt,yk,
	1                     theta,S,pdll,STATE)

	USE dfwin
	INTERFACE
	 SUBROUTINE DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
	 INTEGER ny,nz,nx,nu,ns(6),nt
	 DOUBLE PRECISION theta(nt)
	 DOUBLE PRECISION c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6))
	 END SUBROUTINE
	END INTERFACE
	POINTER (pdll,fittizia)  ! ASSOCIATE  pointer pdll alla DLL ad una varibile fittizia
	POINTER (pdesign,DESIGN) ! IMPORTANT associo il puntatore pdesign alla Interface definita

C INPUT
	INTEGER nobs,d(2),ny,nz,nx,nu,nt,ns(6),S(nobs,6)
	DOUBLE PRECISION yk(nobs,ny+nz),theta(nt)
C OUTPUT
	DOUBLE PRECISION STATE(nobs,nx)
C LOCALS
	INTEGER I,J,K,it,IFAIL,IPIV(nx)
	DOUBLE PRECISION,ALLOCATABLE:: ykP(:,:),XP(:,:),XS(:,:)
      DOUBLE PRECISION,ALLOCATABLE::R(:,:,:),c(:,:,:),H(:,:,:),
	1 G(:,:,:),a(:,:),F(:,:,:)
	DOUBLE PRECISION,ALLOCATABLE:: Xdd(:,:),Pdd(:,:,:),WORK(:),
	1 FP(:,:),WORK1(:),UP(:)
C EXTERNAL SUBROUTINES
      EXTERNAL DGETRF,DGETRI,SETGMN,GENMN,LYAP,KS2
C EXTERNAL FUNCTIONS
      DOUBLE PRECISION GENNOR

	ALLOCATE (R(nx,nu,ns(6)),c(ny,max(nz,1),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),
     2 ykP(nobs,ny+nz),XP(nobs,nx),XS(nobs,nx),
     3 Xdd(MAX(d(1),1),nx),Pdd(MAX(d(1),1),nx,nx),
     3 WORK((nx+2)*(nx+1)/2),FP(nx,nx),WORK1(64*nx),UP(nu) )

	pdesign = getprocaddress(pdll, "design_"C)
      CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)

	ykP(:,:) = 0.D0
C DRAW x(1)+ FROM N[x(1|0),P(1|0)]
	XP(1,1:nx) = 0.D0
	IF (d(2).LT.nx) THEN
	 Xdd(1,1:nx) = 0.D0
	 Pdd(1,1:nx,1:nx) = 0.D0
	 IF(SUM(ABS(a(d(2)+1:nx,S(1,4)))).NE.0.D0) THEN
	   FP(d(2)+1:nx,d(2)+1:nx) = -F(d(2)+1:nx,d(2)+1:nx,S(1,5))
	   DO 20 I = d(2)+1,nx
20	   FP(I,I) = 1.D0+FP(I,I)
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
C	 CALL G05EZF(XP(1,d(2)+1:nx),nx-d(2),WORK,(nx+2)*(nx+1)/2,IFAIL)
       CALL setgmn(Xdd(1,d(2)+1:nx),Pdd(1,d(2)+1:nx,d(2)+1:nx),nx-d(2),
     #             nx-d(2),WORK(1:(nx-d(2)+2)*(nx-d(2)+1)/2))
       CALL genmn(WORK(1:(nx-d(2)+2)*(nx-d(2)+1)/2),STATE(1,d(2)+1:nx),
     #            WORK1(1:nx-d(2)))
      ENDIF

C DRAW u(1)+
	DO 35 J = 1,nu
C35   UP(J) = G05DDF(0.0D0,1.D0)
35    UP(J) = gennor(0.0D0,1.D0)


C COMPUTE y(1)+
	DO 36 K = 1,ny
36	ykP(1,K) = SUM(H(K,1:nx,S(1,2))*XP(1,1:nx))
     #         + SUM(G(K,1:nu,S(1,3))*UP(1:nu))
     #         + SUM(c(K,1:nz,S(1,1))*yk(1,ny+1:ny+nz))

	DO 100 it = 2,nobs
C DRAW u+ ~ N(0,I)
	 DO 40 J = 1,nu
C40    UP(J) = G05DDF(0.0D0,1.D0)
40     UP(J) = gennor(0.0D0,1.D0)

C COMPUTE x+
	 DO 50 K = 1,nx
50	 XP(it,K) = a(K,S(it,4))+SUM(F(K,1:nx,S(it,5))*XP(it-1,1:nx))
     #          + SUM(R(K,1:nu,S(it,6))*UP(1:nu))

C COMPUTE yk+
	 DO 60 K = 1,ny
60	 ykP(it,K) = SUM(H(K,1:nx,S(it,2))*XP(it,1:nx))
     #           + SUM(G(K,1:nu,S(it,3))*UP(1:nu))
     #           + SUM(c(K,1:nz,S(it,1))*yk(it,ny+1:ny+nz))

100   CONTINUE

C KALMAN SMOOTHING RECURSIONS
      a(:,:)   = 0.D0
      c(:,:,:) = 0.D0
	CALL KS2(nobs,d,ny,nz,nx,nu,ns,S,yk-ykP,c,H,G,a,F,R,XS)

C SETTING THE STATE
	STATE(1:nobs,:) = XP(1:nobs,:) + XS(1:nobs,:)

      DEALLOCATE (R,c,H,G,a,F,ykP,XP,XS,Xdd,Pdd,WORK,FP,WORK1,UP)

	RETURN
	END

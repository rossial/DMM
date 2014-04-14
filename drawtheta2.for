C ----------------------------------------------------------------------
C DRAWTHETA2 draws model parameters theta (no missing values)
C Developed by A.Rossi, C.Planas and G.Fiorentini
C p(theta|Y,S,z)  propto  p(theta|Y,S,z)
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
C  OUTPUT:
C   theta (nt x 1)
C   NEVAL # of runs of the Kalman Filter (KF) / MH acceptance prob
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
C ----------------------------------------------------------------------
	SUBROUTINE DRAWTHETA2(nobs,d,ny,nz,nx,nu,nv,ns,nt,INFOS,INDT,
	1                      MEDT,SIGT,yk,Z,thetaprior,tipo,pdll,
     2                      theta0,theta,NEVAL)
C INPUT
	CHARACTER*1 fittizia
	POINTER (pdll,fittizia)
	INTEGER nobs,d(2),ny,nz,nx,nu,nv,ns(6),nt,Z(nobs),
	1 INFOS(9,6),INDT(nt+2)
	DOUBLE PRECISION yk(nobs,ny+nz),thetaprior(nt,4),theta0(nt)
      DOUBLE PRECISION MEDT(max(1,INDT(nt+1))),
	1 SIGT(max(1,INDT(nt+1)),max(1,INDT(nt+1)))
	CHARACTER*2 tipo(nt)
C OUTPUT
	INTEGER NEVAL(nt)
	DOUBLE PRECISION theta(nt)
C LOCALS
	INTEGER I,J,JJ,NN,IFAIL,SEQ(nv)
	INTEGER it,S(nobs,6)
	DOUBLE PRECISION PN,PO,QN,QO,PA,v,AG,WORK((nt+2)*(nt+1)/2),
	1 SEGA(nt),WORK2(nt)
    	DOUBLE PRECISION genunf,ptheta2,prior,mvnpdf
C     DOUBLE PRECISION G05CAF

	IF (nv.GT.0) THEN
	 DO 1 I = 1,nobs
1	 CALL INT2SEQ(Z(I),nv,INFOS,SEQ(1:nv),S(I,1:6))
	ELSE
	 S(1:nobs,:) = 1
      ENDIF

C --------------------------------------------
C DRAW theta1 from p(theta1|Y,S)
C p(theta1|Y,S) prop p(Y|theta1,S) x p(theta1)
C --------------------------------------------
	NN = 0
	JJ = INDT(nt+1)
      IF (JJ.GT.1) THEN ! Draw a set of theta via AMH
C MH proposal:
C theta ~ (1-beta)*N(MEDT,2.38^2*SIGT/JJ) + beta*N(MEDT,.1^2*I/JJ)
C7777	 v = G05CAF(v) ! Sampling U(0,1)
7777   v = genunf(0.D0,1.D0)    ! Sampling U(0,1)
	 IF (v.LE..95) THEN
	  IFAIL = -1
c	  CALL G05EAF(MEDT,JJ,2.38**2*SIGT/DFLOAT(JJ),JJ,1.D-14,
c     1              WORK,(nt+2)*(nt+1)/2,IFAIL)
c	  CALL G05EZF(SEGA(1:JJ),JJ,WORK,(nt+2)*(nt+1)/2,IFAIL)
        CALL setgmn(MEDT,2.38**2*SIGT/DFLOAT(JJ),JJ,JJ,
     #              WORK(1:(JJ+2)*(JJ+1)/2))
        CALL genmn(WORK(1:(JJ+2)*(JJ+1)/2),SEGA(1:JJ),WORK2(1:JJ))
	 ELSE
	  DO I = 1,JJ
c	   CALL G05EAF(MEDT(I),1,1D-2/DFLOAT(JJ),1,1.D-14,
c	1               WORK,(nt+2)*(nt+1)/2,IFAIL)
c	   CALL G05EZF(SEGA(I),1,WORK,(nt+2)*(nt+1)/2,IFAIL)
         CALL setgmn(MEDT(I),1D-2/DFLOAT(JJ),1,1,WORK(1:3))
         CALL genmn(WORK(1:3),SEGA(I),WORK2(1))
        ENDDO
       ENDIF
	 NN = NN + 1
C CHEK theta
	 DO I = 1,JJ
        IF ((SEGA(I).LT.thetaprior(INDT(I),3)).OR.
     #     (SEGA(I).GT.thetaprior(INDT(I),4))) THEN
	   IF (NN.LE.1000) THEN
	     GOTO 7777
	   ELSE
#ifdef DYNARE
          WRITE(*,*) ' '
          WRITE(*,*) 'Reduce skcriterium or use Slice sampling'
          WRITE(*,*) 'Program aborting'
#else
	    type *, ' '
	    type *, 'Reduce skcriterium or use Slice sampling'
          type *, 'Program aborting'
#endif
	    PAUSE
	    STOP
	   ENDIF
	  ENDIF
	 END DO
C f(theta1(new)|theta2,S,y)
	 theta(INDT(1:JJ)) = SEGA(1:JJ)
	 PN = PTHETA2(INDT(1),nobs,d,ny,nz,nx,nu,ns,nt,S,yk,theta,
	1              thetaprior(INDT(1),:),tipo(INDT(1)),pdll)
	 DO I =2,JJ
	  PN = PN + PRIOR(theta(INDT(I)),thetaprior(INDT(I),:),
     #                  tipo(INDT(I)))
	 ENDDO
C q(theta1(new))
	 AG = 1.D0
	 DO I = 1,JJ
	  AG = AG*mvnpdf(SEGA(I),MEDT(I),1D-2/DFLOAT(JJ),1)
	 END DO
	 QN = .95D0*mvnpdf(SEGA(1:JJ),MEDT,2.38**2*SIGT/DFLOAT(JJ),JJ)
     +    + .05D0*AG

C f(theta1(old)|theta2,S,y)
	 PO = PTHETA2(INDT(1),nobs,d,ny,nz,nx,nu,ns,nt,S,yk,theta0,
	1              thetaprior(INDT(1),:),tipo(INDT(1)),pdll)
	 DO I =2,JJ
	  PO = PO + PRIOR(theta0(INDT(I)),thetaprior(INDT(I),:),
     #                  tipo(INDT(I)))
	 ENDDO
C q(theta1(old))
	 AG = 1.D0
	 DO I = 1,JJ
	  AG = AG*mvnpdf(theta0(INDT(1:JJ)),MEDT(I),1D-2/DFLOAT(JJ),1)
	 END DO
	 QO = .95D0*mvnpdf(theta0(INDT(1:JJ)),MEDT,2.38**2*SIGT/DFLOAT(JJ)
     #	  ,JJ) + .05D0*AG

c	 v = G05CAF(v) ! Sampling from U(0,1)
       v = genunf(0.D0,1.D0)    ! Sampling U(0,1)
	 PA = DEXP(PN-PO)*QO/QN
	 IF (v.GT.MIN(1.D0,PA)) THEN
	  theta(INDT(1:JJ)) = theta0(INDT(1:JJ))
	  NEVAL(INDT(1)) = 0
	 ELSE
	  theta0(INDT(1:JJ)) = theta(INDT(1:JJ))
	  NEVAL(INDT(1)) = 1
	 ENDIF
	ENDIF

C The rest of theta by SLICE
	DO J = 1,INDT(nt+2)-JJ
	 it = INDT(JJ+J)
	 CALL SLICE2(it,nobs,d,ny,nz,nx,nu,ns,nt,S,yk,theta0,
	1          thetaprior(it,:),tipo(it),pdll,NEVAL(it),theta(it))
       theta0(it) = theta(it)
	END DO

	RETURN
	END

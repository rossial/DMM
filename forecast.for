C ----------------------------------------------------------------------------
C FORECAST simulates y(T+1),...,y(T+nf),x(T+1),...,x(T+nf),S(T+1),...,S(T+nf)
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
C (nf x (ny+nx+1)) FORE: y(T+k),x(T+k),S(T+k)
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
C ----------------------------------------------------------------------------	
	SUBROUTINE FORECAST(zk,nf,ny,nz,nx,nu,nv,ns,nstot,nt,np,
	1                    theta,psi,INFOS,Z,STATE,pdll,FORE)

	USE dfwin
	INTERFACE
	 SUBROUTINE DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R) 
	 INTEGER ny,nz,nx,nu,ns(6),nt
	 DOUBLE PRECISION theta(nt)
	 DOUBLE PRECISION c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6))
	 END SUBROUTINE
	END INTERFACE
	CHARACTER*1 fittizia
	POINTER (pdll,fittizia)  ! ASSOCIATE  pointer pdll alla DLL ad una varibile fittizia
	POINTER (pdesign,DESIGN) 

C INPUT
	INTEGER nf,ny,nz,nx,nu,nv,nt,np(3),ns(6),nstot,Z,INFOS(9,6)
	DOUBLE PRECISION zk(nf,nz),theta(nt),psi(np(1)),STATE(nx)

C OUTPUT
	DOUBLE PRECISION FORE(nf,ny+nx+1)
C LOCALS
	INTEGER ISEQ,SEQ(nv),inf,I,IFAIL
	INTEGER S(6)
      DOUBLE PRECISION,ALLOCATABLE::R(:,:,:),c(:,:,:),H(:,:,:),
	1 G(:,:,:),a(:,:),F(:,:,:)
	DOUBLE PRECISION,ALLOCATABLE:: P1(:,:),P2(:,:),P3(:,:),P4(:,:),
	1 P5(:,:),P6(:,:),PMAT(:,:),PE(:)
	DOUBLE PRECISION U,AUX
      DOUBLE PRECISION MED(nu)
C	DOUBLE PRECISION WORKU((nu+2)*(nu+1)/2)

C EXTERNAL FUNCTIONS      
      DOUBLE PRECISION genunf,gennor
C EXTERNAL SUBROUTINES      
      EXTERNAL DESIGNZ,PPROD,ERGODIC,INT2SEQ

      
	ALLOCATE(R(nx,nu,ns(6)),c(ny,max(nz,1),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)))
 
C Call DESIGN
	pdesign = getprocaddress(pdll, "design_"C)
	CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)

	IF (nv.GT.0) THEN
	 ALLOCATE (P1(INFOS(8,1),INFOS(8,1)),
	1 P2(INFOS(8,2),INFOS(8,2)),P3(INFOS(8,3),INFOS(8,3)),
     2 P4(INFOS(8,4),INFOS(8,4)),P5(INFOS(8,5),INFOS(8,5)),
     3 P6(INFOS(8,6),INFOS(8,6)),PMAT(nstot,nstot),PE(nstot))

	 CALL DESIGNZ(nv,np(1),psi,INFOS,P1,P2,P3,P4,P5,P6)
C PALL(i,j) = Pr[Z(t+1)=i|Z(t)=j], Z = S1 x S2 x ... x Snv  
	 CALL PPROD(nv,nstot,INFOS,P1,P2,P3,P4,P5,P6,PMAT)
C ERGODIC solves PE: PE*(I-P') = 0
       CALL ERGODIC(nstot,PMAT,PE)	     	  
	ENDIF

C DRAW Z(T+1) ~ Pr[Z(T+1)|Z(T)] 	
	S(1:6) = 1
	IF (nv.GT.0) THEN
C	 U = G05CAF(U) ! Sampling U(0,1)  
       U = genunf(0.D0,1.D0)
	 ISEQ = 1
	 AUX  = PMAT(ISEQ,Z)
	 DO 10 WHILE (AUX.LT.U) 
	 ISEQ = ISEQ + 1 
10	 AUX  = AUX  + PMAT(ISEQ,Z)
	 FORE(1,nx+ny+1) = ISEQ
	 CALL INT2SEQ(ISEQ,nv,INFOS,SEQ,S)	  
	ENDIF

C DRAW x(T+1) ~ f(x(T+1)|x(T),Z(T+1))
	IFAIL = -1
	DO I = 1,nu
       MED(I) = gennor(0.D0,1.D0)          
	END DO
	DO 20 I=1,nx
20	FORE(1,ny+I) = a(I,S(4)) + SUM(F(I,1:nx,S(5))*STATE(1:nx))
     +             + SUM(R(I,1:nu,S(6))*MED(1:nu))

C DRAW yk(T+1) ~ f(yk(T+1)|x(T+1),Z(T+1),zk(T+1))
	DO 30 I=1,ny
30	 FORE(1,I) = SUM(H(I,1:nx,S(2))*FORE(1,ny+1:ny+nx))
     +           + SUM(c(I,1:nz,S(1))*zk(1,1:nz))
     +           + SUM(G(I,1:nu,S(3))*MED(1:nu))

	DO 100 inf = 2,nf
C DRAW Z(T+inf) ~ Pr[Z(T+inf)|Z(T+inf-1)] 	
	 IF (nv.GT.0) THEN
C	  U = G05CAF(U) ! Sampling U(0,1)  
        U = genunf(0.D0,1.D0)
	  ISEQ = 1
	  AUX  = PMAT(ISEQ,FORE(inf-1,nx+ny+1))
	  DO 40 WHILE (AUX.LT.U) 
	  ISEQ = ISEQ + 1 
40	  AUX  = AUX  + PMAT(ISEQ,FORE(inf-1,nx+ny+1))
	  FORE(inf,nx+ny+1) = ISEQ
	  CALL INT2SEQ(ISEQ,nv,INFOS,SEQ,S)	  
	 ENDIF

C DRAW x(T+inf) ~ f(x(T+inf)|x(T+inf-1),Z(T+inf))
	 DO I = 1,nu
        MED(I) = gennor(0.D0,1.D0)   
	 END DO
	 DO 50 I=1,nx
50	 FORE(inf,ny+I) = a(I,S(4)) 
     +                + SUM(F(I,1:nx,S(5))*FORE(inf-1,ny+1:ny+nx))
     +                + SUM(R(I,1:nu,S(6))*MED(1:nu))

C DRAW y(T+inf) ~ f(y(T+1)|x(T+inf),Z(T+inf),zk(T+inf))
	 DO 60 I=1,ny
60	 FORE(inf,I) = SUM(H(I,1:nx,S(2))*FORE(inf,ny+1:ny+nx))
     +             + SUM(c(I,1:nz,S(1))*zk(inf,1:nz))
     +             + SUM(G(I,1:nu,S(3))*MED(1:nu))

100   CONTINUE

      DEALLOCATE (R,c,H,G,a,F)
	IF(nv.GT.0) DEALLOCATE(P1,P2,P3,P4,P5,P6,PMAT,PE) 
	
	RETURN
	END

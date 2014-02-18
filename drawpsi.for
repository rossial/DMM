C ---------------------------------------------------------------------- 
C DRAWS parameters PSI conditionally on S            
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C  
C p(theta,psi|Y,S,z)  propto  p(theta|Y,S,z) x p(psi|S)   
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
C FURTHER INPUT: 
C
C    nobs: # of observations
C       d: order of integration of the system
C      nv: # of discrete latent variables (S1,S2,...)              
C      nt: dimension of parameter theta
C      Z: discrete latent var
C  theta0: previous value of theta 
C INFOS (9 x 6):
C by cols: S1,S2,...,Snv; with nv <=6
C by row: the 1st contains the # of matrices affected by Si  
C         the 2nd-3rd etc point to c (1),H (2),G (3),a (4),F (5),R (6)
C         the 8-th row contains the # of states
C         the 9-th row spec. the dynamic of S
C
C  OUTPUT:     
C   PSI (np(1) x 1) 
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
	SUBROUTINE DRAWPSI(nobs,nv,np,INFOS,Z,psiprior,psi0,psi)
	           
C INPUT
	INTEGER nobs,nv,np(3),Z(nobs),INFOS(9,6)
	DOUBLE PRECISION psiprior(np(2),np(3)),psi0(np(1))
C OUTPUT
	DOUBLE PRECISION psi(np(1))
C LOCALS 
	INTEGER I,K,ii,jj,NN,NSI,IFAIL,SEQ(nobs,nv)
	INTEGER,ALLOCATABLE:: NIJ(:,:)
	INTEGER S(nobs,6)
	DOUBLE PRECISION, ALLOCATABLE:: P1(:,:),P2(:,:),P3(:,:),P4(:,:),	
     1 P5(:,:),P6(:,:),PENEW(:),PEOLD(:),GAM(:)
	DOUBLE PRECISION uv,v,AG
      DOUBLE PRECISION ranf,gengam
C	DOUBLE PRECISION G05CAF
	
	ALLOCATE(P1(INFOS(8,1),INFOS(8,1)),
	1 P2(INFOS(8,2),INFOS(8,2)),P3(INFOS(8,3),INFOS(8,3)),
     2 P4(INFOS(8,4),INFOS(8,4)),P5(INFOS(8,5),INFOS(8,5)),
     3 P6(INFOS(8,6),INFOS(8,6)))
	
	DO 1 I = 1,nobs
1	CALL INT2SEQ(Z(I),nv,INFOS,SEQ(I,1:nv),S(I,1:6))

	psi(:) = 0.D0
	NN = 0
	K  = 0
	DO 100 I = 1,nv 	   
	  NSI = INFOS(8,I) ! # of states for SI
	  ALLOCATE(GAM(NSI))
	  IF (INFOS(9,I).EQ.1) THEN      ! S~IID
C SAMPLING FROM DIRICHLET
	   DO 5 ii = 1,NSI	   
         AG = SUM(ABS((SEQ(1:nobs,I).EQ.ii)))+psiprior(K+1,ii) 
	   IFAIL = -1	      
C 5      CALL G05FFF(AG,1.D0,1,GAM(ii),IFAIL) 	    
5        GAM(ii) = gengam(1.D0,AG)           
	   psi(NN+1:NN+NSI-1)  = GAM(1:NSI-1)/SUM(GAM(1:NSI))	     
         psi0(NN+1:NN+NSI-1) = psi(NN+1:NN+NSI-1)
	   K  = K + 1       
         NN = NN + NSI-1	    	  	   
	  ELSEIF (INFOS(9,I).EQ.2) THEN  ! S~Markov	   
	   ALLOCATE(NIJ(NSI,NSI),PENEW(NSI),PEOLD(NSI))
	   DO 50 jj = 1,NSI 
	    DO 10 ii = 1,NSI
10        NIJ(ii,jj) = SUM(ABS((SEQ(2:nobs,I).EQ.ii).AND.
     #		          (SEQ(1:nobs-1,I).EQ.jj)))
C SAMPLING FROM DIRICHLET
		DO 20 ii = 1,NSI
		AG = NIJ(ii,jj) + psiprior(K+1,ii) 
	    IFAIL = -1	      
C20        CALL G05FFF(AG,1.D0,1,GAM(ii),IFAIL) 
20        GAM(ii) = gengam(1.D0,AG)          
	    psi(NN+1:NN+NSI-1) = GAM(1:NSI-1)/SUM(GAM(1:NSI))	     
          K  = K + 1       
50	   NN = NN + NSI-1	    
C METROPOLIS TO ADJUST INITIAL CONDITION P(S(1)=0|p11,p12,...)			
	   CALL DESIGNZ(nv,np(1),psi,INFOS,P1,P2,P3,P4,P5,P6) !new           
         IF (I.EQ.1) THEN 
	    CALL ERGODIC(NSI,P1,PENEW(1:NSI)) 
	   ELSEIF (I.EQ.2) THEN
	    CALL ERGODIC(NSI,P2,PENEW(1:NSI)) 
	   ELSEIF (I.EQ.3) THEN
	    CALL ERGODIC(NSI,P3,PENEW(1:NSI)) 
	   ELSEIF (I.EQ.4) THEN 
	    CALL ERGODIC(NSI,P4,PENEW(1:NSI)) 
	   ELSEIF (I.EQ.5) THEN 
	    CALL ERGODIC(NSI,P5,PENEW(1:NSI)) 
	   ELSE
	    CALL ERGODIC(NSI,P6,PENEW(1:NSI)) 
	   ENDIF
	   CALL DESIGNZ(nv,np(1),psi0,INFOS,P1,P2,P3,P4,P5,P6)  !old
         IF (I.EQ.1) THEN 
	    CALL ERGODIC(NSI,P1,PEOLD(1:NSI)) 
	   ELSEIF (I.EQ.2) THEN
	    CALL ERGODIC(NSI,P2,PEOLD(1:NSI)) 
	   ELSEIF (I.EQ.3) THEN
	    CALL ERGODIC(NSI,P3,PEOLD(1:NSI)) 
	   ELSEIF (I.EQ.4) THEN 
	    CALL ERGODIC(NSI,P4,PEOLD(1:NSI)) 
	   ELSEIF (I.EQ.5) THEN 
	    CALL ERGODIC(NSI,P5,PEOLD(1:NSI)) 
	   ELSE
	    CALL ERGODIC(NSI,P6,PEOLD(1:NSI)) 
	   ENDIF
	   uv = min(1.D0,PENEW(SEQ(1,I))/PEOLD(SEQ(1,I))) 
C	   v = G05CAF(v) 
         v = ranf()  ! U(0,1)
	   IF (v.GT.uv) THEN
	    psi(NN-NSI*(NSI-1)+1:NN) = psi0(NN-NSI*(NSI-1)+1:NN)
	   ENDIF
	   psi0(NN-NSI*(NSI-1)+1:NN) = psi(NN-NSI*(NSI-1)+1:NN)
         DEALLOCATE(NIJ,PENEW,PEOLD) 
	  ENDIF	   
100   DEALLOCATE(GAM) 
	DEALLOCATE(P1,P2,P3,P4,P5,P6)
	
	RETURN
	END
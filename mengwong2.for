C -------------------------------------------------------------------
C MENGWONG2 (no missing values) computes the Marginal Lilkelihood estimates as deteiled 
C by Meng and Wong, Statistica Sinica, 1996 
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C
C OUTPUT:      
C  MLMW(:,1) all parameters, 
C  MLMW(:,2) non-var params
C  MLMW(1,:) no iteration, 
C  MLMW(2,:) SD, 
C  MLMW(3,:) 10 iterations
C
C Remarks:
C  NPAR is total # of params, 
C  NPARD = NPAR - #Variances 
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
C      
C In addition, as a special exception, the copyright holders give
C permission to link the code of portions of this program with the
C NAG Fortran library under certain conditions as described in each
C individual source file, and distribute linked combinations including
C the two.
C
C You must obey the GNU General Public License in all respects for all
C of the code used other than NAG Fortran library. If you modify file(s)
C with this exception, you may extend this exception to your
C version of the file(s), but you are not obligated to do so. If
C you do not wish to do so, delete this exception statement from
C your version. If you delete this exception statement from all
C source files in the program, then also delete it here.      
C -------------------------------------------------------------------
	SUBROUTINE MENGWONG2(G,nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np,
	1                     INFOS,yk,gibpar,gibZ,thetaprior,psiprior,
     2                     tipo,pdll,MLSTART,MLMW)

C INPUT
      INTEGER G,nobs,d(2),ny,nz,nx,nu,nv,ns(6),nstot,nt,np(3),
	1 INFOS(9,6),gibZ(G,nobs) 
	DOUBLE PRECISION yk(nobs,ny+nz),gibpar(G,nt+np(1)),
	1 thetaprior(nt,4),psiprior(np(2),np(3)),MLSTART
	CHARACTER*2 tipo(nt) 
	POINTER (pdll,fittizia) 

C OUTPUT
	DOUBLE PRECISION MLMW(2,2)

C LOCALS
	INTEGER NPAR,I,J,K,IG,NPOS(nt+np(1)),IFAIL,NQ,ISEQ,ISEQ0,SEQ(nv),
	1 IS(nobs,6),NIM,NI,IND(1),NPARTH,NN,NSI,II,JJ
	DOUBLE PRECISION,ALLOCATABLE::MAT(:,:),VQN(:,:),VQD(:,:),
	1 VHN(:,:),VHD(:,:)
	DOUBLE PRECISION parm(nt),SIGM(nt,nt),
	1 COM(nt+1,nt),ISIGM(nt,nt),par(nt+np(1)),SEGA(nt+np(1)),
     2 ub(nt),lb(nt),R3((nt+1)*(nt+2)/2),WORK(nt)
	DOUBLE PRECISION,ALLOCATABLE:: PTR(:,:,:),PMAT(:,:),PE(:),GAM(:),
	1 ALPHA(:,:),MOM(:,:)
	DOUBLE PRECISION P1(INFOS(8,1),INFOS(8,1)),
	1 P2(INFOS(8,2),INFOS(8,2)),P3(INFOS(8,3),INFOS(8,3)),
     2 P4(INFOS(8,4),INFOS(8,4)),P5(INFOS(8,5),INFOS(8,5)),
     3 P6(INFOS(8,6),INFOS(8,6))
	DOUBLE PRECISION PTHETA2,PRIOR,PRIORDIR,genunf,gengam
	DOUBLE PRECISION Ppar(nt+np(1)),Fpar,PS,QS,QPSI,C,DET,TRC,A0
	DOUBLE PRECISION ERRM,ERR,U,AUX,INDC(1),MUC,SS(2,2),MWNUM,MWDEN
	DOUBLE PRECISION ZERO,ONE,PI
	DATA ZERO/0.0D0/,ONE/1.0D0/,PI/3.141592653589793D0/

      PAR(:) = GIBPAR(1,:) ! set constant values

	NPARTH = 0
	DO I = 1,nt  
	 IF (GIBPAR(1,I).NE.GIBPAR(2,I)) THEN
	  NPARTH = NPARTH + 1
        NPOS(NPARTH) = I
	 ENDIF
      ENDDO
	DO I = 1,np(1)
       NPOS(NPARTH+I) = nt+I
      ENDDO
	NPAR=NPARTH+np(1) 
      Ppar(:) = 0.D0
	parm(:) = ZERO
	DO I = 1,NPARTH
	 parm(I) = SUM(gibpar(:,NPOS(I)))/DFLOAT(G)    
      ENDDO

	NQ = 0 
	CALL NEWEYWESTCOV2(G,NPARTH,NQ,gibpar(:,NPOS(1:NPARTH)),
	1  parm(1:NPARTH),SIGM(1:NPARTH,1:NPARTH)) ! THETA Var-covar

	IF (nv.GT.0) THEN	  
	 ALLOCATE(PTR(nobs,nstot,nstot),PMAT(nstot,nstot),PE(nstot))

C Transition prob for QS
	 DO I = 1,nstot-1
	  PTR(1,I,1) = SUM(ABS(gibZ(1:G,1).EQ.I))/DFLOAT(G)
       ENDDO
       PTR(1,nstot,1) = ONE-SUM(PTR(1,1:nstot-1,1))
	
	 DO 50 K = 2,nobs
	  DO 50 I = 1,nstot-1
	   DO 50 J = 1,nstot
	    COM(1,1) = SUM(ABS(gibZ(1:G,K-1).EQ.J))
	    IF (COM(1,1).GT.ZERO) THEN
	     PTR(K,I,J) = SUM(ABS((gibZ(1:G,K).EQ.I).AND.
     #		          (gibZ(1:G,K-1).EQ.J)))/COM(1,1)
	    ELSE
	     PTR(K,I,J) = ONE/DFLOAT(nstot)
	    ENDIF
50	   PTR(K,nstot,J) = ONE-SUM(PTR(K,1:nstot-1,J))

C Mean and Var of PSI 
	 ALLOCATE (ALPHA(np(2),np(3)),MOM(np(1),2))
	 DO I=1,np(1)
	  MOM(I,1) = SUM(gibpar(:,nt+I))/DFLOAT(G)
	  MOM(I,2) = SUM(gibpar(:,nt+I)**2)/DFLOAT(G)
	  MOM(I,2) = MOM(I,2)-MOM(I,1)**2
	 ENDDO
C Hyperparameters of Dirichelt for Q(PSI)	
C Mothod of Moments: a0 = m1(1-m1)/V1+1, ai = mi*a0, i=1,2,..,N
	 NN = 0
	 K  = 0
	 DO I = 1,nv 	   
	  NSI = INFOS(8,I)           ! # of states for S
	  IF (INFOS(9,I).EQ.1) THEN  ! S~IID	   
	   A0 = MOM(NN+1,1)*(1.D0-MOM(NN+1,1))/MOM(NN+1,2)+1.D0 !alpha0
	   DO  ii = 1,NSI-1	   
          ALPHA(K+1,ii) = MOM(NN+ii,1)*A0
         ENDDO
	   ALPHA(K+1,NSI) = A0-SUM(ALPHA(K+1,1:NSI-1)) 
	   K  = K + 1       
         NN = NN + NSI-1	    	  	   
	  ELSEIF (INFOS(9,I).EQ.2) THEN  ! S~Markov	   
	   DO jj = 1,NSI 
	    A0 = MOM(NN+1,1)*(1.D0-MOM(NN+1,1))/MOM(NN+1,2)+1.D0 !alpha0
	    DO  ii = 1,NSI-1	   
           ALPHA(K+1,ii) = MOM(NN+ii,1)*A0
          ENDDO
	    ALPHA(K+1,NSI) = A0-SUM(ALPHA(K+1,1:NSI-1))
          K  = K + 1       
	    NN = NN + NSI-1	    
         ENDDO
	  ENDIF
	 ENDDO
	ENDIF
C Importance sampling
C Sample THIS from N(THHAT,SIGHAT) with boundaries 
C Evaluate p(THIS) ~ N(THHAT,SIGHAT)
	NIM  = 1000000
	ERRM = 1.D-8

C Normalization constants TRC and TRCD
      lb(1:NPARTH) = thetaprior(NPOS(1:NPARTH),3)
      ub(1:NPARTH) = thetaprior(NPOS(1:NPARTH),4)
	CALL mvncdf(lb(1:NPARTH),ub(1:NPARTH),parm(1:NPARTH),
	1   SIGM(1:NPARTH,1:NPARTH),NPARTH,ERRM,NIM,TRC,ERR,NI) 

C Inverse SIGM  & det for NPARTH 
	COM(1:NPARTH,1:NPARTH) = SIGM(1:NPARTH,1:NPARTH)
	IFAIL = -1
	CALL F01ADF(NPARTH,COM(1:NPARTH+1,1:NPARTH),NPARTH+1,IFAIL) ! Inverse var-covar
	DO 60 I=1,NPARTH
	 ISIGM(I,I) = COM(I+1,I)
	 DO 60 J=1,I-1	
	  ISIGM(I,J) = COM(I+1,J)
60	 ISIGM(J,I) = ISIGM(I,J)

	COM(1:NPARTH,1:NPARTH) = SIGM(1:NPARTH,1:NPARTH)
      IFAIL = -1
      CALL F03ABF(COM(1:NPARTH,1:NPARTH),NPARTH,NPARTH,DET,
	1 WORK(1:NPARTH),IFAIL)
      C = (2.D0*PI)**(-.5D0*NPARTH)/DSQRT(DET) ! constant
	
	ALLOCATE (MAT(G,2),VHN(G,2),VHD(G,2),VQN(G,2),VQD(G,2))
	QS = ONE
	PS = ONE
	IS(:,:) = 1
	DO 200 IG = 1,G
C SAMPLING THETA 
       SEGA(:)  = -1.D0
 	 IFAIL    = -1
	 IND(1)   =  0
	 INDC(1)  = -1.D0
	 DO WHILE (INDC(1).LT.ZERO)
	  INDC(1) = ZERO
	  IND(1)  = IND(1) + 1
	  IF (IND(1).GT.G) EXIT
c	  CALL G05EAF(parm(1:NPARTH),NPARTH,SIGM(1:NPARTH,1:NPARTH),
c	1              NPARTH,EPS,R3,(NPARTH+1)*(NPARTH+2)/2,IFAIL)	   
c	  CALL G05EZF(SEGA(1:NPARTH),NPARTH,R3,(NPARTH+1)*(NPARTH+2)/2,
c	1              IFAIL)
        COM(1:NPARTH,1:NPARTH) = SIGM(1:NPARTH,1:NPARTH)
        CALL setgmn(parm(1:NPARTH),COM(1:NPARTH,1:NPARTH),NPARTH,
     1              NPARTH,R3(1:(NPARTH+2)*(NPARTH+1)/2))
        CALL genmn(R3(1:(NPARTH+2)*(NPARTH+1)/2),SEGA(1:NPARTH),
     1             WORK(1:NPARTH))
	  DO  I=1,NPARTH
	   IF (SEGA(I).LT.thetaprior(NPOS(I),3)) INDC(1)=-1
         IF (SEGA(I).GT.thetaprior(NPOS(I),4)) INDC(1)=-2
        ENDDO
       END DO
C SAMPLING PSI from Dirichlet(ALPHA)
	 NN = NPARTH
	 K  = 0
	 DO 70 I = 1,nv 	   
	  NSI = INFOS(8,I) ! # of states for SI
	  ALLOCATE(GAM(NSI))
	  IF (INFOS(9,I).EQ.1) THEN ! S~IID
	   DO  ii = 1,NSI	              
	    IFAIL = -1	      
C	    CALL G05FFF(ALPHA(K+1,ii),1.D0,1,GAM(ii),IFAIL) 
          GAM(ii) = gengam(1.D0,ALPHA(K+1,ii))
         ENDDO
	   SEGA(NN+1:NN+NSI-1) = GAM(1:NSI-1)/SUM(GAM(1:NSI))	     
	   K  = K + 1       
         NN = NN + NSI-1	    	  	   
	  ELSEIF (INFOS(9,I).EQ.2) THEN  ! S~Markov	   
	   DO jj = 1,NSI
	    DO ii = 1,NSI
		 IFAIL = -1	      
C	     CALL G05FFF(ALPHA(K+1,ii),1.D0,1,GAM(ii),IFAIL) 	    
           GAM(ii) = gengam(1.D0,ALPHA(K+1,ii))
          ENDDO
	    SEGA(NN+1:NN+NSI-1) = GAM(1:NSI-1)/SUM(GAM(1:NSI))	     
          K  = K + 1       
	    NN = NN + NSI-1
	   ENDDO
        ENDIF
70     DEALLOCATE(GAM)

C SAMPLING S
	 IF (nv.GT.0) THEN
	  CALL DESIGNZ(nv,np(1),SEGA(NPARTH+1:NPAR),INFOS,
	1               P1,P2,P3,P4,P5,P6)  
C PMAT(i,j) = Pr[Z(t+1)=i|Z(t)=j], Z = S1 x S2 x ... x Snv  
	  CALL PPROD(nv,nstot,INFOS,P1,P2,P3,P4,P5,P6,PMAT)
C ERGODIC solves PE: PE*(I-P') = 0
        CALL ERGODIC(nstot,PMAT,PE)
C S(1)
C	  U = G05CAF(U) ! Sampling from U(0,1)  
        U = genunf(0.D0,1.D0)
	  ISEQ = 1
	  AUX  = PTR(1,ISEQ,1)
	  DO 80 WHILE (AUX.LT.U) 
	   ISEQ = ISEQ + 1 
80	  AUX  = AUX  + PTR(1,ISEQ,1)
	  CALL INT2SEQ(ISEQ,nv,INFOS,SEQ,IS(1,:))
	  QS = PTR(1,ISEQ,1)
	  PS = PE(ISEQ) ! P(S1)
	  ISEQ0 = ISEQ
C S(2),...,S(nobs)	 
	  DO 90 K = 2,nobs
C	   U = G05CAF(U) ! Sampling from U(0,1)  
         U = genunf(0.D0,1.D0)
	   ISEQ = 1
	   AUX  = PTR(K,ISEQ,ISEQ0)
	   DO 85 WHILE (AUX.LT.U) 
	   ISEQ = ISEQ + 1 
85	   AUX  = AUX  + PTR(K,ISEQ,ISEQ0)
	   CALL INT2SEQ(ISEQ,nv,INFOS,SEQ,IS(K,:))
	   QS = QS*PTR(K,ISEQ,ISEQ0)
	   PS = PS*PMAT(ISEQ,ISEQ0)
90	  ISEQ0 = ISEQ 
	 ENDIF
	
C QUADRATIC FORM FOR for THETA
       DO 91 I = 1,NPARTH
91     COM(I,1) = SUM((SEGA(1:NPARTH)-parm(1:NPARTH))*ISIGM(1:NPARTH,I))  	   
	 MUC = SUM(COM(1:NPARTH,1)*(SEGA(1:NPARTH)-parm(1:NPARTH)))

C	 VQN(IG,1) = QS*C*DEXP(-.5D0*MUC)/TRC       
	 
	 par(NPOS(1:NPARTH+np(1))) = SEGA(1:NPARTH+np(1)) ! (THETA,PSI)

C PRIOR for THETA
	 DO 92 I = 2,NPARTH
92	 Ppar(I) = PRIOR(par(NPOS(I)),thetaprior(NPOS(I),:),tipo(NPOS(I)))  

C PRIOR for PSI and Q(PSI)~Dirichlet(a1,a2,...,aN)
	 QPSI = 0.D0
	 NN = NPARTH
	 K  = 0
	 DO 100 J = 1,nv
	  NSI = INFOS(8,J)  
	  IF(INFOS(9,J).EQ.1) THEN ! S~IID
	   Ppar(NPARTH+K+1) = PRIORDIR(par(NPOS(NN+1:NN+NSI-1)),
	1                               psiprior(K+1,1:NSI),NSI)
	   QPSI = QPSI+PRIORDIR(par(NPOS(NN+1:NN+NSI-1)),
	1                        ALPHA(K+1,1:NSI),NSI)
         K  = K + 1       
         NN = NN + NSI-1	    	  	   
	  ELSEIF(INFOS(9,J).EQ.2) THEN ! S~Markov
	   DO 99 I = 1,NSI
	    Ppar(NPARTH+K+1) = PRIORDIR(par(NPOS(NN+1:NN+NSI-1)),
	1                                psiprior(K+1,1:NSI),NSI)
  	    QPSI = QPSI+PRIORDIR(par(NPOS(NN+1:NN+NSI-1)),
	1                         ALPHA(K+1,1:NSI),NSI)
		K  = K + 1       
99	    NN = NN + NSI-1	    
	  ENDIF
100    CONTINUE

	 Fpar = PTHETA2(NPOS(1),nobs,d,ny,nz,nx,nu,ns,nt,IS,yk,
	1                par(1:nt),thetaprior(NPOS(1),:),
     2                tipo(NPOS(1)),pdll)
       Fpar = Fpar + SUM(Ppar(2:NPARTH+K)) ! log f(y|par,S)f(par,S)
	 
	 VQN(IG,1) = DEXP(QPSI)*QS*C*DEXP(-.5D0*MUC)/TRC       

200	 VHN(IG,1) = Fpar + DLOG(PS)

C ---------------------	
C Meng-Wong denominator 
C ---------------------	
	QS = ONE
	PS = ONE
	DO 400 IG = 1,G   
   	 IF (nv.GT.0) THEN
	  CALL DESIGNZ(nv,np(1),gibpar(IG,nt+1:nt+np(1)),INFOS,
	1               P1,P2,P3,P4,P5,P6)  
C PMAT(i,j) = Pr[Z(t+1)=i|Z(t)=j], Z = S1 x S2 x ... x Snv  
	  CALL PPROD(nv,nstot,INFOS,P1,P2,P3,P4,P5,P6,PMAT)
C ERGODIC solves PE: PE*(I-P') = 0
        CALL ERGODIC(nstot,PMAT,PE)

	  QS = PTR(1,gibZ(IG,1),1)
	  PS = PE(gibZ(IG,1))
	  CALL INT2SEQ(gibZ(IG,1),nv,INFOS,SEQ,IS(1,:))
	  DO 210 K = 2,nobs
	   QS = QS*PTR(K,gibZ(IG,K),gibZ(IG,K-1))
	   PS = PS*PMAT(gibZ(IG,K),gibZ(IG,K-1))
210	   CALL INT2SEQ(gibZ(IG,K),nv,INFOS,SEQ,IS(K,:)) 
       ENDIF

C PRIOR for THETA
       DO 310 I = 2,NPARTH
310	 Ppar(I) = PRIOR(gibpar(IG,NPOS(I)),thetaprior(NPOS(I),:),
     1                 tipo(NPOS(I)))  

C PRIOR for PSI and Q(PSI)~Dirichlet(a1,a2,...,aN)
	 QPSI = 0.D0
	 NN = NPARTH
	 K  = 0
	 DO 305 J = 1,nv
	  NSI = INFOS(8,J)  
	  IF(INFOS(9,J).EQ.1) THEN ! S~IID
	   Ppar(NPARTH+K+1) = PRIORDIR(gibpar(IG,NPOS(NN+1:NN+NSI-1)),
	1                               psiprior(K+1,1:NSI),NSI)
	   QPSI = QPSI+PRIORDIR(gibpar(IG,NPOS(NN+1:NN+NSI-1)),
	1                        ALPHA(K+1,1:NSI),NSI)
         K  = K + 1       
         NN = NN + NSI-1	    	  	   
	  ELSEIF(INFOS(9,J).EQ.2) THEN ! S~Markov
	   DO 304 I = 1,NSI
	    Ppar(NPARTH+K+1) = PRIORDIR(gibpar(IG,NPOS(NN+1:NN+NSI-1)),
	1                                psiprior(K+1,1:NSI),NSI)
	    QPSI = QPSI+PRIORDIR(gibpar(IG,NPOS(NN+1:NN+NSI-1)),
	1                         ALPHA(K+1,1:NSI),NSI)
		K  = K + 1       
304	    NN = NN + NSI-1	    
	  ENDIF
305    CONTINUE

	 Fpar = PTHETA2(NPOS(1),nobs,d,ny,nz,nx,nu,ns,nt,IS,yk,
	1                gibpar(IG,1:nt),thetaprior(NPOS(1),:),
     2                tipo(NPOS(1)),pdll)
       Fpar = Fpar + SUM(Ppar(2:NPARTH+K)) ! log f(y|par,S)f(par,S)

	 VHD(IG,1) = Fpar + DLOG(PS)

	 COM(:,1) = ZERO
	 DO 320 I = 1,NPARTH
320	 COM(I,1) = SUM((gibpar(IG,NPOS(1:NPARTH))-parm(1:NPARTH))
     #          * ISIGM(1:NPARTH,I))  	   
	 MUC = SUM(COM(1:NPARTH,1)*(gibpar(IG,NPOS(1:NPARTH))
     #	 - parm(1:NPARTH)))	   
	 
	 VQD(IG,1) = DEXP(QPSI)*QS*DEXP(-.5D0*MUC)*C/TRC
400   CONTINUE

	IND = MAXLOC(VHN(:,1))
	DET = VHN(IND(1),1) 

      MAT(:,1) = DEXP(VHN(:,1)-DET)/(DEXP(VHN(:,1)-MLSTART)+VQN(:,1))
      MAT(:,2) = VQD(:,1)/(DEXP(VHD(:,1)-MLSTART)+VQD(:,1))

      CALL NEWEYWESTCOV(G,2,1,MAT(:,1:2),SS)
      MLMW(2,1)   = SUM(MAT(:,1))/SUM(MAT(:,2))
      MLMW(1,1)   = DLOG(MLMW(2,1)) + DET
      MLMW(1:2,2) = SS(1,1)*G/SUM(MAT(:,1))**2 + 
     +              SS(2,2)*G/SUM(MAT(:,2))**2 + 
     +     	    - 2.D0*SS(1,2)*G/(SUM(MAT(:,1))*SUM(MAT(:,2)))

      MLMW(2,1) = MLMW(1,1) ! log scale	 
      DO 500 I=1,10
       MWNUM     = SUM(DEXP(VHN(:,1)-DET)
	1           / (DEXP(VHN(:,1)-MLMW(2,1))+VQN(:,1)))
	 MWDEN     = SUM(VQD(:,1)/(DEXP(VHD(:,1)-MLMW(2,1))+VQD(:,1)))
	 MLMW(2,1) = DLOG(MWNUM/MWDEN) + DET	! log-scale
500   CONTINUE

	DEALLOCATE (MAT,VHN,VHD,VQN,VQD)
	IF (nv.GT.0) DEALLOCATE (PTR,PMAT,PE,ALPHA,MOM)
	
	RETURN
	END
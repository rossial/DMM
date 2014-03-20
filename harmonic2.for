C --------------------------------------------------------------------------
C HARMONIC2 (no missing values) computes the harmonic mean estimates of 
C the Marginal Lilkelihood 
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C      
C 1 Modified HME (Geweke, 1999)
C 2 Modified and stabilized HME (Geweke, 1999)
C 1/ML = sum[f(S)f(THETA)/p(Y|THETA,S)P(S|THETA)p(THETA)], 
C        {S,THETA}~p(S,THETA|Y)
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
C --------------------------------------------------------------------------
	SUBROUTINE HARMONIC2(G,nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np,
	1                     INFOS,yk,gibpar,gibZ,thetaprior,psiprior,
     2                     tipo,pdll,MLH)

C INPUT
      INTEGER G,nobs,d(2),ny,nz,nx,nu,nv,ns(6),nstot,nt,np(3),
	1 INFOS(9,6),gibZ(G,nobs) 
	DOUBLE PRECISION yk(nobs,ny+nz),gibpar(G,nt+np(1)),
	1 thetaprior(nt,4),psiprior(np(2),np(3))
	CHARACTER*2 tipo(nt) 
	POINTER (pdll,fittizia)  !  ASSOCIATE  pointer P alla DLL ad una varibile fittizia

C OUTPUT
	DOUBLE PRECISION MLH(11,2)

C LOCALS
  	INTEGER NPAR,I,J,K,IG,NPOS(nt+np(1)),IFAIL,NQ,SEQ(nv),IS(nobs,6),
	1 NPVAL,IND1(1),NPARTH,NN,NSI,II,JJ
	DOUBLE PRECISION, ALLOCATABLE:: DEN2(:),ub(:) 
	INTEGER,ALLOCATABLE::IND(:,:),NIND(:)
	DOUBLE PRECISION parm(nt),SIGM(nt,nt),pval(11),
	1 COM(nt+1,nt),ISIGM(nt,nt),DET,CON(11)
	DOUBLE PRECISION,ALLOCATABLE:: PTR(:,:,:),PMAT(:,:),PE(:),
	1 ALPHA(:,:),MOM(:,:)
	DOUBLE PRECISION P1(INFOS(8,1),INFOS(8,1)),
	1 P2(INFOS(8,2),INFOS(8,2)),P3(INFOS(8,3),INFOS(8,3)),
     2 P4(INFOS(8,4),INFOS(8,4)),P5(INFOS(8,5),INFOS(8,5)),
     3 P6(INFOS(8,6),INFOS(8,6))
	DOUBLE PRECISION Ppar(nt+np(1)),Fpar,QTHETA,QPSI,PS,QS,QF,A0,SS(1,1)
	DOUBLE PRECISION ZERO,ONE,PI
	DATA ZERO/0.0D0/,ONE/1.0D0/,PI/3.141592653589793D0/	
C EXTERNAL FUNCTIONS
      DOUBLE PRECISION PTHETA2,PRIOR,PRIORDIR,CHI2INV !PPCHI2,G01FCF 
C EXTERNAL SUBROUTINES      
      EXTERNAL NEWEYWESTCOV2,DPOTRF,DPOTRI,DESIGNZ,PPROD,ERGODIC,INT2SEQ
      
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
	NPAR = NPARTH + np(1) 
	parm(:) = ZERO
	DO I = 1,NPARTH
	 parm(I) = SUM(gibpar(:,NPOS(I)))/DFLOAT(G)    
      ENDDO

	NQ = 0 
	CALL NEWEYWESTCOV2(G,NPARTH,NQ,gibpar(:,NPOS(1:NPARTH)),
	1       parm(1:NPARTH),SIGM(1:NPARTH,1:NPARTH)) ! THETA Var-covar

	COM(1:NPARTH,1:NPARTH) = SIGM(1:NPARTH,1:NPARTH)
	IFAIL = -1
C	CALL F01ADF(NPARTH,COM(1:NPARTH+1,1:NPARTH),NPARTH+1,IFAIL) ! Inverse var-covar
      CALL DPOTRF('L',NPARTH,COM(1:NPARTH,1:NPARTH),NPARTH,IFAIL) ! COM = L*L'
      DET = 1.D0 ! det(SIGM)
      DO I=1,NPARTH
	 DET = DET*COM(I,I)**2
      ENDDO
      CALL DPOTRI('L',NPARTH,COM(1:NPARTH,1:NPARTH),NPARTH,IFAIL) ! COM = VV^-1

	DO 30 I=1,NPARTH
	ISIGM(I,I) = COM(I,I)
	DO 30 J=1,I-1	
	ISIGM(I,J) = COM(I,J)
30	ISIGM(J,I) = ISIGM(I,J)

c	COM(1:NPARTH,1:NPARTH) = SIGM(1:NPARTH,1:NPARTH)
c     IFAIL = -1
c     CALL F03ABF(COM(1:NPARTH,1:NPARTH),NPARTH,NPARTH,DET,
c	1            WORK(1:NPARTH),IFAIL) ! det(SIGM)
	
C p = .05, .55, ...,1
      NPVAL = 11
      ALLOCATE (DEN2(G),IND(G,NPVAL),NIND(NPVAL),ub(NPVAL))
	pval(NPVAL) = 1.D0
	DO 40 I  = 1,NPVAL-1
	 pval(I) = 0.05D0 + DFLOAT(I-1)*.1D0 
40     ub(I)   = CHI2INV(I,NPARTH)  ! only tabulated values      
C40	 ub(I)   = PPCHI2(pval(I),DFLOAT(NPARTH),IFAIL) ! G01FCF(pval(I),DFLOAT(NPARTH),IFAIL)
	CON(:)   = (2.D0*PI)**(-NPARTH/2.D0)*DET**(-.5D0)/pval(:)    
	
	IF (nv.GT.0) THEN	  
	 ALLOCATE(PTR(nobs,nstot,nstot),PMAT(nstot,nstot),PE(nstot))
C Transition prob for QS
	 DO 55 I = 1,nstot-1
55	 PTR(1,I,1)     = SUM(ABS(gibZ(1:G,1).EQ.I))/DFLOAT(G)
       PTR(1,nstot,1) = ONE-SUM(PTR(1,1:nstot-1,1))
	
	 DO 57 K = 2,nobs
	 DO 57 I = 1,nstot-1
	 DO 57 J = 1,nstot
	  COM(1,1) = SUM(ABS(gibZ(1:G,K-1).EQ.J))
	  IF (COM(1,1).GT.ZERO) THEN
	   PTR(K,I,J) = SUM(ABS((gibZ(1:G,K).EQ.I).AND.(gibZ(1:G,K-1).EQ.J	   
     #                )))/COM(1,1)
	  ELSE
	   PTR(K,I,J) = ONE/DFLOAT(nstot)
	  ENDIF
57	  PTR(K,nstot,J) = ONE-SUM(PTR(K,1:nstot-1,J))

C Mean and VAR for PSI 
	 ALLOCATE (ALPHA(np(2),np(3)),MOM(np(1),2))
	 DO I=1,np(1)
	  MOM(I,1) = SUM(gibpar(:,nt+I))/DFLOAT(G)
	  MOM(I,2) = SUM(gibpar(:,nt+I)**2)/DFLOAT(G)
	  MOM(I,2) = MOM(I,2)-MOM(I,1)**2
	 ENDDO
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

	NIND(:) = 0
	QS = ZERO
	PS = ZERO
	IS(:,:) = 1
	Ppar(:) = 0.D0
	DO 110 IG = 1,G  
	 IF (nv.GT.0) THEN
	  CALL DESIGNZ(nv,np(1),gibpar(IG,nt+1:nt+np(1)),INFOS,
	1               P1,P2,P3,P4,P5,P6)  
C PMAT(i,j) = Pr[Z(t+1)=i|Z(t)=j], Z = S1 x S2 x ... x Snv  
	  CALL PPROD(nv,nstot,INFOS,P1,P2,P3,P4,P5,P6,PMAT)
C ERGODIC solves PE: PE*(I-P') = 0
        CALL ERGODIC(nstot,PMAT,PE)
	  PS = DLOG(PE(gibZ(IG,1)))      ! log P(S1)
        QS = DLOG(PTR(1,gibZ(IG,1),1)) ! log Q(S1)	  
	  CALL INT2SEQ(gibZ(IG,1),nv,INFOS,SEQ,IS(1,:))
	  DO 60 K = 2,nobs
   	   CALL INT2SEQ(gibZ(IG,K),nv,INFOS,SEQ,IS(K,:))
	   PS = PS + DLOG(PMAT(gibZ(IG,K),gibZ(IG,K-1)))
60	   QS = QS + DLOG(PTR(K,gibZ(IG,K),gibZ(IG,K-1)))
	 ENDIF

C Q(THETA)~Gaussian
	 QF = ZERO 
	 DO 70 I = 1,NPARTH
	 DO 70 J = 1,NPARTH	 				
70	 QF = QF + (gibpar(IG,NPOS(I))-parm(I))
     #    * (gibpar(IG,NPOS(J))-parm(J))*ISIGM(I,J)
	 QTHETA = -.5D0*QF 

C PRIOR for THETA
	 DO 80 I = 1,NPARTH
80	 Ppar(I) = PRIOR(gibpar(IG,NPOS(I)),thetaprior(NPOS(I),:),
     #                 tipo(NPOS(I)))  

C PRIOR for PSI and Q(PSI)~Dirichlet(a1,a2,...,aN)
C Mothod of Moments: a0 = m1(1-m1)/V1+1, ai = mi*a0, i=1,2,..,N
	 QPSI = ZERO	 
	 NN = NPARTH
	 K  = 0
	 DO 100 J = 1,nv
	  NSI = INFOS(8,J)  
	  IF(INFOS(9,J).EQ.1) THEN ! S~IID
	   Ppar(NPARTH+K+1) = PRIORDIR(gibpar(IG,NPOS(NN+1:NN+NSI-1)),
	1                               psiprior(K+1,1:NSI),NSI)
	   QPSI = QPSI+PRIORDIR(gibpar(IG,NPOS(NN+1:NN+NSI-1)),
	1                        ALPHA(K+1,1:NSI),NSI)
         K  = K + 1       
         NN = NN + NSI-1	    	  	            
	  ELSEIF(INFOS(9,J).EQ.2) THEN ! S~Markov
	   DO 90 I = 1,NSI
	    Ppar(NPARTH+K+1) = PRIORDIR(gibpar(IG,NPOS(NN+1:NN+NSI-1)),
	1                                psiprior(K+1,1:NSI),NSI)
	    QPSI = QPSI+PRIORDIR(gibpar(IG,NPOS(NN+1:NN+NSI-1)),
	1                         ALPHA(K+1,1:NSI),NSI)
		K  = K + 1       
90	    NN = NN + NSI-1	    
	  ENDIF
100    CONTINUE

	 Fpar = PTHETA2(NPOS(1),nobs,d,ny,nz,nx,nu,ns,nt,IS,yk,
	1                gibpar(IG,1:nt),thetaprior(NPOS(1),:),
     2                tipo(NPOS(1)),pdll)
	 Fpar = Fpar + SUM(Ppar(2:NPARTH+K))  ! log[f(y|par,S)f(par)]
	  
	 DEN2(IG) = QTHETA+QPSI+QS-Fpar-PS

	 DO 105 I=1,NPVAL-1
	 IF (QF.GT.ub(I)) THEN 
	  NIND(I)        = NIND(I) + 1  ! count where to put 0s
        IND(NIND(I),I) = IG           ! those to discard 
105    ENDIF
110   CONTINUE

      IND1 = MINLOC(DEN2)
	DET  = DEN2(IND1(1)) 
      DEN2(1:G) = DEXP(DEN2(1:G)-DET)

	MLH(NPVAL,1) = SUM(DEN2(1:G))/DFLOAT(G)
	CALL NEWEYWESTCOV2(G,1,1,DEN2(1:G),MLH(NPVAL,1),SS)
	MLH(NPVAL,2) = SS(1,1)/(DFLOAT(G)*MLH(NPVAL,1)**2)
	MLH(NPVAL,1) = MLH(NPVAL,1)*CON(NPVAL)
	DO 120 I=NPVAL-1,1,-1
	 DEN2(IND(1:NIND(I),I)) = 0.D0 	 
       MLH(I,1) = SUM(DEN2(1:G))/DFLOAT(G)
	 CALL NEWEYWESTCOV2(G,1,1,DEN2(1:G),MLH(I,1),SS)	 
	 MLH(I,2) = SS(1,1)/(DFLOAT(G)*MLH(I,1)**2)
120	 MLH(I,1) = MLH(I,1)*CON(I)      

	DO 130 I =1,NPVAL
130	IF (MLH(I,1).GT.ZERO) MLH(I,1) = -DLOG(MLH(I,1))-DET
	
	DEALLOCATE(DEN2,IND,NIND,ub)
	IF (nv.GT.0) DEALLOCATE (PTR,PMAT,PE,ALPHA,MOM)

 	RETURN
	END
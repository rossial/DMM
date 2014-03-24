C ------------------------------------------------------------
C OPGH coumputes the Var-Covar matrix of Parameters inverting 
C the Observed Information matrix calcuted by outer product
C of gradient estimator - e.g. Davidson MacKinnon pp 265-66
C for MS-VAR(1) models      
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
      SUBROUTINE OPGH(nobs,ny,nz,nx,nu,nt,nv,ns,nstot,np,pdll,yk,IYK,
	1                INFOS,theta,psi,thetaprior,HESS,thetase,psise,
     1                SSMOOTH,INN,IFAIL)
      
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
	POINTER (pdesign,DESIGN) ! IMPORTANT associo il puntatore pdesign alla Interface definita

! Input 
	INTEGER nobs,ny,nz,nx,nu,nt,nv,ns(6),nstot,np,pdll,IYK(nobs,ny+1),
     1 INFOS(9,6)
      DOUBLE PRECISION yk(nobs,ny+nz),theta(nt),psi(np),
     1 thetaprior(nt,4),HESS((nt+np)*(nt+np+1)/2)
! Output
      INTEGER IFAIL
      DOUBLE PRECISION thetase(nt),psise(np),SSMOOTH(nobs,nstot),
     1 INN(nobs,ny)
! Locals
      INTEGER I,J,K,IFREE(nt+np),NFREE,NFT
	DOUBLE PRECISION DRI,PAR(nt+np),DLL(nobs),DLLM(nobs),RMAX		
	DOUBLE PRECISION c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6))	
	DOUBLE PRECISION, ALLOCATABLE:: GRAD(:,:),P(:),LTR(:),
	1 W(:),OP(:,:),VC(:,:),pro(:,:)

	thetase(:) = 0.D0
      psise(:)   = 0.D0
	IFAIL = 0	
      DRI   = 1.D-3
	NFREE = 0
	DO 20 I = 1,nt
	 IF ((theta(I).GT.thetaprior(I,3)).AND.
	1     (theta(I).LT.thetaprior(I,4))) THEN
         NFREE        = NFREE + 1 
	   IFREE(NFREE) = I   
         PAR(NFREE)   = theta(I)
20     ENDIF
      NFT = NFREE ! only theta free param
      DO 21 I = 1,np
	 IF ((psi(I).GT..001D0).AND.(psi(I).LT..999D0)) THEN
         NFREE        = NFREE + 1 
	   IFREE(NFREE) = I+nt
         PAR(NFREE)   = psi(I)
21     ENDIF
      
	ALLOCATE (LTR(NFREE*(NFREE+1)/2),W(NFREE))
      DO 25 I=1,NFREE
	DO 25 J=1,I
      K = IFREE(I)*(IFREE(I)+1)/2-IFREE(I)+IFREE(J)              
25    LTR(I*(I+1)/2-I+J) = HESS(K)
      IFAIL = 0	
      CALL SYMINV(LTR,NFREE,LTR,W,J,IFAIL,RMAX)
      
      pdesign = getprocaddress(pdll, "design_"C)
     	CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
      CALL HF(nobs,nx,nstot,nz,nu,ns,nv,np,psi,1,yk,IYK,INFOS,
     1        c,a,F,R,SSMOOTH,INN,DLL)

      IF (IFAIL.GT.0) THEN
       ALLOCATE (GRAD(Nobs,NFREE),P(NFREE),OP(NFREE,NFREE),
     1  VC(NFREE,NFREE),pro(nobs*nx,NFREE))

       LTR(:)   = 0.D0
	 W(:)     = 0.D0
	 GRAD(:,:)= 0.D0
	 P(:)     = 0.D0
       OP(:,:)  = 0.D0
       P(1:NFREE) = PAR(1:NFREE)*DRI  ! delta 
	 DO I=1,NFT
        IF(DABS(P(I)).LT.1.D-13) P(I)=1.D-13
        IF (((theta(IFREE(I))+P(I)).GT.thetaprior(IFREE(I),4)).OR. 
	1      ((theta(IFREE(I))+P(I)).LT.thetaprior(IFREE(I),3))) THEN 
         IF (theta(IFREE(I)).GT.0.D0) THEN 
          P(I) = thetaprior(IFREE(I),3) - theta(IFREE(I)) 
         ELSE 
          P(I) = thetaprior(IFREE(I),4) - theta(IFREE(I))
         ENDIF
        ENDIF 
       END DO
      
       DO I=NFT+1,NFREE
        IF (DABS(P(I)).LT.1.D-13) P(I)=1.D-13
        IF (((psi(I-NFT)+P(I)).GT.999D0).OR.
	1     ((psi(I-NFT)+P(I)).LT..001D0)) THEN 
         P(I) = .001D0 - psi(I-NFT)                 
        ENDIF
       END DO
   
C -----------
C Main Cycle
C -----------     
	 DO 1000 I=1,NFREE
	  PAR(I) = PAR(I) + P(I)
        IF (I.LE.NFT) THEN
         theta(IFREE(I)) = PAR(I)
        ELSE
         psi(I-NFT) = PAR(I)
        ENDIF
	  CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
        CALL HF(nobs,nx,nstot,nz,nu,ns,nv,np,psi,0,yk,IYK,INFOS,
     1        c,a,F,R,SSMOOTH,INN,DLLM)
	  PAR(I) = PAR(I) - P(I)
	  IF (I.LE.NFT) THEN
         theta(IFREE(I)) = PAR(I)
        ELSE
         psi(I-NFT) = PAR(I)
        ENDIF
        GRAD(2:nobs,I) = (DLLM(2:nobs) - DLL(2:nobs))/P(I)
1000   CONTINUE	
      
C Use OPG if Hessian is bad      
        DO 300 I=1,NFREE 
	  DO 300 J=1,I
300	  OP(I,J) = SUM(GRAD(2:Nobs,I)*GRAD(2:Nobs,J))
	
	  DO 150 I=1,NFREE
	  DO 150 J=1,I 
150     LTR(i*(i+1)/2-i+j) = OP(I,J)
	  CALL SYMINV(LTR,NFREE,LTR,W,J,IFAIL,RMAX)
        IF (IFAIL.NE.0) GO TO 1111
        DEALLOCATE (OP,W,P,GRAD,VC,pro)    	 
      ENDIF
       
      DO i=1,NFT
	  thetase(IFREE(I))=dsqrt(LTR(i*(i+1)/2))
      ENDDO
      DO i=NFT+1,NFREE
	  psise(I-NFT)=dsqrt(LTR(i*(i+1)/2))
      ENDDO

1111  DEALLOCATE (LTR)  
	RETURN 
	END

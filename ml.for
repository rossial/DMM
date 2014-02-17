C --------------------------------------------------------------------
C ML estimates model parameters theta by maximum likelihood
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
C --------------------------------------------------------------------	
	SUBROUTINE ML(nobs,d,ny,nz,nx,nu,nt,nv,ns,np,INFOS,pdll,INDT,yk,IYK,
	1              S,thetaprior,theta,psi,IMSVAR,HESS,FOPT)

      USE dfwin
	INTERFACE
	 SUBROUTINE DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R) 
	 INTEGER ny,nz,nx,nu,ns(6),nt
	 DOUBLE PRECISION theta(nt)
	 DOUBLE PRECISION c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6))
	 END SUBROUTINE
      END INTERFACE
C     INTEGER POINTER pdll
	POINTER (pdll,fittizia)  
	POINTER (pdesign,DESIGN)

C INPUT
	INTEGER nobs,d(2),ny,nz,nx,nu,nt,nv,ns(6),np,INFOS(9,6),
     1 INDT(nt+2),IYK(nobs,ny+1),S(nobs,6)
	DOUBLE PRECISION thetaprior(nt,4),yk(nobs,ny+nz)

C OUTPUT
	INTEGER IMSVAR
      DOUBLE PRECISION FOPT,theta(nt),psi(max(1,np)),
     1 HESS((nt+np)*(nt+np+1)/2)      

C LOCALS
	INTEGER NPAR,NTHETA,I,J,K,IFAIL,IFAILNEW,IOPT,NCLIN,NCNLN,LDA,ITER,
     1 LDCJ,LWORK,LIWORK
	DOUBLE PRECISION FOPTNEW
	INTEGER, ALLOCATABLE:: ISTATE(:),IWORK(:)
      INTEGER*8, ALLOCATABLE:: IU(:)
	DOUBLE PRECISION, ALLOCATABLE:: CHINEW(:),CHI(:),LBV(:),UBV(:),
     1 CLAMDA(:),OBJGRD(:),WORK(:),CC(:),U(:),CJAC(:,:),RR(:,:),AA(:,:)
      DOUBLE PRECISION, ALLOCATABLE:: c(:,:,:),H(:,:,:),G(:,:,:),
     2 a(:,:),F(:,:,:),R(:,:,:)
	EXTERNAL FUNCT1,E04UCF,E04UEF,CONFUN
	
	NTHETA = INDT(nt+2)! # of free parameters theta
      NPAR   = NTHETA+np ! # of parameters (total)
            
C Check if the model is an MS-VAR(1):      
C d(1)=d(2)=0, H=I(nx), G=0      
      IMSVAR=0
      ALLOCATE(c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6)))  
      pdesign = getprocaddress(pdll, "design_"C)
	CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
C Check H (ny x nx x ns(2)) and G (ny x nu x ns(3))      
      IF ((nx.EQ.ny).AND.(SUM(d).EQ.0)) THEN
       IMSVAR = 1   
	 DO I = 1,ny
	   IF (H(I,I,1).NE.1.D0) THEN
		IMSVAR = 0
		GO TO 123
         ENDIF 
	   DO J = 1,I-1
		IF (H(I,J,1).NE.0.D0) THEN
		  IMSVAR = 0
		  GO TO 123
		ENDIF
         ENDDO
         DO J = I+1,ny
		IF (H(I,J,1).NE.0.D0) THEN
		  IMSVAR = 0
		  GO TO 123
		ENDIF
         ENDDO
	  IF (SUM(DABS(G(I,1:nu,1))).NE.0.D0) THEN
          IMSVAR = 0
		GO TO 123  
        ENDIF    
       ENDDO            
      ENDIF        
	
123   DEALLOCATE(c,H,G,a,F,R)       
C Shrink theta, LB and UB
	ALLOCATE(CHINEW(NPAR),CHI(NPAR))

C Starting values - set via the 1st entry of th prior when possible 	
      DO I = 1,NTHETA
	 IF ((thetaprior(INDT(I),1).LT.thetaprior(INDT(I),3)).OR.
     #     (thetaprior(INDT(I),1).GT.thetaprior(INDT(I),4))) THEN
        CHINEW(I) = theta(INDT(I))   
       ELSE
        CHINEW(I) = thetaprior(INDT(I),1)
	 ENDIF
      ENDDO
C Starting values for PSI      
      CHINEW(NTHETA+1:NPAR) = psi(1:np)  
      
C Linear constraints for transition probabilities
C LBV <= AA*theta <= UBV
      NCLIN = 0
      DO I = 1,nv
        IF(INFOS(9,I).EQ.1) THEN      ! Independent
            NCLIN = NCLIN+1
        ELSEIF(INFOS(9,I).EQ.2) THEN  ! Markov  
            NCLIN = NCLIN+INFOS(8,I)
        ENDIF    
      ENDDO 

C Lower and upper bounds      
      ALLOCATE(LBV(NPAR+NCLIN),UBV(NPAR+NCLIN),AA(NCLIN,NPAR))
	LBV(1:NTHETA) = thetaprior(INDT(1:NTHETA),3)
	UBV(1:NTHETA) = thetaprior(INDT(1:NTHETA),4)
      LBV(NTHETA+1:NPAR)     = 1.D-3
	LBV(NPAR+1:NPAR+NCLIN) = 1.D-6
      UBV(NTHETA+1:NPAR)     = 1.D0-1.D-3
      UBV(NPAR+1:NPAR+NCLIN) = 1.D0-1.D-6
      K = NTHETA
      J = 1
      AA(:,:) = 0.D0
      DO I = 1,nv
        IF(INFOS(9,I).EQ.1) THEN      ! Independent
          AA(J,K+1:K+INFOS(8,I)-1) = 1.D0
          K = K+INFOS(8,I)-1
          J = J+1
        ELSEIF(INFOS(9,I).EQ.2) THEN  ! Markov  
          DO ITER = 1,INFOS(8,I)
            AA(J,K+1:K+INFOS(8,I)-1) = 1.D0
            K = K+INFOS(8,I)-1
            J = J+1
          ENDDO
        ENDIF    
      ENDDO 
	
C --------------------------------------------------------
C Set E04UCF parameters: IU all integers, U data + bounds
C --------------------------------------------------------
	ALLOCATE (IU(72),U(nobs*(2*ny+nz+7)+3*nt+2))	
	IU(1)   = nobs
	IU(2:3) = d(1:2)
	IU(4)   = ny
	IU(5)   = nz
	IU(6)   = nx
	IU(7)   = nu
	IU(8)   = nt
      IU(9:14)= ns(1:6)
	IU(15)  = pdll
      IU(16)  = nv
      DO J=1,6
	 IU(17+9*(J-1):16+J*9) = INFOS(1:9,J)
      ENDDO
      IU(71) = np
      IU(72) = IMSVAR
      
	DO J=1,ny+nz
	 U(1+nobs*(J-1):J*nobs) = yk(:,J)
	ENDDO
	U(nobs*(ny+nz)+1:nobs*(ny+nz)+nt)      = thetaprior(1:nt,3)
	U(nobs*(ny+nz)+nt+1:nobs*(ny+nz)+2*nt) = thetaprior(1:nt,4)
	I = nobs*(ny+nz)+2*nt+1
	U(I:I+nt+1) = INDT(1:nt+2)
	I = I+nt+2
	DO J=1,ny+1
	 U(I+nobs*(J-1):I-1+nobs*J) = IYK(1:nobs,J)
	ENDDO
	I = I+nobs*(ny+1)
	DO J = 1,6
	 U(I+(J-1)*nobs:I-1+J*nobs) = S(1:nobs,J)
	ENDDO
	
C -----------------------------------------------    
C  Likelihood Maximization via E04UCF (NAG mk.17)
C -----------------------------------------------
 	CALL E04UEF('Derivative level = 0')
	CALL E04UEF('Hessian = Yes')
	CALL E04UEF('Major iteration limit = 400')
	CALL E04UEF('Minor iteration limit = 300')
	CALL E04UEF('Cold start')
	CALL E04UEF('Major print level = 1')
      IOPT  = 0
	FOPT  = 1D100
      NCNLN = 0
      LDA   = max(1,NCLIN)
	LDCJ  = max(1,NCNLN)
	LWORK = 2*NPAR**2+20*NPAR+11*NCLIN      
	LIWORK=3*NPAR+2*NCNLN+NCLIN
	ALLOCATE(ISTATE(NPAR+NCNLN+NCLIN),IWORK(LIWORK),
     1 CLAMDA(NPAR+NCNLN+NCLIN),OBJGRD(NPAR),WORK(LWORK),CC(LDCJ),
     1 CJAC(LDCJ,NPAR),RR(NPAR,NPAR))

1234	IFAILNEW = -1
      CALL E04UCF(NPAR,NCLIN,NCNLN,LDA,LDCJ,NPAR,AA,LBV,UBV,CONFUN,
     1 FUNCT1,ITER,ISTATE,CC,CJAC,CLAMDA,FOPTNEW,OBJGRD,RR,CHINEW,
     2 IWORK,LIWORK,WORK,LWORK,IU,U,IFAILNEW)

C Hessian matrix for theta and psi
      IF (IOPT.EQ.0) THEN
       HESS(:) = 0.D0
       DO 50 I=1,NTHETA      
       DO 50 J=1,I
       K = INDT(I)*(INDT(I)+1)/2-INDT(I)+INDT(J)    
50     HESS(K) = SUM(RR(:,I)*RR(:,J))

       DO 55 I=NTHETA+1,NTHETA+np      
       DO 55 J=1,NTHETA
       K = (I+nt-NTHETA)*(I+nt-NTHETA+1)/2 - (I+nt-NTHETA) + INDT(J)    
55     HESS(K) = SUM(RR(:,I)*RR(:,J))

       DO 60 I=NTHETA+1,NTHETA+np      
       DO 60 J=NTHETA+1,I
       K = (I+nt-NTHETA)*(I+nt-NTHETA+1)/2 - I + J    
60     HESS(K) = SUM(RR(:,I)*RR(:,J))
      ENDIF 

      IF (FOPTNEW.LT.FOPT) THEN
	 CHI(:) = CHINEW(:)
	 FOPT   = FOPTNEW
	 IFAIL  = IFAILNEW
	ENDIF

	IF (((IFAILNEW.EQ.1).OR.(IFAILNEW.EQ.4).OR.(IFAILNEW.EQ.6))
     * 	.AND.(IOPT.LT.2)) THEN
	 
	 IOPT = IOPT + 1
	 CALL E04UEF('Warm start')
       RR(:,:) = 0.D0
	 DO 30 I=1,NPAR
30	 RR(I,I)=1.D0
	 GO TO 1234  	
	ENDIF
	theta(1:nt)           = thetaprior(1:nt,3)
	theta(INDT(1:NTHETA)) = CHI(1:NTHETA)
      IF (NPAR.GT.NTHETA) psi(1:np) = CHI(NTHETA+1:NPAR)
      
	DEALLOCATE (ISTATE,IWORK,CLAMDA,OBJGRD,WORK,CC,CJAC,RR,AA)
      DEALLOCATE (CHINEW,CHI,LBV,UBV,IU,U)
	RETURN
	END
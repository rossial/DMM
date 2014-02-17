C -------------------------------------------------------------------------
C FUNCT1 computes -loglikelihood and the numerical derivatives for E04UCF      
C Developed by A.Rossi, C.Planas and G.Fiorentini      
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
      SUBROUTINE FUNCT1(MODE,NPAR,CHI,DLL,OBJGRD,NNN,IU,U)
    
	USE dfwin
	INTERFACE
	 SUBROUTINE DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R) 
	 INTEGER ny,nz,nx,nu,ns(6),nt
	 DOUBLE PRECISION theta(nt)
	 DOUBLE PRECISION c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6))
	 END SUBROUTINE
      END INTERFACE
      INTEGER POINTER pdll
	POINTER (pdll,fittizia)  ! ASSOCIATE  pointer pdll alla DLL ad una varibile fittizia
	POINTER (pdesign,DESIGN)

C INPUT
      INTEGER MODE,NPAR,NNN
      INTEGER*8 IU(72)
      DOUBLE PRECISION U(IU(1)*(2*IU(4)+IU(5)+7)+3*IU(8)+2),
	1 CHI(NPAR),OBJGRD(NPAR)

C OUTPUT 
      DOUBLE PRECISION DLL

C LOCALS
	INTEGER nobs,d(2),ny,nz,nx,nu,nv,ns(6),nt,np,nmis,I,J,nstot,
     1 INFOS(9,6),IMSVAR
      INTEGER, ALLOCATABLE:: INDT(:),IYK(:,:),S(:,:)
      DOUBLE PRECISION, ALLOCATABLE:: theta(:),thetaprior(:,:),psi(:),
	1 yk(:,:)
      DOUBLE PRECISION, ALLOCATABLE:: c(:,:,:),H(:,:,:),G(:,:,:),
     2 a(:,:),F(:,:,:),R(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE:: LIKE(:),XT(:,:),PT(:,:,:),	
	1 Xdd(:,:),Pdd(:,:,:)
      DOUBLE PRECISION, ALLOCATABLE:: XSMOOTH(:,:),XSSE(:,:),
     1 SSMOOTH(:,:),INN(:,:)
      CHARACTER*200 NMLNAME,PATH

C Retrive metainformation
      nobs   = IU(1)
      d(1:2) = IU(2:3) 
      ny     = IU(4)    
      nz     = IU(5)    
      nx     = IU(6)    
      nu     = IU(7)    
      nt     = IU(8)    
      ns(1:6)= IU(9:14) 
      pdll   = IU(15)   
      nv     = IU(16)
      np     = IU(71)
      IMSVAR = IU(72)
      
      ALLOCATE(INDT(nt+2),IYK(nobs,ny+1),S(nobs,6))
      ALLOCATE(theta(nt),thetaprior(nt,4),psi(max(1,np)),yk(nobs,ny+nz))
      ALLOCATE(LIKE(nobs),XT(0:nobs,nx),PT(0:nobs,nx,nx),	
	1 Xdd(max(d(1),1),nx),Pdd(max(d(1),1),nx,nx))	

      DO J=1,6
	 INFOS(1:9,J) = IU(17+9*(J-1):16+J*9)
      ENDDO      
	DO J=1,ny+nz
	 yk(:,J) = U(1+nobs*(J-1):J*nobs)
	ENDDO
	thetaprior(1:nt,3) = U(nobs*(ny+nz)+1:nobs*(ny+nz)+nt)
	thetaprior(1:nt,4) = U(nobs*(ny+nz)+nt+1:nobs*(ny+nz)+2*nt)
	I = nobs*(ny+nz)+2*nt+1
	INDT(1:nt+2) = U(I:I+nt+1)
	I = I+nt+2
	DO J=1,ny+1
	 IYK(1:nobs,J) = U(I+nobs*(J-1):I+nobs*J-1)
	ENDDO
	I = I+nobs*(ny+1)
	DO J = 1,6
	 S(1:nobs,J) = U(I+(J-1)*nobs:I-1+J*nobs)
	ENDDO

C Expand theta and psi
      theta(1:nt)            = thetaprior(1:nt,3) 
      theta(INDT(1:NPAR-np)) = CHI(1:NPAR-np)
      IF (np.GT.0) psi(1:np) = CHI(NPAR-np+1:NPAR)

C Evaluate the likelihood
      ALLOCATE(c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6)))  
      pdesign = getprocaddress(pdll, "design_"C)
	CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
	nmis = ny*nobs-SUM(IYK(1:nobs,ny+1))	
	IF (nv.EQ.0) THEN
       IF (nmis.GT.0) THEN
	  CALL IKF(d,ny,nz,nx,nu,ns,S(1:max(d(1),1),1:6),
	1           yk(1:max(d(1),1),1:ny+nz),IYK(1:max(d(1),1),1:ny+1),
     2           c,H,G,a,F,R,Xdd,Pdd,LIKE(1:max(d(1),1)))
	  XT(d(1),1:nx)      = Xdd(max(d(1),1),1:nx)   
	  PT(d(1),1:nx,1:nx) = Pdd(max(d(1),1),1:nx,1:nx)
	  CALL KF(nobs,d,ny,nz,nx,nu,ns,S,yk,IYK,c,H,G,a,F,R,XT,PT,LIKE)
       ELSE       
	  CALL IKF2(d,ny,nz,nx,nu,ns,S(1:max(d(1),1),1:6),
	1            yk(1:max(d(1),1),1:ny+nz),
     2            c,H,G,a,F,R,Xdd,Pdd,LIKE(1:max(d(1),1)))
	  XT(d(1),1:nx) = Xdd(max(d(1),1),1:nx)   
	  PT(d(1),1:nx,1:nx) = Pdd(max(d(1),1),1:nx,1:nx)
	  CALL KF2(nobs,d,ny,nz,nx,nu,ns,S,yk,c,H,G,a,F,R,XT,PT,LIKE)
       ENDIF
      ELSE  
        nstot = PRODUCT(INFOS(8,1:nv))
        IF (IMSVAR.EQ.1) THEN   ! Hamilton (Ecoca 1989) filter
         ALLOCATE(SSMOOTH(nobs,nstot),INN(nobs,ny))   
         CALL HF(nobs,nx,nstot,nz,nu,ns,nv,np,psi,0,yk,IYK,INFOS,
     1           c,a,F,R,SSMOOTH,INN,LIKE)
         DEALLOCATE(SSMOOTH,INN) 
        ELSE   ! KIM (JoE 1994) algorithm
         ALLOCATE(XSMOOTH(nobs,nx),XSSE(nobs,nx),SSMOOTH(nobs,nstot),
     1            INN(nobs,ny))
         CALL KIM(nobs,d,ny,nz,nx,nu,ns,nstot,nv,np,INFOS,yk,IYK,
     1            c,H,G,a,F,R,psi,0,XSMOOTH,XSSE,SSMOOTH,INN,LIKE)
         DEALLOCATE(XSMOOTH,XSSE,SSMOOTH,INN) 
        ENDIF
        
      ENDIF      
	
      DLL = -SUM(LIKE(d(1)+1:nobs)) 
      
      DEALLOCATE(LIKE)
      DEALLOCATE(c,H,G,a,F,R)       
      DEALLOCATE(INDT,IYK,S,theta,thetaprior,psi,yk,XT,PT,Xdd,Pdd)

      RETURN
      END

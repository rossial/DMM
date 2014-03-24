C -------------------------------------------------------------
C KIM implements the filtering and smoothing algorithm of 
C Kim, Journal of Econometrics, 1994
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

C OUTPUT:
C   SFILT   = Pr(S(t)|y^t), (nobs x nk)
C   SSMOOTH = Pr(S(t)|y^T), (nobs x nk)
C   XFILT   = E(x(t)|y^t),  (nobs x nx)
C   XSMOOTH = E(x(t)|y^T),  (nobs x nx)
C   INN     = E(y(t)|y^t-1),(nobs x ny
C
C INTERMEDIATE:
C   SP1  = Pr(S(t-1)=i,S(t)=j|y^t)      nk x nk
C   SP0  = Pr(S(t)=i,S(t+1)=j|y^(t+1))  nk x nk
C   X1   = x(t|t-1,i,j)             nobs x nx x nk x nk
C   P1   = E[(x(t)-X1)^2|y^(t-1)]   nobs x nx x nx x nk x nk
C   X0   = x(t|t,i,j)               nx x nk x nk
C   P0   = E[(x(t)-X0)^2|y^t]       nx x nx x nk x nk
C   XI   = x(t|t,i)                 nobs x nx x nk
C   PI   = E[(x(t)-XI)^2|y^t]       nobs x nx x nx x nk
C   V    = Var(INN)                 ny x ny x nk x nk
C
C   SPT1 = Pr(S(t)=j,S(t+1)=k|y^T)              nk x nk
C   XS   = E(x(t)|S(t)=j,S(t+1)=k,y^T)          nx x nk x nk
C   PS   = E[(x(t)-XS)^2|S(t)=j,S(t+1)=k,y^T]   nx x nx x nk x nk
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
C -------------------------------------------------------------------------
      SUBROUTINE KIM(nobs,d,ny,nz,nx,nu,ns,nk,nv,np,INFOS,yk,IYK,
     1 c,H,G,a,F,R,psi,ismoother,XSMOOTH,XSSE,SSMOOTH,INN,LIKE)  
C INPUT
	INTEGER nobs,ny,nz,nx,nu,ns(6),nk,nv,np,ismoother,
     1 d(2),IYK(nobs,ny+1),INFOS(9,6) 
	DOUBLE PRECISION yk(nobs,ny+nz),c(ny,max(nz,1),ns(1)),
	1 H(ny,nx,ns(2)),G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),
     2 R(nx,nu,ns(6)),psi(max(1,np))
C OUTPUT	
	DOUBLE PRECISION XSMOOTH(nobs,nx),XSSE(nobs,nx),SSMOOTH(nobs,nk),
     1 LIKE(nobs),INN(nobs,ny)
      
C LOCALS
	INTEGER imain,iny,I,J,K,L,IFAIL
	INTEGER IPIV(nx),IMAX(1),IS(max(1,d(1)),6),SEQ(1)
	DOUBLE PRECISION DETV,RSS,lfy,maxfyss,fsum
      DOUBLE PRECISION, ALLOCATABLE:: SP1(:,:),SP0(:,:),X0(:,:,:),
	1 P0(:,:,:,:),X1(:,:,:,:),P1(:,:,:,:,:),fyss(:,:),
     1 XI(:,:,:),PI(:,:,:,:),Pdd(:,:,:),Xdd(:,:),INNIJ(:),
     1 V(:,:),AUX(:,:),AUX1(:,:,:,:),XFILT(:,:),SFILT(:,:),
     1 SPT1(:,:),XS(:,:,:),PS(:,:,:,:),XSC(:,:),PSC(:,:,:),PTIL(:,:),
     1 FP(:,:),HP1(:,:),RG(:,:),COM(:,:),VINV(:,:),HPV(:,:),W(:),
     1 PP1(:,:),PP2(:,:),PP3(:,:),PP4(:,:),PP5(:,:),PP6(:,:),P(:,:)

C EXTERNAL SUBROUTINES
      EXTERNAL DESIGNZ,PPROD,ERGODIC,INT2SEQ,IKF2,KF2,DPOTRF,DPOTRI,
     1 DGETRF,DGETRS
      
      LIKE(:)  = 0.D0
      INN(:,:) = 0.D0
      ALLOCATE(PP1(INFOS(8,1),INFOS(8,1)),PP2(INFOS(8,2),INFOS(8,2)),
     1         PP3(INFOS(8,3),INFOS(8,3)),PP4(INFOS(8,4),INFOS(8,4)),
     1         PP5(INFOS(8,5),INFOS(8,5)),PP6(INFOS(8,6),INFOS(8,6)),
     1         P(nk,nk))
      CALL DESIGNZ(nv,np,psi,INFOS,PP1,PP2,PP3,PP4,PP5,PP6)  
      CALL PPROD(nv,nk,INFOS,PP1,PP2,PP3,PP4,PP5,PP6,P)          
      DEALLOCATE(PP1,PP2,PP3,PP4,PP5,PP6)
      
C S-filter initialization Pr[S(d-1),S(d)|y^d]
C X-filter initialization      
     	ALLOCATE(SP1(nk,nk),SP0(nk,nk),X0(nx,nk,nk),P0(nx,nx,nk,nk),
	1 X1(nobs,nx,nk,nk),P1(nobs,nx,nx,nk,nk),INNIJ(ny),
     1 XI(nobs,nx,nk),PI(nobs,nx,nx,nk),fyss(nk,nk),V(ny,ny),
     1 AUX(nx,nx),AUX1(nx,nx,nk,nk),XFILT(nobs,nx),SFILT(nobs,nk),
     1 Pdd(MAX(d(1),1),nx,nx),Xdd(MAX(d(1),1),nx),FP(nx,nx),HP1(ny,nx),
     1 RG(nx,ny),COM(ny+1,ny),VINV(ny,ny),HPV(nx,ny))	

      CALL ERGODIC(nk,P,fyss(1:nk,1)) 
      IF (d(1).LE.1) THEN          
        RSS = 0.D0  
        DO J = 1,nk   ! S(d)
	   DO I = 1,nk  ! S(d-1)                     
          CALL INT2SEQ(J,nv,INFOS,SEQ,IS(max(1,d(1)),:))             
          CALL IKF2(d,ny,nz,nx,nu,ns,IS(1,:),yk(1,:),
	1              c,H,G,a,F,R,Xdd(1,:),Pdd(1,:,:),LIKE(1))  
          
          IF (d(1).EQ.0) THEN
           X1(1,:,1,1)   = Xdd(1,:)   ! 0|0
           P1(1,:,:,1,1) = Pdd(1,:,:)
           CALL KF2(1,d,ny,nz,nx,nu,ns,IS,yk(1,:),c,H,G,a,F,R,
	1              X1(1:2,:,1,1),P1(1:2,:,:,1,1),LIKE(1))
           SP0(I,J) = P(J,I)*fyss(I,1)*LIKE(1)            
         
           X0(:,I,J)   = X1(2,:,1,1)
	     P0(:,:,I,J) = P1(2,:,:,1,1)       
           XI(1,:,J)   = X1(2,:,1,1)
	     PI(1,:,:,J) = P1(2,:,:,1,1)       
          ELSE
           SP0(I,J) = P(J,I)*fyss(I,1)
           X0(:,I,J)   = Xdd(1,:)   ! 1|1
	     P0(:,:,I,J) = Pdd(1,:,:)       
           XI(1,:,J)   = Xdd(1,:) 
	     PI(1,:,:,J) = Pdd(1,:,:)
          ENDIF
          RSS = RSS+SP0(I,J)
         ENDDO
         
        ENDDO        
      ELSE
        WRITE(*,*) 'WARNING: d(1)=2 not implemeted yet'
        PAUSE
        STOP
      ENDIF    
      SP0(:,:) = SP0(:,:)/RSS           
      DO I =1,nk
       SFILT(max(1,d(1)),I) = SUM(SP0(:,I))
      ENDDO
      DEALLOCATE(Xdd,Pdd)

C ---------
C Filtering
C ---------
	XFILT(:,:) = 0.D0	
	DO 1000 imain=max(2,d(1)+1),nobs 
       iny = IYK(imain,ny+1)   
	 DO I = 1,nk
	 DO J = 1,nk
        CALL INT2SEQ(J,nv,INFOS,SEQ,IS(1,:))    
C ------------------------------------
C Prediction       x1 = a+F*x0 
C Prediction var.  P1 = F*p0*F'+ R*R'
C ------------------------------------
	  DO 10 K=1,nx
10      X1(imain,K,I,J)=a(K,IS(1,4))+SUM(F(K,:,IS(1,5))*XI(imain-1,:,I))

	  DO 20 K=1,nx
        DO 20 L=1,nx
20      FP(K,L) = SUM(F(K,:,IS(1,5))*PI(imain-1,:,L,I))

        DO 30 K=1,nx
        P1(imain,K,K,I,J) = SUM(FP(K,:)*F(K,:,IS(1,5))) 
     +          + SUM(R(K,1:nu,IS(1,6))*R(K,1:nu,IS(1,6)))
	  DO 30 L=1,K-1	
        P1(imain,K,L,I,J) = SUM(FP(K,:)*F(L,:,IS(1,5))) 
     +          + SUM(R(K,1:nu,IS(1,6))*R(L,1:nu,IS(1,6)))
30      P1(imain,L,K,I,J) = P1(imain,K,L,I,J)

C -------------------------------
C Innovations: INN = yk-H*X1-c*z
C -------------------------------
	  DO 40 K=1,iny
40      INNIJ(K) = yk(imain,IYK(imain,K))
     +         - SUM(H(IYK(imain,K),1:nx,IS(1,2))*X1(imain,:,I,J))
     +         - SUM(c(IYK(imain,K),1:nz,IS(1,1))*yk(imain,ny+1:ny+nz)) 
        INN(imain,1:ny) = INN(imain,1:ny)
     +                  + INNIJ(1:ny)*P(J,I)*SFILT(imain-1,I)
        
C ---------------------------------------------------------
C Innovation variance V = H*P1*H' + G*G' + H*R*G' + G*R'*H'
C ---------------------------------------------------------
        DO 50 K=1,iny
        DO 50 L=1,nx
50      HP1(K,L)=SUM(H(IYK(imain,K),1:nx,IS(1,2))*P1(imain,1:nx,L,I,J))

	  DO 55 K=1,nx
	  DO 55 L=1,iny
55	  RG(K,L) = SUM(R(K,1:nu,IS(1,6))*G(IYK(imain,L),1:nu,IS(1,3)))  ! R*G'	 

	  DO 56 K=1,iny
	  DO 56 L=1,iny
56	  COM(K,L)=SUM(H(IYK(imain,K),1:nx,IS(1,2))*RG(1:nx,L)) ! H*R*G'	 

	  DO 60 K=1,iny
        V(K,K) = SUM(HP1(K,1:nx)*H(IYK(imain,K),1:nx,IS(1,2)))
     #         + SUM(G(IYK(imain,K),1:nu,IS(1,3))
     #         * G(IYK(imain,K),1:nu,IS(1,3)))+2.*COM(K,K)
     
        DO 60 L=1,K-1
        V(K,L) = SUM(HP1(K,1:nx)*H(IYK(imain,L),1:nx,IS(1,2)))
     #         + SUM(G(IYK(imain,K),1:nu,IS(1,3))
     #         * G(IYK(imain,L),1:nu,IS(1,3)))+COM(K,L)+COM(L,K)
60      V(L,K) = V(K,L)
        
C -------------------------------------------------------------------
C  Updating equations:
C  x0 = x1 + (P1*H'+R*G')*Vinv*INN 
C  p0 = p1 - (P1*H'+R*G')*Vinv*(P1*H'+R*G')' 
C -------------------------------------------------------------------
        IF (iny.GT.0) THEN
         COM(1:iny,1:iny) = V(1:iny,1:iny)
	   IFAIL = -1
C	   CALL F01ADF(iny,COM(1:iny+1,1:iny),iny+1,IFAIL) 	 
         CALL DPOTRF('L',iny,COM(1:iny,1:iny),iny,IFAIL) ! COM = L*L'
         DETV = 1.D0 ! det(L)
         DO K=1,iny
	    DETV = DETV*COM(K,K)
         ENDDO
         CALL DPOTRI('L',iny,COM(1:iny,1:ny),iny,IFAIL) ! COM = VV^-1
         
	   DO 70 K=1,iny
	   VINV(K,K) = COM(K,K)
	   DO 70 L=1,K-1
	   VINV(K,L) = COM(K,L)
70	   VINV(L,K) = VINV(K,L) 

	   DO 90 K=1,nx
	   DO 90 L=1,iny	 
90       HPV(K,L) = SUM((HP1(1:iny,K)+RG(K,1:iny))*VINV(1:iny,L))

     	   DO 100 K=1,nx
100      X0(K,I,J)=X1(imain,K,I,J)+SUM(HPV(K,1:iny)*INNIJ(1:iny))

	   DO 110 K=1,nx
	   P0(K,K,I,J) = P1(imain,K,K,I,J)
     +               - SUM(HPV(K,1:iny)*(HP1(1:iny,K)+RG(K,1:iny)))
	   DO 110 L=1,K-1 	  
	   P0(K,L,I,J) = P1(imain,K,L,I,J)
     +               - SUM(HPV(K,1:iny)*(HP1(1:iny,L)+RG(L,1:iny)))
110	   P0(L,K,I,J) = P0(K,L,I,J)

C ---------------------------------------------                 
C  log f(y(t)|S(t-1)=i,S(t)=j,y^(t-1))
C  Log-Likelihood = -(RSS + ln(det(V))/2
C	        RSS = INN'*VINV*INN
C ---------------------------------------------
C	   IFAIL=-1
C        CALL F03ABF(V(1:iny,1:iny),iny,iny,DETV,COM(1:iny,1),IFAIL)
	   RSS = 0.D0
	   DO 120 K=1,iny
	   DO 120 L=1,iny
120      RSS = RSS + INNIJ(K)*VINV(K,L)*INNIJ(L)
       
         lfy=-.5D0*(RSS+2.*DLOG(DETV)+iny*DLOG(2.*3.141592653589793D0))
         fyss(I,J) = lfy + DLOG(P(J,I)) + DLOG(SUM(SP0(:,I))) 
        
        ELSE

         X0(:,I,J)   = X1(imain,:,I,J)
         P0(:,:,I,J) = P1(imain,:,:,I,J)
            
        ENDIF    
        
       ENDDO
	 ENDDO

C log f(y(t)|y^(t-1))
	 IF (iny.GT.0) THEN 
        IMAX = MAXLOC(fyss(:,1))
	  maxfyss = fyss(IMAX(1),1)
	  DO 130 K = 2,nk	
	  IMAX    = MAXLOC(fyss(:,K))
130	  maxfyss = MAX(fyss(IMAX(1),K),maxfyss)

	  fyss(:,:) = fyss(:,:) - maxfyss
	 
	  fsum = 0.D0
	  DO 140 I=1,nk
	  DO 140 J=1,nk
140     fsum = fsum + dexp(fyss(I,J))

        LIKE(imain) = maxfyss+DLOG(fsum)
        
C  Compute SP1 = Pr(S(t-1)=i,S(t)=j|y^t)
        SP1(:,:) = DEXP(fyss(:,:)-DLOG(fsum))
        
       ELSE

C  Compute SP1 = Pr(S(t-1)=i,S(t)=j|y^t-1)=Pr(S(t)=j|S(t-1)=i)*Pr(S(t-1)=i|y^t-1)
       DO I = 1,nk
        DO J = 1,nk   
         SP1(I,J) = P(J,I)*SFILT(imain-1,I)
        ENDDO
       ENDDO 
       ENDIF    
       
C  Compute Pr(S(t)=j|y^t)
       DO 145 J = 1,nk
145    SFILT(imain,J) = SUM(SP1(:,J))

C Shrinking XJ
	 DO 150 J = 1,nk
	 DO 150 K = 1,nx
150	 XI(imain,K,J) = SUM(SP1(:,J)*X0(K,:,J))/SFILT(imain,J)

C Shrinking PJ
	 DO 160 I=1,nk
	 DO 160 J=1,nk
       DO 160 K=1,nx
	 DO 160 L=1,nx
160    AUX1(K,L,I,J) = P0(K,L,I,J) 
     #     + (XI(imain,K,J)-X0(K,I,J))*(XI(imain,L,J)-X0(L,I,J))

       DO 170 J=1,nk
	 DO 170 K=1,nx
	 DO 170 L=1,nx
170    PI(imain,K,L,J) = SUM(SP1(:,J)*AUX1(K,L,:,J))/SFILT(imain,J) 

C Computing E(x(t)|y^t)
	 DO 175 J=1,nk
175    XFILT(imain,:)=XFILT(imain,:)+XI(imain,:,J)*SFILT(imain,J)

1000  SP0(:,:) = SP1(:,:)
      
      DEALLOCATE(SP1,SP0,X0,P0,AUX1,fyss,VINV,INNIJ)

C Smoothing S and X
	IF (ismoother.EQ.1) THEN
       ALLOCATE(SPT1(nk,nk),XS(nx,nk,nk),PS(nx,nx,nk,nk),XSC(nx,nk),
     1 PSC(nx,nx,nk),PTIL(nx,nx),W(nk))

	 SSMOOTH(nobs,:) = SFILT(nobs,:)
	 XSMOOTH(nobs,:) = XFILT(nobs,:)
       DO I=1,nx       
        XSSE(nobs,I) = dsqrt(SUM(SFILT(nobs,1:nk)*PI(nobs,I,I,1:nk)))
       ENDDO
       XSC(:,:)   = XI(nobs,:,:)
	 PSC(:,:,:) = PI(nobs,:,:,:)
	 DO 2000 imain = nobs-1,1,-1           
C SPT1:= Pr(S(t)=j,S(t+1)=k|y^T)=Pr(S(t+1)=k|y^T)*Pr(S(t)=j|y^t)*Pr(S(t+1)=k|S(t)=j)/Pr(S(t+1)=k|y^t) - (2.20')
C Pr(S(t+1)=k|y^t) = sum_i Pr(S(t+1)=k|S(t)=i)*Pr(S(t)=i|y^t) 
C XS:= E(x(t)|S(t)=j,S(t+1)=k,y^T)          nx x nk x nk       (2.24)
C PS:= E[(x(t)-XS)^2|S(t)=j,S(t+1)=k,y^T]   nx x nx x nk x nk  (2.25)
        DO J=1,nk         
         DO K=1,nk
          CALL INT2SEQ(K,nv,INFOS,SEQ,IS(1,:))        
          SPT1(J,K) = SSMOOTH(imain+1,K)*SFILT(imain,J)*P(K,J)
     #	          / SUM(P(K,1:nk)*SFILT(imain,1:nk))
C PTIL=(PI(imain,:,:,J)*PHI((K-1)*nx+1:K*nx,(K-1)*nx+1:K*nx)')/P1(imain+1,:,:,J,K)
C the traspose is stored         
	    DO 180 I=1,nx
          DO 180 L=1,nx
180       PTIL(L,I) = SUM(PI(imain,I,:,J)*F(L,:,IS(1,5))) 

C  P1(imain+1,:,:,J,K) * PTIL'  = PHI((K-1)*nx+1:K*nx,(K-1)*nx+1:K*nx)*PI(imain,:,:,J)'  
C         nx,nx         nx,nx     nx,nx 
	    AUX(:,:) = P1(imain+1,:,:,J,K)
C	    CALL F07ADF(nx,nx,AUX,nx,IPIV(1:nx),IFAIL)
C	    CALL F07AEF('N',nx,nx,AUX,nx,IPIV(1:nx),PTIL,nx,IFAIL) ! this gives PTIL'
          CALL DGETRF(nx,nx,AUX,nx,IPIV(1:nx),IFAIL)
          CALL DGETRS('N',nx,nx,AUX,nx,IPIV(1:nx),PTIL,nx,IFAIL) ! this gives PTIL'
          
C  XS(:,J,K) = XI(imain,:,J) + PTIL*(XSC(:,K)-X1(imain+1,:,J,K)')
	    DO 190 I=1,nx
190	    XS(I,J,K) = XI(imain,I,J) 
     #	          + SUM(PTIL(:,I)*(XSC(:,K)-X1(imain+1,:,J,K)))
         
C  PS(:,:,J,K) = PI(imain,:,:,J) + PTIL*(PSC(:,:,K)-P1(imain+1,:,:,J,K))*PTIL';           
	    DO 200 I=1,nx
          DO 200 L=1,nx
200       AUX(I,L) = SUM(PTIL(:,I)*(PSC(:,L,K)-P1(imain+1,:,L,J,K)))    
	    
          DO 210 I=1,nx
	    PS(I,I,J,K) = PI(imain,I,I,J)+SUM(AUX(I,:)*PTIL(:,I))       
          DO 210 L=1,I-1
          PS(I,L,J,K) = PI(imain,I,L,J)+SUM(AUX(I,:)*PTIL(:,L))
210	    PS(L,I,J,K) = PS(I,L,J,K)                
         ENDDO        

C  SSMOOTH Kim eqn (2.21)        
	   SSMOOTH(imain,J) = SUM(SPT1(J,:)) 
         
         IF (SSMOOTH(imain,J).GT.10.D-12) THEN
          DO I=1,nk
           W(I)=SPT1(J,I)/SSMOOTH(imain,J)
          ENDDO
         ELSE
          W(1:nk)=1.D0/DFLOAT(nk)   
         ENDIF

C  XSC(:,J) = XS(:,J,1:nk)*SPT1(J,1:nk)'/SSMOOTH(imain,J)  Kim eqn (2.26)
         DO 220 I=1,nx       
220      XSC(I,J) = SUM(XS(I,J,1:nk)*W(1:nk))

C  PSC(:,:,J) = PSC(:,:,J)+SPT1(J,K)*(PS(:,:,J,K)+AUX)/SSMOOTH(imain,J) Kim eqn (2.27)
C  AUX = (XSC(:,J)-X1(imain,:,J,K))*(XSC(:,J)-X1(imain,:,J,K))'        
         PSC(:,:,J) = 0.D0
         DO 235 K = 1,nk
          HPV(1:nx,1) = XSC(1:nx,J)-X1(imain,1:nx,J,K)	   
          DO 230 I=1,nx
	     AUX(I,I) = HPV(I,1)*HPV(I,1)
           DO 230 L=1,I-1
	     AUX(I,L) = HPV(I,1)*HPV(L,1)
230	     AUX(L,I) = AUX(I,L)
235        PSC(:,:,J) = PSC(:,:,J) + W(K)*(PS(:,:,J,K) + AUX(:,:))
        ENDDO

C  Kim eqn (2.28)       
	  DO 240 I=1,nx
        XSMOOTH(imain,I) = SUM(SSMOOTH(imain,1:nk)*XSC(I,1:nk))   
240     XSSE(imain,I)    = dsqrt(SUM(SSMOOTH(imain,1:nk)*PSC(I,I,1:nk)))
	
2000   CONTINUE
	 DEALLOCATE(X1,P1,XI,PI,AUX,SPT1,XS,PS,XSC,PSC,PTIL,XFILT,SFILT,P,W)

      ENDIF
	RETURN
	END

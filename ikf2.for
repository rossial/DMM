C -----------------C --------------------------------------------------------------------
C IKF2 (no missing values) COMPUTES INITIAL VALUES FOR THE KALMAN RECURSIONS FOR STATIONARY
C AND NON-STATIONARY TIME SERIES MODELS. 
C Stationary: Unconditional mean and variance of x1 as in 
C             A.C.Harvey 1989,"Forecasting strucural time series 
C             models and the Kalman filter", p.121 
C
C Non-stationary: Filtered state estimates and covaiance matrices 
C                 as in S.J.Koopman (1997), "Exact initial Kalman 
C                 Filtering and Smoothing for non-statonary Time 
C                 Series models", JASA, 92, pp.1630-38 
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
C H(t) (ny x nx x ns2)   ns2 = # of states for H(t)  
C G(t) (ny x nu x ns3)   ns3 = # of states for G(t)  
C a(t) (nx x ns4)        ns4 = # of states for a(t)  
C F(t) (nx x nx x ns5)   ns5 = # of states for F(t)  
C R(t) (nx x nu x ns6)   ns6 = # of states for R(t)  
C
C d(1): order of integration of the system
C d(2): number of non-stationary elements  
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
      SUBROUTINE IKF2(d,ny,nz,nx,nu,ns,S,yk,c,H,G,a,F,R,
	1                Xdd,Pdd,LIKE)  
C INPUT
	INTEGER d(2),ny,nz,nx,nu,ns(6),S(MAX(d(1),1),6)
	DOUBLE PRECISION yk(MAX(d(1),1),ny+nz),c(ny,max(1,nz),ns(1)),
	1 H(ny,nx,ns(2)),G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),
     2 R(nx,nu,ns(6))
C OUTPUT
	DOUBLE PRECISION Pdd(MAX(d(1),1),nx,nx),Xdd(MAX(d(1),1),nx),
	1 LIKE(MAX(d(1),1))
C LOCALS
	INTEGER I,J,IFAIL,imain,FiRANK
	INTEGER IPIV(nx)
	DOUBLE PRECISION aa(nx),Ps(nx,nx),Pi(nx,nx)
	DOUBLE PRECISION HPs(ny,nx),HPi(ny,nx),
	1 Fi(ny,ny),Fs(ny,ny),Fim(ny,ny),Fsm(ny,ny),
     2 PHFs(nx,ny),PHFi(nx,ny),FFF(ny,ny),Mi(nx,ny),Ci(nx,nx)
      DOUBLE PRECISION W1(ny),WORK(64*nx),WORK1(64*ny),
	1 PFP(nx,nx),APPO(nx,nx),APPO1(nx,ny),COM(ny+1,ny),RG(nx,ny)
     	DOUBLE PRECISION ONE,ZERO,DETV,RSS,SUMW1
	DATA ONE/1.0D0/,ZERO/0.0D0/
	
C Unconditional mean and variance
	LIKE(:) = ZERO
      IF (d(1).EQ.0) THEN ! stationary models
	 IF(SUM(ABS(a(:,S(1,4)))).EQ.ZERO) THEN
	  Xdd(1,:) = ZERO
	 ELSE
	  APPO = -F(:,:,S(1,5))	 
	  DO 1 I = 1,nx
1	  APPO(I,I) = 1.D0+APPO(I,I)	 
        IFAIL = -1
C	  CALL F07ADF(nx,nx,APPO,nx,IPIV,IFAIL)
C	  CALL F07AJF(nx,APPO,nx,IPIV,WORK,64*nx,IFAIL)        
	  CALL DGETRF(nx,nx,APPO,nx,IPIV,IFAIL)
	  CALL DGETRI(nx,APPO,nx,IPIV,WORK,64*nx,IFAIL)        
	  DO 3 I =1,nx
3	  Xdd(1,I) = SUM(APPO(I,:)*a(:,S(1,4))) ! inv(I-F)*a
       ENDIF

C Pdd - F*Pdd*F' = R*R' 
	 CALL LYAP(nx,nu,1.D-3,F(:,:,S(1,5)),R(:,:,S(1,6)),Pdd)
	ELSE  
C -----------------------------------------------------------
C Non-stationary models
C Define X(1) = aa + A*eta + B*delta (A*B' = 0) 
C               eta~N(0,I), delta~N(0,k*I) k -> +inf  
C               X(1)~N(aa,P), P=Ps+k*Pi, Ps=AA', Pi=BB'. 
C -----------------------------------------------------------
	 aa(1:nx)      = ZERO  
	 Ps(1:nx,1:nx) = ZERO
	 IF (d(2).LT.nx) THEN
	  IF(SUM(ABS(a(d(2)+1:nx,S(1,4)))).NE.0.D0) THEN
	   APPO(d(2)+1:nx,d(2)+1:nx) = -F(d(2)+1:nx,d(2)+1:nx,S(1,5))	 
	   DO 5 I = d(2)+1,nx
5	   APPO(I,I) = 1.D0+APPO(I,I)	 
C	   CALL F07ADF(nx-d(2),nx-d(2),APPO(d(2)+1:nx,d(2)+1:nx),nx-d(2),
C	1               IPIV(d(2)+1:nx),IFAIL)
C	   CALL F07AJF(nx-d(2),APPO(d(2)+1:nx,d(2)+1:nx),nx-d(2),
C	1               IPIV(d(2)+1:nx),WORK,64*nx,IFAIL)
	   CALL DGETRF(nx-d(2),nx-d(2),APPO(d(2)+1:nx,d(2)+1:nx),nx-d(2),
	1               IPIV(d(2)+1:nx),IFAIL)
	   CALL DGETRI(nx-d(2),APPO(d(2)+1:nx,d(2)+1:nx),nx-d(2),
	1               IPIV(d(2)+1:nx),WORK,64*nx,IFAIL)
	   DO 6 I = d(2)+1,nx
6	   aa(I) = SUM(APPO(I,d(2)+1:nx)*a(d(2)+1:nx,S(1,4))) ! inv(I-F)*a
        ENDIF

C Lyapunov eqn
	  CALL LYAP(nx-d(2),nu,1.D-3,F(d(2)+1:nx,d(2)+1:nx,S(1,5)),
	1            R(d(2)+1:nx,1:nu,S(1,6)),Ps(d(2)+1:nx,d(2)+1:nx))

       ENDIF

	 Pi(1:nx,1:nx) = ZERO 
	 DO 10 I = 1,d(2)
10	 Pi(I,I) = ONE

       Xdd(:,:)   = ZERO
	 Pdd(:,:,:) = ZERO
	 DO 1000 imain = 1,d(1)
	  DO 30 I=1,ny
        DO 30 J=1,nx
30      HPs(I,J) = SUM(H(I,:,S(imain,2))*Ps(:,J)) 

	  DO 40 I=1,ny
	   Fs(I,I) = SUM(HPs(I,:)*H(I,:,S(imain,2)))
     +   +SUM(G(I,:,S(imain,3))*G(I,:,S(imain,3))) 
         DO 40 J=1,I-1
         Fs(I,J) = SUM(HPs(I,:)*H(J,:,S(imain,2)))
     +   +SUM(G(I,:,S(imain,3))*G(J,:,S(imain,3))) 
40       Fs(J,I) = Fs(I,J)

	  DO 50 I=1,ny
        DO 50 J=1,nx
50      HPi(I,J) = SUM(H(I,:,S(imain,2))*Pi(:,J))  

	  DO 60 I=1,ny
         Fi(I,I) = SUM(HPi(I,:)*H(I,:,S(imain,2)))
	   DO 60 J=1,I-1
         Fi(I,J) = SUM(HPi(I,:)*H(J,:,S(imain,2)))
60       Fi(J,I) = Fi(I,J)

C --------------------------------------------------------------------------
C Computes inverse of the innovation variance matrix
C Cases: ny = 1, Fi is scalar >0 (or 0 not considered)
C        ny > 1, Fi is full rank or singular (or 0 matrix not considered)
C --------------------------------------------------------------------------
        IF (ny.EQ.1) THEN
	   Fsm = ZERO
	   Fim = 1.D0/Fi
	   FFF = Fim*Fs*Fim 
	  ELSE
	   
	  IFAIL = -1
	  COM(1:ny,1:ny) = Fi(1:ny,1:ny)  
C	  CALL F02FAF('N','U',ny,COM(1:ny,1:ny),ny,W1(1:ny),WORK1,64*ny,IFAIL)
        CALL DSYEV('N','U',ny,COM(1:ny,1:ny),ny,W1(1:ny),WORK1,
     1             64*ny,IFAIL)
	  FiRANK = 0
	  SUMW1  = SUM(ABS(W1(1:ny)))
	  DO 70 I=1,ny
	  W1(I) = W1(I)/SUMW1
70      IF (W1(I).GT.1.D-10) FiRANK=FiRANK+1
	  FiRANK = min(FiRANK,d(2))

	  IF(FiRANK.EQ.ny) THEN
	    Fsm = ZERO
	    COM(1:ny,1:ny) = Fi(1:ny,1:ny)
		IFAIL = -1
C		CALL F01ADF(ny,COM(1:ny+1,1:ny),ny+1,IFAIL) 
          CALL DPOTRF('L',ny,COM(1:ny,1:ny),ny,IFAIL) ! COM = L*L'
          CALL DPOTRI('L',ny,COM(1:ny,1:ny),ny,IFAIL) ! COM = VV^-1

		DO 80 I=1,ny
       	Fim(I,I) = COM(I,I)
		DO 80 J=1,I-1
		Fim(I,J) = COM(I,J)
80		Fim(J,I) = Fim(I,J)           	   
          
		DO 81 I=1,ny
		DO 81 J=1,ny
81		COM(I,J) = SUM(Fim(I,1:ny)*Fs(1:ny,J))  ! Fim x Fs
 
		DO 82 I=1,ny
		 FFF(I,I) = SUM(COM(I,1:ny)*Fim(1:ny,I))
		 DO 82 J=1,I-1
	     FFF(I,J) = SUM(COM(I,1:ny)*Fim(1:ny,J))  ! Fim x Fs x Fim
82		 FFF(J,I) = FFF(I,J)

	   ELSE 
	    SUMW1=0.D0
		DO I=Firank+1,ny
		 SUMW1 = SUMW1 + Fi(I,I)
          ENDDO
          IF (SUMW1.GT.0.D0) THEN
		 CALL INVFBIS(Fs(1:ny,1:ny),Fi(1:ny,1:ny),ny,FiRANK,
	1                  Fsm(1:ny,1:ny),Fim(1:ny,1:ny),FFF(1:ny,1:ny))
          ELSE
           CALL INVF(Fs(1:ny,1:ny),Fi(1:ny,1:ny),ny,FiRANK,
	1               Fsm(1:ny,1:ny),Fim(1:ny,1:ny),FFF(1:ny,1:ny)) 
          ENDIF
	   ENDIF 	
	  ENDIF
C ------------------------------------------------------------------
C X(d|d) = X(d|d-1)+((Ps*H'+R*G')*Fsm+Pi*H'*Fim)*(Y(d)-H*X(d|d-1)-c)  		
C ------------------------------------------------------------------
	  DO 85 I = 1,nx
	  DO 85 J = 1,ny        
 	  RG(I,J) = 
     #	 SUM(R(I,1:nu,S(imain,6))*G(J,1:nu,S(imain,3)))
85      HPs(J,I) = HPs(J,I) + RG(I,J) ! HPs = (Ps*H'+R*G')'

	  DO 90 I = 1,nx
	  DO 90 J = 1,ny        
 	  PHFs(I,J) = SUM(HPs(1:ny,I)*Fsm(1:ny,J))
90	  PHFi(I,J) = SUM(HPi(1:ny,I)*Fim(1:ny,J))
	
C Innovations
	  DO 100 I=1,ny
100	  COM(I,1) = yk(imain,I)
     +    - SUM(H(I,1:nx,S(imain,2))*aa(1:nx))
     +    - SUM(c(I,1:nz,S(imain,1))*yk(imain,ny+1:ny+nz)) 
      	  
	  DO 110 I=1,nx
110     Xdd(imain,I) = aa(I)
     +               + SUM((PHFs(I,1:ny)+PHFi(I,1:ny))*COM(1:ny,1))

C P(d|d) = P(d|d-1) + Pi*H'*Fim*Fs*Fim*H*Pi - Ps*H'*Fsm*H*Ps 
C        - Ps*H'*Fim*H*Pi - (Ps*H'*Fim*H*Pi)' 

C- Ps*H'*Fsm*H*Ps
	  DO 120 I = 1,nx
	   APPO(I,I) = -SUM(PHFs(I,1:ny)*HPs(1:ny,I))
	   DO 120 J = 1,I-1
	   APPO(I,J) = -SUM(PHFs(I,1:ny)*HPs(1:ny,J))
120	   APPO(J,I) = APPO(I,J)

C - Ps*H'*Fim*H*Pi - (Ps*H'*Fim*H*Pi)' 
	  DO 130 I = 1,nx
	   APPO(I,I) = APPO(I,I) - SUM(HPs(1:ny,I)*PHFi(I,1:ny))
     +             - SUM(PHFi(I,1:ny)*HPs(1:ny,I))
	   DO 130 J = 1,I-1
	   APPO(I,J) = APPO(I,J) - SUM(HPs(1:ny,I)*PHFi(J,1:ny))
     +             - SUM(PHFi(I,1:ny)*HPs(1:ny,J))
130	   APPO(J,I) = APPO(I,J)
	
C Pi*H'*Fim*Fs*Fim*H*Pi
	  DO 140 I = 1,nx
	  DO 140 J = 1,ny
140	  APPO1(I,J) = SUM(HPi(1:ny,I)*FFF(1:ny,J))

	  DO 150 I = 1,nx
	   PFP(I,I) = SUM(APPO1(I,1:ny)*HPi(1:ny,I))
	   DO 150 J = 1,I-1
	   PFP(I,J) = SUM(APPO1(I,1:ny)*HPi(1:ny,J))
150      PFP(J,I) = PFP(I,J)

	  Pdd(imain,:,:) = Ps(:,:) + PFP(:,:) + APPO(:,:)  

C ----------------------------------------------
C CONTRIBUTE TO THE LIKELIHOOD 1ST d INNOVATIONS
C ----------------------------------------------
	  IFAIL = -1
C       CALL F03ABF(Fsm(1:ny,1:ny)+Fim(1:ny,1:ny),ny,ny,
C	1              DETV,WORK1(1:ny),IFAIL)
        FFF(1:ny,1:ny) = Fsm(1:ny,1:ny)+Fim(1:ny,1:ny)
        CALL DPOTRF('L',ny,FFF(1:ny,1:ny),ny,IFAIL) ! FFF = L*L'       
	  DETV = 1.D0 
	  RSS = ZERO
	  DO 155 I=1,ny
        DETV = DETV*FFF(I,I)    
	  DO 155 J=1,ny
155     RSS = RSS + COM(I,1)*Fsm(I,J)*COM(J,1)
	
	  LIKE(imain) = -.5D0*(RSS - 2.D0*DLOG(DETV))
	  IF (LIKE(imain).NE.0.D0) THEN 
	   LIKE(imain)=LIKE(imain)-ny/2.D0*DLOG(2.*3.141592653589793D0)	    
	  ENDIF
C ----------------------------------
C Predictions X(d+1|d) and P(d+1|d)
C ----------------------------------
	  IF (imain.LT.d(1)) THEN
C aa = a + F*Xdd
	   DO 160 I=1,nx
160	   aa(I) = a(I,S(imain+1,4))+SUM(F(I,:,S(imain+1,5))*Xdd(imain,:))

C Pi = F*Pi*F'-Ci
C Ps = F*PddF'+R*R'
	   DO 170 I = 1,nx
	   DO 170 J = 1,nx	   
	   PFP(I,J)  = SUM(F(I,:,S(imain+1,5))*Pi(:,J))        ! F*Pi
170	   APPO(I,J) = SUM(F(I,:,S(imain+1,5))*Pdd(imain,:,J)) ! F*Pdd

C Mi = F*Pi*H'  ! H to be checked 
	   DO 172 I = 1,nx
	   DO 172 J = 1,ny	   
172	   Mi(I,J) = SUM(PFP(I,1:nx)*H(J,1:nx,S(imain,2)))!S(+1,..) before

C Ci = Mi*Fim*Mi'
	   DO 174 I = 1,nx
	   DO 174 J = 1,ny	   
174	   RG(I,J) = SUM(Mi(I,1:ny)*Fim(1:ny,J))  ! Mi*Fim
	   
	   DO 176 I = 1,nx
	    Ci(I,I) = SUM(RG(I,1:ny)*Mi(I,1:ny))
	    DO 176 J = 1,I-1
176		Ci(I,J) = SUM(RG(I,1:ny)*Mi(J,1:ny))
	   
	   DO 180 I = 1,nx
	    Pi(I,I) = SUM(PFP(I,1:nx)*F(I,1:nx,S(imain+1,5)))-Ci(I,I) 
	    Ps(I,I) = SUM(APPO(I,1:nx)*F(I,1:nx,S(imain+1,5))) 
     +            + SUM(R(I,1:nu,S(imain+1,6))*R(I,1:nu,S(imain+1,6)))
	    DO 180 J = 1,I-1	   
	    Pi(I,J) = SUM(PFP(I,1:nx)*F(J,1:nx,S(imain+1,5)))-Ci(I,J) 
	    Ps(I,J) = SUM(APPO(I,1:nx)*F(J,1:nx,S(imain+1,5))) 
     +            + SUM(R(I,1:nu,S(imain+1,6))*R(J,1:nu,S(imain+1,6)))
          Pi(J,I) = Pi(I,J)
180	    Ps(J,I) = Ps(I,J)  

        ENDIF
1000   CONTINUE
      ENDIF
      
      RETURN
      END
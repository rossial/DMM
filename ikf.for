C --------------------------------------------------------------------
C IKF COMPUTES INITIAL VALUES FOR THE KALMAN RECURSIONS FOR STATIONARY
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
      SUBROUTINE IKF(d,ny,nz,nx,nu,ns,S,yk,IYK,c,H,G,a,F,R,
	1               Xdd,Pdd,LIKE)  
C INPUT
	INTEGER d(2),ny,nz,nx,nu,ns(6),S(MAX(d(1),1),6),
	1 IYK(max(d(1),1),ny+1)
	DOUBLE PRECISION yk(MAX(d(1),1),ny+nz),c(ny,max(1,nz),ns(1)),
	1 H(ny,nx,ns(2)),G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),
     2 R(nx,nu,ns(6))
C OUTPUT
	DOUBLE PRECISION Pdd(MAX(d(1),1),nx,nx),Xdd(MAX(d(1),1),nx),
	1 LIKE(MAX(d(1),1))
C LOCALS
	INTEGER I,J,IFAIL,imain,iny,FiRANK,NULLITY
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
	  CALL F07ADF(nx,nx,APPO,nx,IPIV,IFAIL)
	  CALL F07AJF(nx,APPO,nx,IPIV,WORK,64*nx,IFAIL)        
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
	   CALL F07ADF(nx-d(2),nx-d(2),APPO(d(2)+1:nx,d(2)+1:nx),nx-d(2),
	1               IPIV(d(2)+1:nx),IFAIL)
	   CALL F07AJF(nx-d(2),APPO(d(2)+1:nx,d(2)+1:nx),nx-d(2),
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
	  iny = IYK(imain,ny+1)
	  DO 30 I=1,iny
        DO 30 J=1,nx
30      HPs(I,J) = SUM(H(IYK(imain,I),:,S(imain,2))*Ps(:,J)) 

	  DO 40 I=1,iny
	   Fs(I,I) = SUM(HPs(I,:)*H(IYK(imain,I),:,S(imain,2)))+
     +   SUM(G(IYK(imain,I),:,S(imain,3))*G(IYK(imain,I),:,S(imain,3))) 
         DO 40 J=1,I-1
         Fs(I,J) = SUM(HPs(I,:)*H(IYK(imain,J),:,S(imain,2)))+
     +   SUM(G(IYK(imain,I),:,S(imain,3))*G(IYK(imain,J),:,S(imain,3))) 
40       Fs(J,I) = Fs(I,J)

	  DO 50 I=1,iny
        DO 50 J=1,nx
50      HPi(I,J) = SUM(H(IYK(imain,I),:,S(imain,2))*Pi(:,J))  

	  DO 60 I=1,iny
         Fi(I,I) = SUM(HPi(I,:)*H(IYK(imain,I),:,S(imain,2)))
	   DO 60 J=1,I-1
         Fi(I,J) = SUM(HPi(I,:)*H(IYK(imain,J),:,S(imain,2)))
60       Fi(J,I) = Fi(I,J)

C --------------------------------------------------------------------------
C Computes inverse of the innovation variance matrix
C Cases: ny = 1, Fi is scalar >0 (or 0 not considered)
C        ny > 1, Fi is full rank or singular (or 0 matrix not considered)
C --------------------------------------------------------------------------
        IF (iny.EQ.1) THEN
	   Fsm = ZERO
	   Fim = 1.D0/Fi
	   FFF = Fim*Fs*Fim 
	  ELSE

	  IFAIL = -1
	  COM(1:iny,1:iny) = Fi(1:iny,1:iny)  
	  CALL F02FAF('N','U',iny,COM(1:iny,1:iny),iny,W1(1:iny),
	1              WORK1,64*iny,IFAIL)
	  FiRANK = 0
	  SUMW1  = SUM(ABS(W1(1:iny)))
	  DO 70 I=1,iny
	  W1(I) = W1(I)/SUMW1
70      IF (W1(I).GT.1.D-10) FiRANK=FiRANK+1
	  FiRANK = min(FiRANK,d(2))

	  IF(FiRANK.EQ.iny) THEN
	    Fsm = ZERO
	    COM(1:iny,1:iny) = Fi(1:iny,1:iny)
		IFAIL = -1
		CALL F01ADF(iny,COM(1:iny+1,1:iny),iny+1,IFAIL) 
		DO 80 I=1,iny
		 Fim(I,I) = COM(I+1,I)
		 DO 80 J=1,I-1
		 Fim(I,J) = COM(I+1,J)
80		 Fim(J,I) = Fim(I,J)           	   
          
		DO 81 I=1,iny
		DO 81 J=1,iny
81		COM(I,J) = SUM(Fim(I,1:iny)*Fs(1:iny,J))  ! Fim x Fs
 
		DO 82 I=1,iny
		 FFF(I,I) = SUM(COM(I,1:iny)*Fim(1:iny,I))
		 DO 82 J=1,I-1
	     FFF(I,J) = SUM(COM(I,1:iny)*Fim(1:iny,J))  ! Fim x Fs x Fim
82		 FFF(J,I) = FFF(I,J)

	   ELSE 
	    SUMW1=0.D0
		DO I=Firank+1,iny
		 SUMW1 = SUMW1 + Fi(I,I)
          ENDDO
          IF (SUMW1.GT.0.D0) THEN
		 CALL INVFBIS(Fs(1:iny,1:iny),Fi(1:iny,1:iny),iny,FiRANK,
	1          Fsm(1:iny,1:iny),Fim(1:iny,1:iny),FFF(1:iny,1:iny))
		ELSE
		 CALL INVF(Fs(1:iny,1:iny),Fi(1:iny,1:iny),iny,FiRANK,
	1          Fsm(1:iny,1:iny),Fim(1:iny,1:iny),FFF(1:iny,1:iny)) 
          ENDIF	 
	   ENDIF 	
	  ENDIF
C ------------------------------------------------------------------
C X(d|d) = X(d|d-1)+((Ps*H'+R*G')*Fsm+Pi*H'*Fim)*(Y(d)-H*X(d|d-1)-c)  		
C ------------------------------------------------------------------
	  DO 85 I = 1,nx
	  DO 85 J = 1,iny        
 	  RG(I,J) = 
     #	 SUM(R(I,1:nu,S(imain,6))*G(IYK(imain,J),1:nu,S(imain,3)))
85      HPs(J,I) = HPs(J,I) + RG(I,J) ! HPs = (Ps*H'+R*G')'

	  DO 90 I = 1,nx
	  DO 90 J = 1,iny        
 	  PHFs(I,J) = SUM(HPs(1:iny,I)*Fsm(1:iny,J))
90	  PHFi(I,J) = SUM(HPi(1:iny,I)*Fim(1:iny,J))
	
C Innovations
	  DO 100 I=1,iny
100	  COM(I,1) = yk(imain,IYK(imain,I))
     +    - SUM(H(IYK(imain,I),1:nx,S(imain,2))*aa(1:nx))
     +    - SUM(c(IYK(imain,I),1:nz,S(imain,1))*yk(imain,ny+1:ny+nz)) 
      	  
	  DO 110 I=1,nx
110     Xdd(imain,I) = aa(I)
     +               + SUM((PHFs(I,1:iny)+PHFi(I,1:iny))*COM(1:iny,1))

C P(d|d) = P(d|d-1) + Pi*H'*Fim*Fs*Fim*H*Pi - Ps*H'*Fsm*H*Ps - Ps*H'*Fim*H*Pi - (Ps*H'*Fim*H*Pi)'
C - Ps*H'*Fsm*H*Ps
	  DO 120 I = 1,nx
	   APPO(I,I) = -SUM(PHFs(I,1:iny)*HPs(1:iny,I))
	   DO 120 J = 1,I-1
	   APPO(I,J) = -SUM(PHFs(I,1:iny)*HPs(1:iny,J))
120	   APPO(J,I) = APPO(I,J)

C - Ps*H'*Fim*H*Pi - (Ps*H'*Fim*H*Pi)' 
	  DO 130 I = 1,nx
	   APPO(I,I) = APPO(I,I) - SUM(HPs(1:iny,I)*PHFi(I,1:iny))
     +             - SUM(PHFi(I,1:iny)*HPs(1:iny,I))
	   DO 130 J = 1,I-1
	   APPO(I,J) = APPO(I,J) - SUM(HPs(1:iny,I)*PHFi(J,1:iny))
     +             - SUM(PHFi(I,1:iny)*HPs(1:iny,J))
130	   APPO(J,I) = APPO(I,J)
	
C Pi*H'*Fim*Fs*Fim*H*Pi
	  DO 140 I = 1,nx
	  DO 140 J = 1,iny
140	  APPO1(I,J) = SUM(HPi(1:iny,I)*FFF(1:iny,J))

	  DO 150 I = 1,nx
	   PFP(I,I) = SUM(APPO1(I,1:iny)*HPi(1:iny,I))
	   DO 150 J = 1,I-1
	   PFP(I,J) = SUM(APPO1(I,1:iny)*HPi(1:iny,J))
150      PFP(J,I) = PFP(I,J)

	  Pdd(imain,:,:) = Ps(:,:) + PFP(:,:) + APPO(:,:)  

C ----------------------------------------------
C CONTRIBUTE TO THE LIKELIHOOD 1ST d INNOVATIONS
C ----------------------------------------------
	  IFAIL = -1
        CALL F03ABF(Fsm(1:iny,1:iny)+Fim(1:iny,1:iny),iny,iny,
	1              DETV,WORK1(1:iny),IFAIL)
	  
	  RSS = ZERO
	  DO 155 I=1,iny
	  DO 155 J=1,iny
155     RSS = RSS + COM(I,1)*Fsm(I,J)*COM(J,1)
	
	  LIKE(imain) = -.5D0*(RSS - DLOG(DETV))
	  IF (LIKE(imain).NE.0.D0) THEN 
	   LIKE(imain)=LIKE(imain)-iny/2.D0*DLOG(2.*3.141592653589793D0)	    
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
	   DO 172 J = 1,iny	   
172	   Mi(I,J) = SUM(PFP(I,1:nx)*H(IYK(imain,J),1:nx,S(imain,2)))!S(+1,..) before

C Ci = Mi*Fim*Mi'
	   DO 174 I = 1,nx
	   DO 174 J = 1,iny	   
174	   RG(I,J) = SUM(Mi(I,1:iny)*Fim(1:iny,J))  ! Mi*Fim
	   
	   DO 176 I = 1,nx
	    Ci(I,I) = SUM(RG(I,1:iny)*Mi(I,1:iny))
	    DO 176 J = 1,I-1
176		Ci(I,J) = SUM(RG(I,1:iny)*Mi(J,1:iny))
	   
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

c	  IF (IFAIL.NE.0) THEN ! pseudo-inverse
c        IFAIL = -1
c	   RSS = 1.D-13
c	   CALL F01BLF(nx,nx,RSS,Ci,nx,WORK(1:nx),IRANK,INC,
c	1   WORK(nx+1:2*nx),APPO,nx,WORK(2*nx+1:3*nx),IFAIL)
c	   DO 2 I = 1,nx
c	   DO 2 J = 1,nx
c2        APPO(I,J) = Ci(J,I)
c	  ELSE                 ! normal-inverse

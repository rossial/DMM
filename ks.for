C --------------------------------------------------------------------
C KS IMPLEMENTS THE KALMAN SMMOOTHER RECURSIONS in 
C Koopman (1997), JASA, 92, 440, 1630-38
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C
C XS = E[x(t)|y(1),...,y(nobs)] 
C PS = V[x(t)|y(1),...,y(nobs)], t = 1,2,...,nobs
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
      SUBROUTINE KS(nobs,d,ny,nz,nx,nu,ns,S,yk,IYK,c,H,G,a,F,R,XS,PS)  
C INPUT
	INTEGER nobs,d(2),ny,nz,nx,nu,ns(6),S(nobs,6),IYK(nobs,ny+1)
	DOUBLE PRECISION yk(nobs,ny+nz),c(ny,max(nz,1),ns(1)),
	1 H(ny,nx,ns(2)),G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),
     2 R(nx,nu,ns(6))

C OUTPUT
	DOUBLE PRECISION XS(nobs,nx),PS(nobs,nx,nx)	

C LOCALS
	INTEGER imain,iny,I,J,IFAIL,FiRANK,ITIME
	INTEGER IPIV(nx)
	DOUBLE PRECISION TOL,SUMW1
	DOUBLE PRECISION,ALLOCATABLE:: Pi(:,:,:),HPs(:,:),HPi(:,:),
	1 Fi(:,:),Fs(:,:),Fim(:,:),Fsm(:,:),PHFs(:,:),PHFi(:,:),FFF(:,:),
     2 Mi(:,:),Ms(:,:),Ci(:,:),Kst(:,:,:),Kit(:,:,:),W1(:),WORK(:),
     3 WORK1(:),PFP(:,:),APPO(:,:),APPO1(:,:),COM(:,:),RG(:,:),
	4 FP(:,:),HP1(:,:),V(:,:),CC(:,:),HPV(:,:),X0(:),P0(:,:),
     5 RECR(:),RECRI(:),RECN(:,:),XT(:,:),PT(:,:,:),INN(:,:),
	1 Vinv(:,:,:),Vis(:,:,:),Vii(:,:,:)
	
	ALLOCATE(Pi(d(1),nx,nx),HPs(ny,nx),HPi(ny,nx),
	1 Fi(ny,ny),Fs(ny,ny),Fim(ny,ny),Fsm(ny,ny),
     2 PHFs(nx,ny),PHFi(nx,ny),FFF(ny,ny),Mi(nx,ny),Ms(nx,ny),Ci(nx,nx),
	3 Kst(d(1),nx,ny),Kit(d(1),nx,ny))
      
	ALLOCATE(W1(ny),WORK(64*nx),WORK1(64*ny),
	1 PFP(nx,nx),APPO(nx,nx),APPO1(nx,ny),COM(ny+1,ny),RG(nx,ny))
	
	ALLOCATE(FP(nx,nx),HP1(ny,nx),V(ny,ny),CC(nx,nx),
	1 HPV(nx,ny),X0(nx),P0(nx,nx),RECR(nx),RECRI(nx),RECN(nx,nx))
  
	ALLOCATE(XT(nobs,nx),PT(nobs,nx,nx),INN(nobs,ny),
	1 Vinv(nobs,ny,ny),Vis(d(1),ny,ny),Vii(d(1),ny,ny))

	TOL = 1.D-3
C Unconditional mean and variance
      IF (d(1).EQ.0) THEN ! stationary models
	 IF(SUM(ABS(a(:,S(1,4)))).EQ.0.D0) THEN
	  XT(1,:) = 0.D0  ! X(1|0)
	 ELSE
	  APPO = -F(:,:,S(1,5))	 
	  DO 1 I = 1,nx
1	  APPO(I,I) = 1.D0+APPO(I,I)	 
	  CALL F07ADF(nx,nx,APPO,nx,IPIV,IFAIL)
	  CALL F07AJF(nx,APPO,nx,IPIV,WORK,64*nx,IFAIL)
	  DO 2 I =1,nx
2	  XT(1,I) = SUM(APPO(I,:)*a(:,S(1,4))) ! inv(I-F)*a
       ENDIF

C P(1|0) - F*P(1|0)*F' = R*R'
	 CALL LYAP(nx,nu,TOL,F(:,:,S(1,5)),R(:,:,S(1,6)),PT(1,:,:))
	ELSE  
C -----------------------------------------------------------
C Non-stationary models
C Define X(1) = aa + A*eta + B*delta (A*B' = 0) 
C               eta~N(0,I), delta~N(0,k*I) k -> +inf  
C               X(1)~N(aa,P), P=Ps+k*Pi, Ps=AA', Pi=BB'. 
C CARE!! aa (uncond. mean),Ps, and Pi to be filled by users 
C -----------------------------------------------------------
	 XT(1,1:nx)      = 0.D0  ! X(1|0)
	 PT(1,1:nx,1:nx) = 0.D0  ! P(1|0)
	 IF (d(2).LT.nx) THEN
	  IF(SUM(ABS(a(d(2)+1:nx,S(1,4)))).NE.0.D0) THEN
	   APPO(d(2)+1:nx,d(2)+1:nx) = -F(d(2)+1:nx,d(2)+1:nx,S(1,5))	 
	   DO 3 I = d(2)+1,nx
3	   APPO(I,I) = 1.D0+APPO(I,I)	 
	   CALL F07ADF(nx-d(2),nx-d(2),APPO(d(2)+1:nx,d(2)+1:nx),nx-d(2),
	1               IPIV(d(2)+1:nx),IFAIL)
	   CALL F07AJF(nx-d(2),APPO(d(2)+1:nx,d(2)+1:nx),nx-d(2),
	1               IPIV(d(2)+1:nx),WORK,64*nx,IFAIL)
	   DO 4 I = d(2)+1,nx
4	   XT(1,I) = SUM(APPO(I,d(2)+1:nx)*a(d(2)+1:nx,S(1,4))) ! inv(I-F)*a
        ENDIF

C Lyapunov eqn
	  CALL LYAP(nx-d(2),nu,TOL,F(d(2)+1:nx,d(2)+1:nx,S(1,5)),
	1            R(d(2)+1:nx,1:nu,S(1,6)),PT(1,d(2)+1:nx,d(2)+1:nx))
       ENDIF

	 Pi(:,:,:) = 0.D0 
	 DO 5 I = 1,d(2)
5	 Pi(1,I,I) = 1.D0

	 DO 200 imain = 1,d(1)
	  iny = IYK(imain,ny+1)
	  DO 30 I=1,iny
        DO 30 J=1,nx
30      HPs(I,J) = SUM(H(IYK(imain,I),:,S(imain,2))*PT(imain,:,J)) 

	  DO 40 I=1,iny
	   Fs(I,I) = SUM(HPs(I,:)*H(IYK(imain,I),:,S(imain,2)))
     +   +SUM(G(IYK(imain,I),:,S(imain,3))*G(IYK(imain,I),:,S(imain,3))) 
         DO 40 J=1,I-1
         Fs(I,J) = SUM(HPs(I,:)*H(IYK(imain,J),:,S(imain,2)))
     +   +SUM(G(IYK(imain,I),:,S(imain,3))*G(IYK(imain,J),:,S(imain,3))) 
40       Fs(J,I) = Fs(I,J)

  	  DO 50 I=1,iny
        DO 50 J=1,nx
50      HPi(I,J) = SUM(H(IYK(imain,I),:,S(imain,2))*Pi(imain,:,J))  

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
	   Fsm = 0.D0
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
	    Fsm = 0.D0
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
	  Vis(imain,1:iny,1:iny) = Fsm(1:iny,1:iny)
	  Vii(imain,1:iny,1:iny) = Fim(1:iny,1:iny)  

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
100	  INN(imain,I) = yk(imain,IYK(imain,I)) 
     +   - SUM(H(IYK(imain,I),1:nx,S(imain,2))*XT(imain,1:nx)) 
     +   - SUM(c(IYK(imain,I),1:nz,S(imain,1))*yk(imain,ny+1:ny+nz)) 
                
	  DO 110 I=1,nx
110     X0(I) = XT(imain,I)
     +        + SUM((PHFs(I,1:iny)+PHFi(I,1:iny))*INN(imain,1:iny))

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

	  P0(:,:) = PT(imain,:,:) + PFP(:,:) + APPO(:,:)  

C Mi = F*Pi*H'
	  DO 151 I = 1,nx
	  DO 151 J = 1,nx	   
151	  PFP(I,J) = SUM(F(I,1:nx,S(imain,5))*Pi(imain,1:nx,J)) ! F*Pi

  	  DO 152 I = 1,nx
	  DO 152 J = 1,iny	   
152	  Mi(I,J) = SUM(PFP(I,1:nx)*H(IYK(imain,J),1:nx,S(imain,2)))

C Ms = F*Ps*H' + R*G'
	  DO 153 I = 1,nx
	  DO 153 J = 1,nx	   
153	  FP(I,J) = SUM(F(I,:,S(imain,5))*PT(imain,:,J)) ! F*Ps

	  DO 154 I = 1,nx
	  DO 154 J = 1,iny	   
154	  Ms(I,J)=RG(I,J)+SUM(FP(I,1:nx)*H(IYK(imain,J),1:nx,S(imain,2)))

C Kit = Ms*Fim - Mi*FFF
C Kst = Ms*Fsm + Mi*Fim 
  	  DO 155 I = 1,nx
	  DO 155 J = 1,iny	   
	  Kit(imain,I,J) = SUM(Ms(I,1:iny)*Fim(1:iny,J))
     +                 - SUM(Mi(I,1:iny)*FFF(1:iny,J)) 
155	  Kst(imain,I,J) = SUM(Ms(I,1:iny)*Fsm(1:iny,J))
     +                 + SUM(Mi(I,1:iny)*Fim(1:iny,J)) 

C ----------------------------------
C Predictions X(t+1|t) and P(t+1|t)
C ----------------------------------
	  IF (imain.LT.d(1)) THEN
C Pi+1 = F*Pi*F'-Ci

C Ci = Mi*Fim*Mi'
	   DO 164 I = 1,nx
	   DO 164 J = 1,iny	   
164	   RG(I,J) = SUM(Mi(I,1:iny)*Fim(1:iny,J))  ! Mi*Fim
	   
	   DO 166 I = 1,nx
	   Ci(I,I) = SUM(RG(I,1:iny)*Mi(I,1:iny))
	   DO 166 J = 1,I-1
166	   Ci(I,J) = SUM(RG(I,1:iny)*Mi(J,1:iny))

	   DO 168 I = 1,nx
	   Pi(imain+1,I,I)=SUM(PFP(I,1:nx)*F(I,1:nx,S(imain+1,5)))-Ci(I,I) 
	   DO 168 J = 1,I-1	   
	   Pi(imain+1,I,J)=SUM(PFP(I,1:nx)*F(J,1:nx,S(imain+1,5)))-Ci(I,J) 
168      Pi(imain+1,J,I) = Pi(imain+1,I,J)

        ENDIF

C X(t+1|t) = a + F*X(t|t)
	  DO 170 I=1,nx
170	  XT(imain+1,I) = a(I,S(imain+1,4))
     +                + SUM(F(I,1:nx,S(imain+1,5))*X0(1:nx))

C P(t+1|t) = F*PddF' + R*R'
	  DO 172 I = 1,nx
	  DO 172 J = 1,nx	   
172	  APPO(I,J) = SUM(F(I,:,S(imain+1,5))*P0(:,J)) ! F*Pdd

	  DO 180 I = 1,nx
	   PT(imain+1,I,I) = SUM(APPO(I,1:nx)*F(I,1:nx,S(imain+1,5))) 
     +           + SUM(R(I,1:nu,S(imain+1,6))*R(I,1:nu,S(imain+1,6)))
	   DO 180 J = 1,I-1	   
	   PT(imain+1,I,J) = SUM(APPO(I,1:nx)*F(J,1:nx,S(imain+1,5))) 
     +           + SUM(R(I,1:nu,S(imain+1,6))*R(J,1:nu,S(imain+1,6)))
180	   PT(imain+1,J,I) = PT(imain+1,I,J)  

200    CONTINUE
	ENDIF

	DO 400 imain = d(1)+1,nobs
	 iny = IYK(imain,ny+1)
C -------------------------------
C Innovations: INN = yk-H*X1-c*z
C -------------------------------
       DO 210 I=1,iny
210    INN(imain,I) = yk(imain,IYK(imain,I)) 
     +       - SUM(H(IYK(imain,I),1:nx,S(imain,2))*XT(imain,:))
     +       - SUM(c(IYK(imain,I),1:nz,S(imain,1))*yk(imain,ny+1:ny+nz)) 

C ----------------------------------------------------------
C Innovation variance V = H*P1*H' + G*G' + H*R*G' + G*R'*H'
C ----------------------------------------------------------
       DO 220 I=1,iny
       DO 220 J=1,nx
220    HP1(I,J) = SUM(H(IYK(imain,I),1:nx,S(imain,2))*PT(imain,1:nx,J))
 
	 DO 221 I=1,nx
	 DO 221 J=1,iny
221	 RG(I,J)=SUM(R(I,1:nu,S(imain,6))*G(IYK(imain,J),1:nu,S(imain,3)))  	  

	 DO 222 I=1,iny
	 DO 222 J=1,iny
222	 COM(I,J)=SUM(H(IYK(imain,I),1:nx,S(imain,2))*RG(1:nx,J)) ! H*R*G'	 

	 DO 230 I=1,iny
        V(I,I) = SUM(HP1(I,1:nx)*H(IYK(imain,I),1:nx,S(imain,2))) 
     #         + SUM(G(IYK(imain,I),1:nu,S(imain,3))*
     #           G(IYK(imain,I),1:nu,S(imain,3))) + 2.*COM(I,I)
        DO 230 J=1,I-1
        V(I,J) = SUM(HP1(I,1:nx)*H(IYK(imain,J),1:nx,S(imain,2)))+
     #           SUM(G(IYK(imain,I),1:nu,S(imain,3))*
     #           G(IYK(imain,J),1:nu,S(imain,3)))+COM(I,J)+COM(J,I)  
230     V(J,I) = V(I,J)

C -------------------------------------------------------------------
C  Updating equations:
C  x0 = x1 + (P1*H'+R*G')*Vinv*INN 
C  p0 = p1 - (P1*H'+R*G')*Vinv*(P1*H'+R*G')' 
C -------------------------------------------------------------------
       IF (iny.GT.0) THEN
	  COM(1:iny,1:iny) = V(1:iny,1:iny)
	  IFAIL = -1
	  CALL F01ADF(iny,COM(1:iny+1,1:iny),iny+1,IFAIL) 	 
	  DO 240 I=1,iny
	   Vinv(imain,I,I) = COM(I+1,I)
	   DO 240 J=1,I-1
	   Vinv(imain,I,J) = COM(I+1,J)
240	   Vinv(imain,J,I) = Vinv(imain,I,J) 

	  DO 260 I=1,nx
	  DO 260 J=1,iny
260	  HPV(I,J) = SUM((HP1(1:iny,I)+RG(I,1:iny))*Vinv(imain,1:iny,J))	 

     	  DO 270 I=1,nx 	  
270	  X0(I) = XT(imain,I)+SUM(HPV(I,1:iny)*INN(imain,1:iny))

	  DO 280 I=1,nx
	  P0(I,I) = PT(imain,I,I)
     +          - SUM(HPV(I,1:iny)*(HP1(1:iny,I)+RG(I,1:iny)))
	  DO 280 J=1,I-1 	  
	  P0(I,J) = PT(imain,I,J)
     +	     - SUM(HPV(I,1:iny)*(HP1(1:iny,J)+RG(J,1:iny)))
280	  P0(J,I) = P0(I,J)
       
	 ELSE

	  X0(1:nx)      = XT(imain,1:nx)
	  P0(1:nx,1:nx) = PT(imain,1:nx,1:nx)
	  
	 ENDIF

C ------------------------------------
C Prediction       x1 = c+F*x0 
C Prediction var.  P1 = F*p0*F'+ R*R'
C ------------------------------------
	 IF (imain.LT.nobs) THEN
	  DO 290 I=1,nx
290     XT(imain+1,I) = a(I,S(imain+1,4))
     +                + SUM(F(I,1:nx,S(imain+1,5))*X0(1:nx))

	  DO 300 I=1,nx
        DO 300 J=1,nx
300     FP(I,J) = SUM(F(I,1:nx,S(imain+1,5))*P0(1:nx,J))

        DO 310 I=1,nx
        PT(imain+1,I,I) = SUM(FP(I,:)*F(I,:,S(imain+1,5))) 
     +           + SUM(R(I,1:nu,S(imain+1,6))*R(I,1:nu,S(imain+1,6)))
	  DO 310 J=1,I-1	
        PT(imain+1,I,J) = SUM(FP(I,:)*F(J,:,S(imain+1,5))) 
     +          + SUM(R(I,1:nu,S(imain+1,6))*R(J,1:nu,S(imain+1,6)))
310     PT(imain+1,J,I) = PT(imain+1,I,J) 
	 ENDIF
400   CONTINUE

C **** SMOOTHING BAKWARD RECURSIONS ****
	RECR(1:nx)      = 0.D0
	RECN(1:nx,1:nx) = 0.D0
	DO 600 ITIME = nobs,d(1)+1,-1
	 iny = IYK(ITIME,ny+1)

C R*G' and H'*Vinv
	 DO 420 J=1,iny
	 DO 420 I=1,nx
	 RG(I,J)=SUM(R(I,1:nu,S(ITIME,6))*G(IYK(ITIME,J),1:nu,S(ITIME,3))) 
420    PHFs(I,J) = 
     #  SUM(H(IYK(ITIME,1:iny),I,S(ITIME,2))*Vinv(ITIME,1:iny,J)) 

C F(t+1)*P(t|t-1)
	 DO 430 I=1,nx
       DO 430 J=1,nx
430    FP(I,J) = SUM(F(I,1:nx,S(min(nobs,ITIME+1),5))*PT(ITIME,1:nx,J))

C H'*Vinv*H
     	 DO 440 I=1,nx
       APPO(I,I) = SUM(PHFs(I,1:iny)*H(IYK(ITIME,1:iny),I,S(ITIME,2)))
	 DO 440 J=1,I-1
       APPO(I,J) = SUM(PHFs(I,1:iny)*H(IYK(ITIME,1:iny),J,S(ITIME,2)))
440    APPO(J,I) = APPO(I,J)

C L(t) =  F(t+1)-(F(t+1)*P(t|t-1)*H(t)'+R(t)*G(t)')*Vinv(t)*H(t)
    	 DO 450 I=1,nx
       DO 450 J=1,nx
450    PFP(I,J) = F(I,J,S(min(nobs,ITIME+1),5)) 
     +          - SUM(FP(I,1:nx)*APPO(1:nx,J)) 
     +          - SUM(RG(I,1:iny)*PHFs(J,1:iny)) 

C r(t-1) = H(t)'*Vinv(t)*INN(t)+L(t)'*r(t) 
	 DO 460 I=1,nx
460    WORK(I) = SUM(PFP(1:nx,I)*RECR(1:nx)) 
	 
	 DO 470 I=1,nx
470	 RECR(I) = WORK(I) + SUM(PHFs(I,1:iny)*INN(ITIME,1:iny))  

C N(t-1) = H(t)'*Vinv(t)*H(t)+L(t)'*N(t)*L(t) 
     	 DO 480 I=1,nx
       DO 480 J=1,nx
480    CC(I,J) = SUM(PFP(1:nx,I)*RECN(1:nx,J)) ! L(t)'*N(t)
     	 
	 DO 490 I=1,nx
	  RECN(I,I) = APPO(I,I) + SUM(CC(I,1:nx)*PFP(1:nx,I)) 
        DO 490 J=1,I-1
	  RECN(I,J) = APPO(I,J) + SUM(CC(I,1:nx)*PFP(1:nx,J)) 
490     RECN(J,I) = RECN(I,J)

C X(t|T) = X(t|t-1) + P(t|t-1)*r(t-1)
	 DO 500 I = 1,nx
500	 XS(ITIME,I) = XT(ITIME,I) + SUM(PT(ITIME,I,1:nx)*RECR(1:nx))

C P(t|T) = P(t|t-1) - P(t|t-1)*N(t-1)*P(t|t-1)
    	 DO 510 I=1,nx
       DO 510 J=1,nx
510    CC(I,J) = SUM(PT(ITIME,I,1:nx)*RECN(1:nx,J)) ! P(t|t-1)*N(t-1)

	 DO 520 I=1,nx
	 PS(ITIME,I,I) = PT(ITIME,I,I) - SUM(CC(I,1:nx)*PT(ITIME,1:nx,I)) 
       DO 520 J=1,I-1
	 PS(ITIME,I,J) = PT(ITIME,I,J) - SUM(CC(I,1:nx)*PT(ITIME,1:nx,J)) 
520    PS(ITIME,J,I) = PS(ITIME,I,J)

600   CONTINUE

C INITIAL KALMAN SAMOOTING
	RECRI(1:nx) = 0.D0
	DO 800 ITIME = d(1),1,-1 
	 iny = IYK(ITIME,ny+1) 
C L(t) =  F(t)-Kst(t)*H(t)     
    	 DO 610 I=1,nx
       DO 610 J=1,nx
610    PFP(I,J) = F(I,J,S(ITIME,5)) 
     +   - SUM(Kst(ITIME,I,1:iny)*H(IYK(ITIME,1:iny),J,S(ITIME,2))) 

C r(t-1)  = H(t)'*Fsm*INN(t) + L(t)'*r(t) 
C ri(t-1) = H(t)'*Fim*INN(t) + L(t)'*ri(t) + Li(t)'*r(t) 
	 DO 620 I=1,nx
       WORK(I)    = SUM(PFP(1:nx,I)*RECR(1:nx))  ! L(t)'*r(t)
620    WORK(nx+I) = SUM(PFP(1:nx,I)*RECRI(1:nx)) ! L(t)'*ri(t)

C Li(t) = -Kit(t)*H(t)     
   	 DO 621 I=1,nx
       DO 621 J=1,nx
621    PFP(I,J) = 
     #  - SUM(Kit(ITIME,I,1:iny)*H(IYK(ITIME,1:iny),J,S(ITIME,2))) 

C  L(t)'*ri(t) + Li(t)'*r(t) 
	 DO 622 I=1,nx
622    WORK(nx+I) = WORK(nx+I)+SUM(PFP(1:nx,I)*RECR(1:nx))
	 
C H'*Fsm
	 DO 625 I=1,nx
       DO 625 J=1,iny
625    PHFs(I,J) = 
     #  SUM(H(IYK(ITIME,1:iny),I,S(ITIME,2))*Vis(ITIME,1:iny,J)) 

	 DO 630 I=1,nx
630	 RECR(I) = WORK(I) + SUM(PHFs(I,1:iny)*INN(ITIME,1:iny))  

C H'*Fim
	 DO 631 I=1,nx
       DO 631 J=1,iny
631    PHFs(I,J) = 
     #  SUM(H(IYK(ITIME,1:iny),I,S(ITIME,2))*Vii(ITIME,1:iny,J)) 

	 DO 632 I=1,nx
632	 RECRI(I) = WORK(nx+I) + SUM(PHFs(I,1:iny)*INN(ITIME,1:iny)) 

C X(d|T) = X(d|d-1) + Psd*r(t-1) + Pid*ri(t-1)
	 DO 660 I = 1,nx
660	 XS(ITIME,I) = XT(ITIME,I) + SUM(PT(ITIME,I,1:nx)*RECR(1:nx))
     +             + SUM(Pi(ITIME,I,1:nx)*RECRI(1:nx))

800   CONTINUE

	DEALLOCATE(Pi,HPs,HPi,Fi,Fs,Fim,Fsm,PHFs,PHFi,FFF,Mi,Ms,Ci,
	1 Kst,Kit,W1,WORK,WORK1,PFP,APPO,APPO1,COM,RG,FP,HP1,V,CC,HPV,
     2 X0,P0,RECR,RECRI,RECN,XT,PT,INN,Vinv,Vis,Vii)

	RETURN
	END

C This is the variance and must be completed!!
C N(t-1) = H(t)'*Vinv(t)*H(t)+L(t)'*N(t)*L(t) 
c     	 DO 640 I=1,nx
c       DO 640 J=1,nx
c640    CC(I,J) = SUM(PFP(1:nx,I)*RECN(1:nx,J)) ! L(t)'*N(t)
c	 DO 650 I=1,nx
c	  RECN(I,I) = APPO(I,I) + SUM(CC(I,1:nx)*PFP(1:nx,I)) 
c        DO 650 J=1,I-1
c	  RECN(I,J) = APPO(I,J) + SUM(CC(I,1:nx)*PFP(1:nx,J)) 
c650     RECN(J,I) = RECN(I,J)
C P(t|T) = P(t|t-1) - P(t|t-1)*N(t-1)*P(t|t-1)
c    	 DO 670 I=1,nx
c       DO 670 J=1,nx
c670    CC(I,J) = SUM(PT(ITIME,I,1:nx)*RECN(1:nx,J)) ! P(t|t-1)*N(t-1)
c
c	 DO 680 I=1,nx
c	  PS(ITIME,I,I) = PT(ITIME,I,I) - SUM(CC(I,1:nx)*PT(ITIME,1:nx,I)) 
c       DO 680 J=1,I-1
c	  PS(ITIME,I,J) = PT(ITIME,I,J) - SUM(CC(I,1:nx)*PT(ITIME,1:nx,J)) 
c680    PS(ITIME,J,I) = PS(ITIME,I,J)


C -------------------------------------------------------------
C AMH DRAWS THE DISCRETE LATENT VARIABLE IN BLOCKS 
C (see Fiorentini, Planas and Rossi, Statistics and Computing, 2014, 24, 77-89) 
C Developed by A.Rossi, C.Planas and G.Fiorentini      
C Pr[Z(t)|Z(\t),Y] pto Pr[Y^(t+1,T)|Y^t,Z] x Pr[Y(t)|Y^(t-1),Z^t]
C                   x  Pr[Z(t)|Z(\t)]
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
C FURTHER INPUT 
C   nobs: # of observatios
C   d(1): order of integration of the system
C   d(2): number of non-stationary elements  
C     nv: # of discrete latent variables (S1,S2,...)              
C     ns: ns1,ns2,...
C  nstot: total # of states (states of S1 x S2 x ...x Snv)
C     nt: dimension of theta
C     np: dimension of psi
C   PMAT: (nstot x nstot) one-step transition probabilities  
C     PE: ergodic distribution of S1 x S2 x ...x Snv
C  INFOS: (9 x 6) set latent variables 
C  nstot: total # of states i.e. ns1 x ns2 x ...x nsv   
C    PTR: transition porbabilities
C     PM: marginal probabilities
C
C OUTPUT     
C       Z:(nobs x 1) takes values {1,2,...,nstot}
C ACCRATE: CUMULATES ACCEPTANS 1 IF DRAW IS ACCEPTED 0 Otherwise 
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
C -------------------------------------------------------------
	SUBROUTINE AMH(HFIX,nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np,
	1               yk,IYK,theta,psi,PTR,PM,INFOS,pdll,Z,S,ACCRATE)

	USE dfwin
	INTERFACE
	 SUBROUTINE DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R) 
	 INTEGER ny,nz,nx,nu,ns(6),nt
	 DOUBLE PRECISION theta(nt)
	 DOUBLE PRECISION c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6))
	 END SUBROUTINE
	END INTERFACE
	POINTER (pdll,fittizia) 
	POINTER (pdesign,DESIGN) 

C INPUT
	INTEGER HFIX,nobs,d(2),ny,nz,nx,nu,nv,ns(6),nstot,nt,np,
	1 IYK(nobs,ny+1),INFOS(9,6)
	DOUBLE PRECISION yk(nobs,ny+nz),theta(nt),psi(np),
	1 PTR(nobs,nstot,nstot),PM(nobs,nstot)

C INPUT/OUTPUT
	INTEGER Z(nobs),S(nobs,6),ACCRATE(nobs)

C LOCALS 
	INTEGER IT,I,II,J,JJ,K,IFAIL,ISEQ,iny,HFIXL,I0,I1,IACC,IC,IR,id1
	INTEGER Z0(nobs),IS(6),IND(HFIX+2,6),SEQ(nv),dn(2)
	DOUBLE PRECISION c(ny,max(nz,1),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6))
	DOUBLE PRECISION P1(INFOS(8,1),INFOS(8,1)),
	1 P2(INFOS(8,2),INFOS(8,2)),P3(INFOS(8,3),INFOS(8,3)),
     2 P4(INFOS(8,4),INFOS(8,4)),P5(INFOS(8,5),INFOS(8,5)),
     3 P6(INFOS(8,6),INFOS(8,6)),PMAT(nstot,nstot),PE(nstot),
	4 PALL(nstot,nstot,HFIX+1)
	DOUBLE PRECISION Xdd(max(d(1),1),nx),Pdd(max(d(1),1),nx,nx)	
	DOUBLE PRECISION OM(nobs,nx,nx),MU(nobs,nx)
	DOUBLE PRECISION HRG(ny,nu),VV(ny,ny),HF(ny,nx),COM(ny+1,ny),
	1 HRGV(nu,ny),B(nx,ny),RR(nx,nx),BH(nx,nx),BHRG(nx,nu),DD(nx,nx),
	2 CS(nx,nx),CC(nx,nx),AA(nx,nx),DI(nx,nx),COM1(nx+1,nx),Ha(ny),
	3 OMC(nx,nx),OMCDIC(nx,nx),AOMCDIC(nx,nx),AOMCDICOM(nx,nx),
	4 VVHF(ny,nx),WORK(3*nx),LAM(nx),PR(nstot),
     5 QN,QO,PA,PN,PO,GN,GO,AUX,delta,v,vd,PS0,PS1,AUX1,AUX0
	DOUBLE PRECISION, ALLOCATABLE:: DLL(:),XT(:,:),PT(:,:,:),
	1 XT0(:,:),PT0(:,:,:),PTR2(:,:,:)
	DOUBLE PRECISION EPS,ONE,ZERO
	DATA EPS/1.D-14/,ONE/1.0D0/,ZERO/0.0D0/
	DOUBLE PRECISION G05CAF,LEMMA4,MARKOVP
	
	dn(1) = 0
	dn(2) = d(2)
	id1   = max(1,d(1)) 
	Z0    = Z
	delta = 1.D-3
	pdesign = getprocaddress(pdll, "design_"C)
	CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
	CALL DESIGNZ(nv,np,psi,INFOS,P1,P2,P3,P4,P5,P6) 
C PALL(i,j) = Pr[Z(t+1)=i|Z(t)=j], Z = S1 x S2 x ... x Snv  
	CALL PPROD(nv,nstot,INFOS,P1,P2,P3,P4,P5,P6,PMAT)
C ERGODIC solves PE: PE*(I-P') = 0
      CALL ERGODIC(nstot,PMAT,PE)	 
C KOLMOGOROV EQNS
	PALL(:,:,1) = PMAT
	DO J = 2,HFIX+1  ! Kolmogorov eqns
	 DO jj = 1,nstot
	  DO ii = 1,nstot-1	 
	   PALL(ii,jj,J) = SUM(PALL(ii,:,1)*PALL(:,jj,J-1))
	  ENDDO
        PALL(nstot,jj,J) = 1.D0-SUM(PALL(1:nstot-1,jj,J))
	 ENDDO
	ENDDO

C  OMEGA and MU RECURSIONS
      OM(:,:,:)= ZERO
      MU(:,:)  = ZERO
	DO 250 IT = nobs-1,1,-1
C INT2SEQ map from Z(IT+1) to IS = (k1 k2 k3 k4 k5 k6)
	 CALL INT2SEQ(Z(IT+1),nv,INFOS,SEQ,IS) 
	 iny = IYK(IT+1,ny+1)
 
	 DO 10 I=1,iny 
	 Ha(I) = SUM(H(IYK(IT+1,I),:,IS(2))*a(:,IS(4))) ! H*a (ny x 1)
	 DO 10 J=1,nu
10	 HRG(I,J) = SUM(H(IYK(IT+1,I),:,IS(2))*R(:,J,IS(6)))
     +          + G(IYK(IT+1,I),J,IS(3))  ! HR+G (ny x nu)

	 DO 20 I=1,iny  
	 DO 20 J=1,iny
20	 VV(I,J) = SUM(HRG(I,1:nu)*HRG(J,1:nu)) ! (HR+G)*(HR+G)' (ny x ny) 

	 DO 30 I=1,iny  
	 DO 30 J=1,nx
30	 HF(I,J)=SUM(H(IYK(IT+1,I),:,IS(2))*F(:,J,IS(5))) ! HF(ny x nx)

	 COM(1:iny,1:iny) = VV(1:iny,1:iny)
	 IFAIL = -1
	 CALL F01ADF(iny,COM(1:iny+1,1:iny), iny+1, IFAIL) 
	 DO 40 I=1,iny
	 DO 40 J=1,I
	 VV(I,J) = COM(I+1,J)
40	 VV(J,I) = VV(I,J) ! inv[(HR+G)*(HR+G)'] (ny x ny) 

C B = R*(H*R+G)'*VV (nx x ny)
       DO 50 I=1,nu
	 DO 50 J=1,iny
50	 HRGV(I,J) = SUM(HRG(1:iny,I)*VV(1:iny,J))  ! (H*R+G)'*VV (nu x ny)
    
       DO 60 I=1,nx
	 DO 60 J=1,iny
60	 B(I,J) = SUM(R(I,1:nu,IS(6))*HRGV(1:nu,J)) ! B (nx x ny)

       DO 70 I=1,nx
	 DO 70 J=1,nx
	 RR(I,J) = SUM(R(I,:,IS(6))*R(J,:,IS(6)))             ! RR' (nx x nx)
70	 BH(I,J) = SUM(B(I,1:iny)*H(IYK(IT+1,1:iny),J,IS(2))) ! BH  (nx x nx)

C FIND CS such that CS*CS' = RR'-B*HRG*R' (nx x nx)
       DO 80 I=1,nx
	 DO 80 J=1,nu
80	 BHRG(I,J) = SUM(B(I,1:iny)*HRG(1:iny,J))

       DO 90 I=1,nx
	 DO 90 J=1,I
90	 CC(I,J) = RR(I,J) - SUM(BHRG(I,1:nu)*R(J,1:nu,IS(6)))
	 
       IFAIL=-1
       CALL F02FAF('V','L',nx,CC,nx,LAM,WORK,3*nx,IFAIL)
	 DO 100 I=1,nx
	 IF (LAM(I).LE.EPS) LAM(I)= ZERO
100	 CS(:,I) = CC(1:nx,I)*DSQRT(LAM(I))

C AA = F - B*HF (nx x nx)
	 DO 110 I=1,nx
	 DO 110 J=1,nx
110    AA(I,J) = F(I,J,IS(5)) - SUM(B(I,1:iny)*HF(1:iny,J))

C OMC = OM(+1)*CS (nx x nx)
	 DO 120 I=1,nx
       DO 120 J=1,nx
120	 OMC(I,J) = SUM(OM(IT+1,I,:)*CS(:,J))

C DD = I + CS'*OM(+1)*CS (nx x nx)  
       DD(:,:) = ZERO
	 DO 130 I=1,nx
       DD(I,I) = ONE
       DO 130 J=1,I
	 DD(I,J) = DD(I,J) + SUM(CS(:,I)*OMC(:,J))
130    DD(J,I) = DD(I,J)
       
C DI = inv(DD) (nx x nx) 	  
	 COM1(1:nx,:) = DD(:,:)
	 IFAIL = -1
	 CALL F01ADF(nx,COM1,nx+1,IFAIL) 
	 DO 135 I=1,nx
	 DO 135 J=1,I
	 DI(I,J) = COM1(I+1,j)
135	 DI(J,I) = DI(I,J)          

C OMCDIC = I - OM(+1)*CS*DI*CS' (nx x nx)
	 DO 140 I=1,nx
       DO 140 J=1,nx
140    COM1(I,J) = SUM(OMC(I,:)*DI(:,J))  ! OM(+1)*CS*DI (nx x nx)

       OMCDIC(:,:) = ZERO
	 DO 145 I=1,nx
	 OMCDIC(I,I) = ONE
       DO 145 J=1,nx
145    OMCDIC(I,J) = OMCDIC(I,J)-SUM(COM1(I,:)*CS(J,:)) 

C AOMCDIC = AA'*(I - OM(+1)*CS*DINV CS') (nx x nx)
    	 DO 150 I=1,nx
       DO 150 J=1,nx
150    AOMCDIC(I,J) = SUM(AA(:,I)*OMCDIC(:,J))  

C AOMCDICOM = AA'*(I - OM(+1)*CS*DINV*CS')*OM(+1) (nx x nx)
    	 DO 160 I=1,nx
       DO 160 J=1,nx
160    AOMCDICOM(I,J) = SUM(AOMCDIC(I,:)*OM(IT+1,:,J))  

C VV*H*F (ny x nx)
       DO 170 I=1,iny
	 DO 170 J=1,nx
170	 VVHF(I,J) = SUM(VV(I,1:iny)*HF(1:iny,J)) 

C OM = AA*(I - OM(+1)*C*DI*C')*OM(+1)*AA' +
C    + F'*H'*VV*H*F
    	 DO 180 I=1,nx   
       DO 180 J=1,nx   
180    OM(IT,I,J) = SUM(AOMCDICOM(I,:)*AA(:,J))
     +            + SUM(HF(1:iny,I)*VVHF(1:iny,J))  

C  MU = AA'*(I - OM(+1)*C*DI* C')*MU(+1) +
C     - AA'*(I - OM C DINV C')*OM(+1)*LAM 
C     + F'*H'*VV*(y(+1) - H*a - c*z)
C LAM = a - B*H*a + B*[y(+1)-c*z] (nx x 1)
       COM(1:iny,1) = 0.D0
	 DO 185 I=1,iny
185	 COM(I,1) = SUM(c(IYK(IT+1,I),1:nz,IS(1))*yk(IT+1,ny+1:ny+nz)) 

	 DO 190 I=1,nx
190	 LAM(I) = a(I,IS(4)) - SUM(BH(I,1:nx)*a(1:nx,IS(4)))
     +        + SUM(B(I,1:iny)*(yk(IT+1,IYK(IT+1,1:iny))
     +        - COM(1:iny,1)))  
	 DO 200 I=1,nx
200	 MU(IT,I) = SUM(AOMCDIC(I,:)*MU(IT+1,:)) 
     +          - SUM(AOMCDICOM(I,:)*LAM(:))
     +          + SUM(VVHF(1:iny,I)*(yk(IT+1,IYK(IT+1,1:iny))
     #          - Ha(1:iny)-COM(1:iny,1))) 
          	 
250	CONTINUE

      ALLOCATE (XT(0:HFIX,nx),PT(0:HFIX,nx,nx),XT0(0:HFIX,nx),
	1 PT0(0:HFIX,nx,nx),PTR2(HFIX+1,nstot,nstot),DLL(HFIX))

C First block
	I0  = 1
	I1  = HFIX
	GN  = 1.D0
	GO  = 1.D0
	QN  = 1.D0
	QO  = 1.D0
	vd = G05CAF(vd)  ! For qn(z) = delta*g0(z) + (1-delta)*gn(z)
	DO 300 IT = I1,I0,-1
	 PR(:) = PTR(IT+1,Z(IT+1),:)*PM(IT,:)/PM(IT+1,Z(IT+1))  ! P[Z(j)|Z(j+1)]
	 IF (vd.GT.delta) THEN  ! sample from gn(x)
	  v = G05CAF(v)       
	  AUX = PR(1) 
	  ISEQ = 1
	  DO 290 WHILE (AUX.LT.v) 
	   ISEQ = ISEQ + 1 
290	   AUX  = AUX  + PR(ISEQ)
       ELSE
	  v = G05CAF(v)       
	  AUX = PM(IT,1) 
	  ISEQ = 1
	  DO 291 WHILE (AUX.LT.v) 
	   ISEQ = ISEQ+1 
291	   AUX  = AUX + PM(IT,ISEQ) 
	 ENDIF             
	 Z(IT) = ISEQ
	 GN = GN*PM(IT,ISEQ)
	 GO = GO*PM(IT,Z0(IT))
	 QN = QN*PR(ISEQ)
300	 QO = QO*PTR(IT+1,Z0(IT+1),Z0(IT))*PM(IT,Z0(IT))/PM(IT+1,Z0(IT+1)) ! P[Z0(j)|Z0(j+1)]	   
       	
	QN = delta*GN + (1.D0-delta)*QN
	QO = delta*GO + (1.D0-delta)*QO

	IF (SUM(ABS(Z(I0:I1)-Z0(I0:I1))).NE.0) THEN
	 IND(1,1) = 0
	 IND(2:HFIX+2,1) = Z0(1:HFIX+1) 
       PS0 = MARKOVP(PALL,PE,nstot,HFIX,1,nobs,IND(1:HFIX+2,1))	       
	 IND(2:HFIX+2,1) = Z(1:HFIX+1)
       PS1 = MARKOVP(PALL,PE,nstot,HFIX,1,nobs,IND(1:HFIX+2,1))

       DO 305 I = 1,max(d(1),HFIX)
305	 CALL INT2SEQ(Z(I),nv,INFOS,SEQ,IND(I,:)) 
	 
	 CALL IKF(d,ny,nz,nx,nu,ns,IND(1:id1,:),yk(1:id1,:),
	1          IYK(1:id1,:),c,H,G,a,F,R,Xdd,Pdd,DLL(1:d(1)))  
       XT(d(1),1:nx)      = Xdd(id1,1:nx)   
	 PT(d(1),1:nx,1:nx) = Pdd(id1,1:nx,1:nx)
	 CALL KF(HFIX,d,ny,nz,nx,nu,ns,IND(1:HFIX,:),yk(1:HFIX,:),
	1         IYK(1:HFIX,:),c,H,G,a,F,R,XT,PT,DLL(1:HFIX)) 
       AUX1=LEMMA4(OM(HFIX,:,:),MU(HFIX,:),XT(HFIX,:),PT(HFIX,:,:),nx)
       PN = AUX1 + SUM(DLL(1:HFIX)) + PS1 ! prior x likelihood

       DO 306 I = 1,max(d(1),HFIX)
306	 CALL INT2SEQ(Z0(I),nv,INFOS,SEQ,IND(I,:)) 
	 
	 CALL IKF(d,ny,nz,nx,nu,ns,IND(1:id1,:),yk(1:id1,:),
	1          IYK(1:id1,:),c,H,G,a,F,R,Xdd,Pdd,DLL(1:d(1)))  
       XT0(d(1),1:nx)      = Xdd(id1,1:nx)   
	 PT0(d(1),1:nx,1:nx) = Pdd(id1,1:nx,1:nx)
	 CALL KF(HFIX,d,ny,nz,nx,nu,ns,IND(1:HFIX,:),yk(1:HFIX,:),
	1         IYK(1:HFIX,:),c,H,G,a,F,R,XT0,PT0,DLL(1:HFIX)) 
       AUX0=LEMMA4(OM(HFIX,:,:),MU(HFIX,:),XT0(HFIX,:),PT0(HFIX,:,:),nx)
	 PO = AUX0 + SUM(DLL(1:HFIX)) + PS0 ! prior x likelihood

c	 PA = DEXP(PN-PO)*QO/QN 
	 PA = MAX(MIN(PN-PO+DLOG(QO)-DLOG(QN),0.D0),-300.D0) 
	 PA = DEXP(PA)
	ELSE
       DO 307 I = 1,max(d(1),HFIX)
307	 CALL INT2SEQ(Z(I),nv,INFOS,SEQ,IND(I,:)) 
	 
 	 CALL IKF(d,ny,nz,nx,nu,ns,IND(1:id1,:),yk(1:id1,:),
	1          IYK(1:id1,:),c,H,G,a,F,R,Xdd,Pdd,DLL(1:d(1)))  
       XT(d(1),1:nx)      = Xdd(id1,1:nx)   
	 PT(d(1),1:nx,1:nx) = Pdd(id1,1:nx,1:nx)
	 CALL KF(HFIX,d,ny,nz,nx,nu,ns,IND(1:HFIX,:),yk(1:HFIX,:),
	1         IYK(1:HFIX,:),c,H,G,a,F,R,XT,PT,DLL(1:HFIX)) 
       XT0(HFIX,:)  = XT(HFIX,:)
	 PT0(HFIX,:,:)= PT(HFIX,:,:)
	 PA = 1.D0 
	ENDIF
	v  = G05CAF(v) 
	IF (v.GT.PA) THEN
	 IACC      = 0 
	 Z(1:HFIX) = Z0(1:HFIX)
	 ACCRATE(1:HFIX) = ACCRATE(1:HFIX) + 1
	ELSE
	 Z0(1:HFIX) = Z(1:HFIX)
	 IACC       = 1
	ENDIF

C Inner blocks
	DO 400 WHILE (I1.LT.(nobs-HFIX))
	 I1 = I1 + HFIX 
	 I0 = I1 - HFIX + 1
	 GN  = 1.D0
	 GO  = 1.D0
	 QN  = 1.D0
	 QO  = 1.D0	 
	 vd = G05CAF(vd)  ! For qn(z) = delta*g0(z) + (1-delta)*gn(z)
	 PTR2(1,:,:) = PTR(I0,:,:) ! q[Z(t+1)|Z(t)]
	 DO 310 J  = 2,HFIX+1      ! q[Z(t+j)|Z(t)] j = 2,...,HFIX+1
	 DO 310 IC = 1,nstot
	 DO 310 IR = 1,nstot
310	 PTR2(J,IR,IC) = SUM(PTR(I0+J-1,IR,:)*PTR2(J-1,:,IC))
	 DO 350 IT = I1,I0,-1  ! t+h,t+h-1,...,t+j,...,t+1
	  K = IT - I0 + 1      ! h,h-1,...,1 
 	  PR(:) = PTR(IT+1,Z(IT+1),:)*PTR2(K,:,Z(I0-1))
     #        / PTR2(K+1,Z(IT+1),Z(I0-1))
	  IF (vd.GT.delta) THEN  ! sample from gn(x)
	   v = G05CAF(v)
	   AUX = PR(1) 
	   ISEQ = 1
	   DO 320 WHILE (AUX.LT.v) 
	    ISEQ = ISEQ+1 
320	    AUX  = AUX + PR(ISEQ)
 	  ELSE
	   v = G05CAF(v)       
	   AUX = PM(IT,1) 
	   ISEQ = 1
	   DO 322 WHILE (AUX.LT.v) 
	   ISEQ = ISEQ+1 
322	   AUX  = AUX + PM(IT,ISEQ)
	  ENDIF             
	  Z(IT) = ISEQ
	  GN = GN*PM(IT,ISEQ) 
	  GO = GO*PM(IT,Z0(IT))	  
	  QN = QN*PR(ISEQ)
350	 QO = QO*PTR(IT+1,Z0(IT+1),Z0(IT))*PTR2(K,Z0(IT),Z0(I0-1))
     #    / PTR2(K+1,Z0(IT+1),Z0(I0-1))  ! P[Z0(j)|Z0(j+1)]	             

	 QN = delta*GN + (1.D0-delta)*QN
	 QO = delta*GO + (1.D0-delta)*QO

	 IF (SUM(ABS(Z(I0:I1)-Z0(I0:I1))).NE.0) THEN
	  IND(1:HFIX+2,1) = Z0(I0-1:I1+1)
        PS0 = MARKOVP(PALL,PE,nstot,HFIX,2,nobs,IND(1:HFIX+2,1))  
        IND(1:HFIX+2,1) = Z(I0-1:I1+1)
        PS1 = MARKOVP(PALL,PE,nstot,HFIX,2,nobs,IND(1:HFIX+2,1))
        
	  XT(0,:)   = IACC*XT(HFIX,:)  + (1-IACC)*XT0(HFIX,:)   ! Xt|t
	  PT(0,:,:) = IACC*PT(HFIX,:,:)+ (1-IACC)*PT0(HFIX,:,:) ! Pt|t
        XT0(0,:)  = XT(0,:)
	  PT0(0,:,:)= PT(0,:,:)
	  
        DO 360 I = I0,I1
360	  CALL INT2SEQ(Z(I),nv,INFOS,SEQ,IND(I-I0+1,:)) 
	  CALL KF(HFIX,dn,ny,nz,nx,nu,ns,IND(1:HFIX,:),yk(I0:I1,:),
	1          IYK(I0:I1,:),c,H,G,a,F,R,XT,PT,DLL(1:HFIX)) 
        AUX1=LEMMA4(OM(I1,:,:),MU(I1,:),XT(HFIX,:),PT(HFIX,:,:),nx)
        PN = AUX1 + SUM(DLL(1:HFIX)) + PS1 ! prior x likelihood

        DO 361 I = I0,I1
361	  CALL INT2SEQ(Z0(I),nv,INFOS,SEQ,IND(I-I0+1,:)) 
	 
	  CALL KF(HFIX,dn,ny,nz,nx,nu,ns,IND(1:HFIX,:),yk(I0:I1,:),
	1          IYK(I0:I1,:),c,H,G,a,F,R,XT0,PT0,DLL(1:HFIX)) 
        AUX0=LEMMA4(OM(I1,:,:),MU(I1,:),XT0(HFIX,:),PT0(HFIX,:,:),nx)
        PO = AUX0 + SUM(DLL(1:HFIX)) + PS0 ! prior x likelihood
C	  PA = DEXP(PN-PO)*QO/QN 
	  PA = MAX(MIN(PN-PO+DLOG(QO)-DLOG(QN),0.D0),-300.D0) 
	  PA = DEXP(PA)
	 ELSE
        XT(0,:)   = IACC*XT(HFIX,:)  + (1-IACC)*XT0(HFIX,:)   ! Xt|t
	  PT(0,:,:) = IACC*PT(HFIX,:,:)+ (1-IACC)*PT0(HFIX,:,:) ! Pt|t

        DO 362 I = I0,I1
362	  CALL INT2SEQ(Z(I),nv,INFOS,SEQ,IND(I-I0+1,:)) 

	  CALL KF(HFIX,dn,ny,nz,nx,nu,ns,IND(1:HFIX,:),yk(I0:I1,:),
	1          IYK(I0:I1,:),c,H,G,a,F,R,XT,PT,DLL(1:HFIX)) 
	  XT0(HFIX,:)  = XT(HFIX,:)
	  PT0(HFIX,:,:)= PT(HFIX,:,:)
        PA = 1.D0 
	 ENDIF
	 v  = G05CAF(v) 
	 IF (v.GT.PA) THEN
	  Z(I0:I1)   = Z0(I0:I1)
	  ACCRATE(I0:I1) = ACCRATE(I0:I1) + 1
	  IACC = 0
	 ELSE
	  Z0(I0:I1) = Z(I0:I1)
	  IACC = 1
	 ENDIF
400	CONTINUE
	
C Last block
	I0 = I1+1
	I1 = nobs
	HFIXL = I1 - I0 + 1
	QN = 1.D0
	QO = 1.D0	 
	vd = G05CAF(vd)  ! For qn(z) = delta*g0(z) + (1-delta)*gn(z)
	DO 500 IT = HFIXL,1,-1
	 IF (vd.GT.delta) THEN  ! sample from gn(x)
	  v = G05CAF(v)
	  AUX = PTR(nobs-IT+1,1,Z(nobs-IT))   !nobs-9: nobs
	  ISEQ = 1
	  DO 450 WHILE (AUX.LT.v) 
	   ISEQ = ISEQ+1 
450	   AUX  = AUX + PTR(nobs-IT+1,ISEQ,Z(nobs-IT))
 	 ELSE
	  v = G05CAF(v)       
	  AUX = PM(nobs-IT+1,1) 
	  ISEQ = 1
	  DO 451 WHILE (AUX.LT.v) 
	   ISEQ = ISEQ+1 
451	   AUX  = AUX + PM(nobs-IT+1,ISEQ) 
	 ENDIF             
       Z(nobs-IT+1) = ISEQ
	 GN = GN*PM(nobs-IT+1,ISEQ)
	 GO = GO*PM(nobs-IT+1,Z0(nobs-IT+1))
	 QN = QN*PTR(nobs-IT+1,ISEQ,Z(nobs-IT))
500	 QO = QO*PTR(nobs-IT+1,Z0(nobs-IT+1),Z0(nobs-IT))

	QN = delta*GN + (1.D0-delta)*QN
	QO = delta*GO + (1.D0-delta)*QO

	IF (SUM(ABS(Z(I0:I1)-Z0(I0:I1))).NE.0) THEN
	 IND(HFIXL+2,1)   = 0
	 IND(1:HFIXL+1,1) = Z0(I0-1:nobs)
       PS0 = MARKOVP(PALL,PE,nstot,HFIXL,nobs,nobs,IND(1:HFIXL+2,1))  
	 IND(1:HFIXL+1,1) = Z(I0-1:nobs)
       PS1 = MARKOVP(PALL,PE,nstot,HFIXL,nobs,nobs,IND(1:HFIXL+2,1))

       XT(0,:)   = IACC*XT(HFIXL,:)  + (1-IACC)*XT0(HFIXL,:)   ! Xt|t
	 PT(0,:,:) = IACC*PT(HFIXL,:,:)+ (1-IACC)*PT0(HFIXL,:,:) ! Pt|t
       XT0(0,:)  = XT(0,:)
	 PT0(0,:,:)= PT(0,:,:)
	 
       DO 510 I = I0,I1
510	 CALL INT2SEQ(Z(I),nv,INFOS,SEQ,IND(I-I0+1,:)) 

	 CALL KF(HFIXL,dn,ny,nz,nx,nu,ns,IND(1:HFIXL,:),yk(I0:I1,:),
	1         IYK(I0:I1,:),c,H,G,a,F,R,XT,PT,DLL(1:HFIXL)) 

       PN = SUM(DLL(1:HFIXL)) + PS1 ! prior x likelihood

       DO 520 I = I0,I1
520	 CALL INT2SEQ(Z0(I),nv,INFOS,SEQ,IND(I-I0+1,:)) 

	 CALL KF(HFIXL,dn,ny,nz,nx,nu,ns,IND(1:HFIXL,:),yk(I0:I1,:),
	1         IYK(I0:I1,:),c,H,G,a,F,R,XT0,PT0,DLL(1:HFIXL)) 

       PO = SUM(DLL(1:HFIXL)) + PS0 ! prior x likelihood
C	 PA = DEXP(PN-PO)*QO/QN
	 PA = MAX(MIN(PN-PO+DLOG(QO)-DLOG(QN),0.D0),-300.D0) 
	 PA = DEXP(PA)
	ELSE
       PA = 1.D0 
	ENDIF
	
	v = G05CAF(v) 
	IF (v.GT.PA) THEN
	 Z(I0:I1)   = Z0(I0:I1)
	 ACCRATE(I0:I1) = ACCRATE(I0:I1) + 1
	ELSE
	 Z0(I0:I1) = Z(I0:I1)
	ENDIF

	DO I=1,nobs
 	 CALL INT2SEQ(Z(I),nv,INFOS,SEQ,S(I,:))
	ENDDO
	
	DEALLOCATE (XT,PT,XT0,PT0,PTR2,DLL)
	
	RETURN
	END
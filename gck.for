C ----------------------------------------------------------------------
C GCK Implements the SINGLE-MOVE Sampler of 
C Gerlach, Carter and Kohn (2000): Efficient Bayesian Inference 
C for Dynamic Mixture Models, JASA, 95,451, pp.819-28
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C
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
C FURTHER INPUT: 
C
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
C
C OUTPUT:     
C
C     Z:(nobs x 1) takes values {1,2,...,nstot}
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
C ----------------------------------------------------------------------
	SUBROUTINE GCK(nobs,d,ny,nz,nx,nu,nv,ns,nstot,nt,np,yk,IYK,
	1               theta,psi,INFOS,pdll,Z,S)

	USE dfwin
	INTERFACE
	 SUBROUTINE DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R) 
	 INTEGER ny,nz,nx,nu,ns(6),nt
	 DOUBLE PRECISION theta(nt)
	 DOUBLE PRECISION c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6))
	 END SUBROUTINE
	END INTERFACE
	POINTER (pdll,fittizia)  !  ASSOCIATE  pointer P alla DLL ad una varibile fittizia
	POINTER (pdesign,DESIGN) ! IMPORTANT associo il puntatore pdesign alla Interface definita

C INPUT
	INTEGER nobs,d(2),ny,nz,nx,nu,nv,ns(6),nstot,nt,np,IYK(nobs,ny+1),
	1 INFOS(9,6)
	DOUBLE PRECISION yk(nobs,ny+nz),theta(nt),psi(np)

C INPUT/OUTPUT
	INTEGER Z(nobs),S(nobs,6)

C LOCALS 
	INTEGER IT,I,J,K,IFAIL,ISEQ,IMAX(1),iny,NIFS,KKK
	INTEGER IS(6),SH(3),IND(max(1,d(1)),6),SEQ(nv),IFS(nstot),dc(2)
	DOUBLE PRECISION c(ny,max(nz,1),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6))
	DOUBLE PRECISION PMAT(nstot,nstot),PE(nstot),
     1 P1(MAX(1,INFOS(8,1)*MIN(INFOS(9,1),1)),
	1    MAX(1,INFOS(8,1)*MIN(INFOS(9,1),1))),
     1 P2(MAX(1,INFOS(8,2)*MIN(INFOS(9,2),1)),
	1    MAX(1,INFOS(8,2)*MIN(INFOS(9,2),1))),
     1 P3(MAX(1,INFOS(8,3)*MIN(INFOS(9,3),1)),
	1    MAX(1,INFOS(8,3)*MIN(INFOS(9,3),1))),
     1 P4(MAX(1,INFOS(8,4)*MIN(INFOS(9,4),1)),
	1    MAX(1,INFOS(8,4)*MIN(INFOS(9,4),1))),
     1 P5(MAX(1,INFOS(8,5)*MIN(INFOS(9,5),1)),
	1    MAX(1,INFOS(8,5)*MIN(INFOS(9,5),1))),
     1 P6(MAX(1,INFOS(8,6)*MIN(INFOS(9,6),1)),
	1    MAX(1,INFOS(8,6)*MIN(INFOS(9,6),1)))
	DOUBLE PRECISION OM(nobs,nx,nx),MU(nobs,nx)
	DOUBLE PRECISION HRG(ny,nu),VV(ny,ny),HF(ny,nx),COM(ny+1,ny),
	1 HRGV(nu,ny),B(nx,ny),RR(nx,nx),BH(nx,nx),BHRG(nx,nu),DD(nx,nx),
	2 CS(nx,nx),CC(nx,nx),AA(nx,nx),DI(nx,nx),COM1(nx+1,nx),Ha(ny),
	3 OMC(nx,nx),OMCDIC(nx,nx),AOMCDIC(nx,nx),AOMCDICOM(nx,nx),
	4 VVHF(ny,nx),Xdd(0:max(1,d(1)),nx),Pdd(0:max(1,d(1)),nx,nx),
     5 WORK(3*nx),LAM(nx),AUX,U
	INTEGER, ALLOCATABLE:: IRANK(:)
	DOUBLE PRECISION, ALLOCATABLE:: DLL(:),PROB(:),PROBL(:),
	1 XT(:,:),PT(:,:,:),PMUL(:,:,:)
	DOUBLE PRECISION EPS,ONE,ZERO
	DATA EPS/1.D-14/,ONE/1.0D0/,ZERO/0.0D0/
	DOUBLE PRECISION genunf,LEMMA4,MARKOVP

	pdesign = getprocaddress(pdll, "design_"C)
	CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
	CALL DESIGNZ(nv,np,psi,INFOS,P1,P2,P3,P4,P5,P6) 
C PALL(i,j) = Pr[Z(t+1)=i|Z(t)=j], Z = S1 x S2 x ... x Snv  
	CALL PPROD(nv,nstot,INFOS,P1,P2,P3,P4,P5,P6,PMAT)
C ERGODIC solves PE: PE*(I-P') = 0
      CALL ERGODIC(nstot,PMAT,PE)	     

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

C ---------------
C START SAMPLING
C ---------------
	ALLOCATE(DLL(nstot),PROB(nstot),PROBL(nstot),XT(0:nstot,nx),
	1 PT(0:nstot,nx,nx),IRANK(nstot),PMUL(nstot,nstot,2))
	PMUL(:,:,1) = PMAT(:,:) ! one-step ahead
	PMUL(:,:,2) = 0.D0      ! two-step ahead
	DO 260 I = 1,nstot
	DO 260 J = 1,nstot
	DO 260 K = 1,nstot
260	PMUL(I,J,2) = PMUL(I,J,2) + PMAT(I,K)*PMAT(K,J)
	
C FEASIBLE Z-STATES
	NIFS   = 0
	IFS(:) = 0
	DO 265 K =1,nstot
 	 IF (PE(K).GT.0.D0) THEN	
	  NIFS = NIFS + 1
        IFS(NIFS) = K
265    ENDIF
	dc(1:2) = 0    
	DO 2000 IT = 1,nobs 
	 DO 1000 KKK = 1,NIFS
	  K = IFS(KKK)
	  CALL INT2SEQ(K,nv,INFOS,SEQ,IS(:))  
	  IF ((IT.LE.d(1)).AND.(d(1).GT.0)) THEN
   	   DO 300 I = 1,d(1) 
300	   CALL INT2SEQ(Z(I),nv,INFOS,SEQ,IND(I,:))
	   IND(IT,:)= IS(:) 
	   CALL IKF(d,ny,nz,nx,nu,ns,IND(1:d(1),:),yk(1:d(1),:),
	1            IYK(1:d(1),:),c,H,G,a,F,R,Xdd(1:d(1),:),
     2            Pdd(1:d(1),:,:),DLL(1:max(1,d(1))))  
	   XT(K,:)   = Xdd(IT,:)    ! xi(t|t) 
	   PT(K,:,:) = Pdd(IT,:,:)  ! P(t|t)		    
	   DLL(:)    = ZERO         ! log likelihood
	  ELSEIF ((IT.GT.d(1)).AND.(d(1).GT.0)) THEN 
C Input XT(0) PT(0), Output XT(K),PT(K),DLL(K)
	   Xdd(0,:)   = XT(0,:)   
	   Pdd(0,:,:) = PT(0,:,:) 	   
	   CALL KF(1,dc,ny,nz,nx,nu,ns,IS,yk(IT,:),IYK(IT,:),c,H,G,a,F,R,
	1           Xdd(0:1,:),Pdd(0:1,:,:),DLL(K))
	   XT(K,:)   = Xdd(1,:)   
	   PT(K,:,:) = Pdd(1,:,:)
	  ELSEIF ((IT.EQ.1).AND.(d(1).EQ.0)) THEN 	   
	   CALL IKF(d,ny,nz,nx,nu,ns,IS(:),yk(1,:),IYK(1,:),c,H,G,a,F,R,
	1            Xdd(0,:),Pdd(0,:,:),DLL(K))  
	   CALL KF(1,dc,ny,nz,nx,nu,ns,IS(:),yk(1,:),IYK(1,:),c,H,G,a,
	1           F,R,Xdd(0:1,:),Pdd(0:1,:,:),DLL(K)) ! log likelihood
	   XT(K,:)   = Xdd(IT,:)    ! xi(t|t) 
	   PT(K,:,:) = Pdd(IT,:,:)  ! P(t|t)		    
	  ELSEIF ((IT.GT.1).AND.(d(1).EQ.0)) THEN 
C Input XT(0) PT(0), Output XT(K),PT(K),DLL(K)
	   Xdd(0,:)   = XT(0,:)   
	   Pdd(0,:,:) = PT(0,:,:)   
	   CALL KF(1,dc,ny,nz,nx,nu,ns,IS(:),yk(IT,:),IYK(IT,:),c,H,G,a,
	1           F,R,Xdd(0:1,:),Pdd(0:1,:,:),DLL(K))
	   XT(K,:)   = Xdd(1,:)   
	   PT(K,:,:) = Pdd(1,:,:)
	  ENDIF
	  
	  SH(1)  = Z(max(1,IT-1))
	  SH(2)  = K  
	  SH(3)  = Z(min(nobs,IT+1))
	  PROBL(K) = DLL(K)
     +           + LEMMA4(OM(IT,:,:),MU(IT,:),XT(K,:),PT(K,:,:),nx)
     +           + MARKOVP(PMUL,PE,nstot,1,IT,nobs,SH)
	  
1000   CONTINUE

C ---------------------------------------------
C SAMPLING Z(t:t+h-1) using PROB
C ISEQ is the sampled sequence - out of nstot
C ---------------------------------------------  
C To prevent exp overflow
	 PROB(:) = 0.D0
	 IMAX = MAXLOC(PROBL(IFS(1:NIFS)))
	 KKK = IFS(IMAX(1))
	 PROBL(IFS(1:NIFS)) = PROBL(IFS(1:NIFS))-PROBL(KKK)
	 PROB(IFS(1:NIFS))  = DEXP(PROBL(IFS(1:NIFS)))
     #                    / SUM(DEXP(PROBL(IFS(1:NIFS))))
	 
C	 U = G05CAF(U) ! Sampling from U(0,1)  
       U = genunf(0.D0,1.D0)
	 ISEQ = 1
	 AUX  = PROB(1)
	 DO 310 WHILE (AUX.LT.U) 
	 ISEQ = ISEQ + 1 
310	 AUX  = AUX  + PROB(ISEQ)

	 Z(IT) = ISEQ

	 XT(0,:)   = XT(ISEQ,:)
	 PT(0,:,:) = PT(ISEQ,:,:)

2000  CONTINUE   

	DO I=1,nobs
 	 CALL INT2SEQ(Z(I),nv,INFOS,SEQ,S(I,:))
	ENDDO
	
	DEALLOCATE(DLL,PROB,PROBL,XT,PT,IRANK,PMUL)
	RETURN
	END

C **** da butta ******************
c	DO IT = 1,nobs-1
c	 CALL INT2SEQ(Z(IT),nv,INFOS,SEQ,SS(IT,:)) 
c	ENDDO

c	SS(nobs,:) = 1 
c	CALL IKF(d,ny,nz,nx,nu,ns,SS(1:max(d(1),1),1:6),
c	1         yk(1:max(d(1),1),1:ny+nz),IYK(1:max(d(1),1),1:ny+1),
c    2         c,H,G,a,F,R,Xdd,Pdd,LIKE(1:max(d(1),1)))
c	XX(d(1),1:nx) = Xdd(max(d(1),1),1:nx)   
c	PP(d(1),1:nx,1:nx) = Pdd(max(d(1),1),1:nx,1:nx)
c	CALL KF(nobs,d,ny,nz,nx,nu,ns,SS,yk,IYK,c,H,G,a,F,R,XX,PP,LIKE)
	
c	R1 = DEXP(SUM(LIKE))*P1(SS(nobs,4),SS(nobs-1,4))
	
c	SS(nobs,4) = 2 
c	CALL IKF(d,ny,nz,nx,nu,ns,SS(1:max(d(1),1),1:6),
c	1         yk(1:max(d(1),1),1:ny+nz),IYK(1:max(d(1),1),1:ny+1),
c     2         c,H,G,a,F,R,Xdd,Pdd,LIKE(1:max(d(1),1)))
c	XX(d(1),1:nx) = Xdd(max(d(1),1),1:nx)   
c	PP(d(1),1:nx,1:nx) = Pdd(max(d(1),1),1:nx,1:nx)
c	CALL KF(nobs,d,ny,nz,nx,nu,ns,SS,yk,IYK,c,H,G,a,F,R,XX,PP,LIKE)
c	R2 = DEXP(SUM(LIKE))*P1(SS(nobs,4),SS(nobs-1,4))
c	R1 = R1/(R1+R2)

C --------------------------------------------------------------------
C INNOV2 (no mising values) RETURNS MODEL INNOVATIONS
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
      SUBROUTINE INNOV2(nobs,d,ny,nz,nx,nu,ns,nt,S,yk,
	1                  theta,pdll,INN)
#if !defined(DYNARE)
	USE dfwin
#endif
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
	INTEGER nobs,d(2),ny,nz,nx,nu,nt,ns(6),S(nobs,6)
	DOUBLE PRECISION yk(nobs,ny+nz),theta(nt)

C OUTPUT
	DOUBLE PRECISION INN(nobs,ny)

C LOCALS
	INTEGER imain,I,J,IFAIL
	DOUBLE PRECISION,ALLOCATABLE::R(:,:,:),c(:,:,:),H(:,:,:),
	1 G(:,:,:),a(:,:),F(:,:,:)
      DOUBLE PRECISION,ALLOCATABLE:: X1(:),P1(:,:),FP(:,:),HP1(:,:),
     1 V(:,:),Vinv(:,:),COM(:,:),RG(:,:),HPV(:,:),Xdd(:,:),Pdd(:,:,:),
     1 LIKE(:),XT(:),PT(:,:)

      ALLOCATE(c(ny,max(nz,1),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6)))
      ALLOCATE(X1(nx),P1(nx,nx),FP(nx,nx),HP1(ny,nx),V(ny,ny),
     1 Vinv(ny,ny),COM(ny+1,ny),RG(nx,ny),HPV(nx,ny),
     1 Xdd(max(d(1),1),nx),Pdd(max(d(1),1),nx,nx),LIKE(max(d(1),1)),
     1 XT(nx),PT(nx,nx))

      pdesign = getprocaddress(pdll, "design_"C)
	INN(:,:) = 0.D0
	CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
	CALL IKF2(d,ny,nz,nx,nu,ns,S(1:max(d(1),1),1:6),
	1 yk(1:max(d(1),1),1:ny+nz),c,H,G,a,F,R,
     2 Xdd,Pdd,LIKE(1:max(d(1),1)))
	XT(1:nx)      = Xdd(max(d(1),1),1:nx)
	PT(1:nx,1:nx) = Pdd(max(d(1),1),1:nx,1:nx)

      DO 1000 imain = d(1)+1,nobs

C ------------------------------------
C Prediction       x1 = c+F*x0
C Prediction var.  P1 = F*p0*F'+ R*R'
C ------------------------------------
	 DO 10 I=1,nx
10     X1(I) = a(I,S(imain,4))+SUM(F(I,:,S(imain,5))*XT(:))

	 DO 20 I=1,nx
       DO 20 J=1,nx
20     FP(I,J) = SUM(F(I,:,S(imain,5))*PT(:,J))

       DO 30 I=1,nx
        P1(I,I) = SUM(FP(I,:)*F(I,:,S(imain,5)))
     +          + SUM(R(I,1:nu,S(imain,6))*R(I,1:nu,S(imain,6)))
	  DO 30 J=1,I-1
        P1(I,J) = SUM(FP(I,:)*F(J,:,S(imain,5)))
     +          + SUM(R(I,1:nu,S(imain,6))*R(J,1:nu,S(imain,6)))
30      P1(J,I) = P1(I,J)

C -------------------------------
C Innovations: INN = yk-H*X1-c*z
C -------------------------------
	 DO 40 I=1,ny
40	 INN(imain,I) = yk(imain,I)
     +       - SUM(H(I,1:nx,S(imain,2))*X1(1:nx))
     +       - SUM(c(I,1:nz,S(imain,1))*yk(imain,ny+1:ny+nz))

C ----------------------------------------------------------
C Innovation variance V = H*P1*H' + G*G' + H*R*G' + G*R'*H'
C ----------------------------------------------------------
       DO 50 I=1,ny
       DO 50 J=1,nx
50     HP1(I,J) = SUM(H(I,1:nx,S(imain,2))*P1(1:nx,J))

	 DO 55 I=1,nx
	 DO 55 J=1,ny
55	 RG(I,J) = SUM(R(I,1:nu,S(imain,6))
     #         * G(J,1:nu,S(imain,3)))  ! R*G'

	 DO 56 I=1,ny
	 DO 56 J=1,ny
56	 COM(I,J)=SUM(H(I,1:nx,S(imain,2))*RG(1:nx,J)) ! H*R*G'

	 DO 60 I=1,ny
        V(I,I) = SUM(HP1(I,1:nx)*H(I,1:nx,S(imain,2))) +
     #           SUM(G(I,1:nu,S(imain,3))*
     #           G(I,1:nu,S(imain,3)))+2.*COM(I,I)
        DO 60 J=1,I-1
        V(I,J) = SUM(HP1(I,1:nx)*H(J,1:nx,S(imain,2)))+
     #           SUM(G(I,1:nu,S(imain,3))*
     #           G(J,1:nu,S(imain,3)))+COM(I,J)+COM(J,I)
60      V(J,I) = V(I,J)

C -------------------------------------------------------------------
C  Updating equations:
C  x0 = x1 + (P1*H'+R*G')*Vinv*INN
C  p0 = p1 - (P1*H'+R*G')*Vinv*(P1*H'+R*G')'
C -------------------------------------------------------------------
       IF (ny.GT.0) THEN
	  COM(1:ny,1:ny) = V(1:ny,1:ny)
	  IFAIL = -1
c	  CALL F01ADF(ny,COM(1:ny+1,1:ny),ny+1,IFAIL)
        CALL DPOTRF('L',ny,COM(1:ny,1:ny),ny,IFAIL) ! COM = L*L'
        CALL DPOTRI('L',ny,COM(1:ny,1:ny),ny,IFAIL) ! COM = VV^-1

	  DO 70 I=1,ny
	   Vinv(I,I) = COM(I,I)
	   DO 70 J=1,I-1
	   Vinv(I,J) = COM(I,J)
70	   Vinv(J,I) = Vinv(I,J)

	  DO 90 I=1,nx
	  DO 90 J=1,ny
90      HPV(I,J) = SUM((HP1(1:ny,I)+RG(I,1:ny))*Vinv(1:ny,J))

     	  DO 100 I=1,nx
100	  XT(I) = X1(I)+SUM(HPV(I,1:ny)*INN(imain,1:ny))

	  DO 110 I=1,nx
	   PT(I,I) = P1(I,I)
     +	             - SUM(HPV(I,1:ny)*(HP1(1:ny,I)+RG(I,1:ny)))
	  DO 110 J=1,I-1
	  PT(I,J) = P1(I,J)
     +                - SUM(HPV(I,1:ny)*(HP1(1:ny,J)+RG(J,1:ny)))
110	  PT(J,I) = PT(I,J)

	 ELSE

	  XT(1:nx)      = X1(1:nx)
	  PT(1:nx,1:nx) = P1(1:nx,1:nx)

1000   ENDIF

      DEALLOCATE(R,c,H,G,a,F,X1,P1,FP,HP1,V,Vinv,COM,RG,HPV,Xdd,Pdd,
     1 LIKE,XT,PT)
	RETURN
	END

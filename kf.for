C --------------------------------------------------------------------
C KF implements the Kalman filter recursions and returns
C XT = E[x(t)|t], PT = V[x(t)|t]
C LIKE(t) = log-likelihood(t), t = 1,...,nobs   
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
      SUBROUTINE KF(nobs,d,ny,nz,nx,nu,ns,S,yk,IYK,c,H,G,a,F,R,
	1              XT,PT,LIKE)  
C INPUT
	INTEGER nobs,d(2),ny,nz,nx,nu,ns(6),S(nobs,6),IYK(nobs,ny+1)
	DOUBLE PRECISION yk(nobs,ny+nz),c(ny,max(nz,1),ns(1)),
	1 H(ny,nx,ns(2)),G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),
     2 R(nx,nu,ns(6))

C OUTPUT
	DOUBLE PRECISION LIKE(nobs)
	DOUBLE PRECISION XT(0:nobs,nx),PT(0:nobs,nx,nx)	

C LOCALS
	INTEGER imain,iny,I,J,IFAIL
	DOUBLE PRECISION RSS,DETV,X1(nx),P1(nx,nx),FP(nx,nx),
	1 INN(ny),HP1(ny,nx),V(ny,ny),Vinv(ny,ny),COM(ny+1,ny),
     2 RG(nx,ny),HPV(nx,ny)	  
	 
      LIKE(d(1)+1:nobs) = 0.D0
      DO 1000 imain = d(1)+1,nobs
	 iny = IYK(imain,ny+1)
	   
C ------------------------------------
C Prediction       x1 = c+F*x0 
C Prediction var.  P1 = F*p0*F'+ R*R'
C ------------------------------------
	 DO 10 I=1,nx
10     X1(I) = a(I,S(imain,4))+SUM(F(I,:,S(imain,5))*XT(imain-1,:))

	 DO 20 I=1,nx
       DO 20 J=1,nx
20     FP(I,J) = SUM(F(I,:,S(imain,5))*PT(imain-1,:,J))

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
	 DO 40 I=1,iny
40	 INN(I)= yk(imain,IYK(imain,I))
     +       - SUM(H(IYK(imain,I),1:nx,S(imain,2))*X1(1:nx))
     +       - SUM(c(IYK(imain,I),1:nz,S(imain,1))*yk(imain,ny+1:ny+nz))

C ---------------------------------------------------------
C Innovation variance V = H*P1*H' + G*G' + H*R*G' + G*R'*H'
C ---------------------------------------------------------
       DO 50 I=1,iny
       DO 50 J=1,nx
50     HP1(I,J) = SUM(H(IYK(imain,I),1:nx,S(imain,2))*P1(1:nx,J))

	 DO 55 I=1,nx
	 DO 55 J=1,iny
55	 RG(I,J) = SUM(R(I,1:nu,S(imain,6)) 
     #         * G(IYK(imain,J),1:nu,S(imain,3)))  ! R*G'	 

	 DO 56 I=1,iny
	 DO 56 J=1,iny
56	 COM(I,J)=SUM(H(IYK(imain,I),1:nx,S(imain,2))*RG(1:nx,J)) ! H*R*G'	 

	 DO 60 I=1,iny
        V(I,I) = SUM(HP1(I,1:nx)*H(IYK(imain,I),1:nx,S(imain,2))) +
     #           SUM(G(IYK(imain,I),1:nu,S(imain,3))*
     #           G(IYK(imain,I),1:nu,S(imain,3))) + 2.*COM(I,I)
     
        DO 60 J=1,I-1
        V(I,J) = SUM(HP1(I,1:nx)*H(IYK(imain,J),1:nx,S(imain,2)))+
     #           SUM(G(IYK(imain,I),1:nu,S(imain,3))*
     #           G(IYK(imain,J),1:nu,S(imain,3)))+COM(I,J)+COM(J,I)
60      V(J,I) = V(I,J)

C -------------------------------------------------------------------
C  Updating equations:
C  x0 = x1 + (P1*H'+R*G')*Vinv*INN 
C  p0 = p1 - (P1*H'+R*G')*Vinv*(P1*H'+R*G')' 
C -------------------------------------------------------------------
       IF (iny.GT.0) THEN
	  COM(1:iny,1:iny) = V(1:iny,1:iny)
	  IFAIL = -1
	  CALL F01ADF(iny,COM(1:iny+1,1:iny),iny+1,IFAIL)     
	  DO 70 I=1,iny
	   Vinv(I,I) = COM(I+1,I)
	   DO 70 J=1,I-1
	   Vinv(I,J) = COM(I+1,J)
70	   Vinv(J,I) = Vinv(I,J) 

	  DO 90 I=1,nx
	  DO 90 J=1,iny	 
90      HPV(I,J) = SUM((HP1(1:iny,I)+RG(I,1:iny))*Vinv(1:iny,J))

     	  DO 100 I=1,nx 	  
100	  XT(imain,I) = X1(I)+SUM(HPV(I,1:iny)*INN(1:iny))

	  DO 110 I=1,nx
	   PT(imain,I,I) = P1(I,I)
     +	             - SUM(HPV(I,1:iny)*(HP1(1:iny,I)+RG(I,1:iny)))
	  DO 110 J=1,I-1 	  
	  PT(imain,I,J) = P1(I,J)
     +                - SUM(HPV(I,1:iny)*(HP1(1:iny,J)+RG(J,1:iny)))
110	  PT(imain,J,I) = PT(imain,I,J)

C ---------------------------------------------                 
C  Log-Likelihood = -(RSS + ln(det(V))/2
C	        RSS = INN'*Vinv*INN
C ---------------------------------------------
	  IFAIL=-1
        CALL F03ABF(V(1:iny,1:iny),iny,iny,DETV,COM(1:iny,1),IFAIL)   
       
	  RSS = 0.D0
	  DO 120 I=1,iny
	  DO 120 J=1,iny
120     RSS = RSS + INN(I)*Vinv(I,J)*INN(J)

	 LIKE(imain) = -.5D0*(RSS + DLOG(DETV)) 
     #             - iny/2.D0*DLOG(2.*3.141592653589793D0)
	 ELSE 
	 
	  XT(imain,1:nx)      = X1(1:nx)
	  PT(imain,1:nx,1:nx) = P1(1:nx,1:nx)

1000	 ENDIF 

	RETURN
	END
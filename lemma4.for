C ------------------------------------------------------------
C LEMMA4 returns log(f(y(t+1),...,y(T)|y(1),...,y(t),theta,S))
C For details see lemma 4 in Gerlach et al. JASA, 2000
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
C ------------------------------------------------------------
	DOUBLE PRECISION FUNCTION LEMMA4(OM,MU,X0,P0,nx)  

! INPUT
	INTEGER nx
	DOUBLE PRECISION OM(nx,nx),MU(nx),X0(nx),P0(nx,nx)

! LOCALS
	INTEGER I,J,NULLITY,IFAIL
      DOUBLE PRECISION COM(nx+1,nx),OMC(nx,nx),C(nx,nx),T(nx,nx),
	1 W(nx),LAM(nx),WORK(3*nx),DETV,AUX,EPS

	EPS   = 1.D-14
c	IFAIL = -1
c      COM(1:nx,1:nx) = P0(:,:)
c	CALL F02FAF('V','L',nx,COM(1:nx,1:nx),nx,LAM,WORK,3*nx,IFAIL) 
c	DO 10 I=1,nx
c 	IF (LAM(I).LE.EPS) LAM(I) = 0.D0
c10	T(:,I) = COM(1:nx,I)*DSQRT(LAM(I))  
	T(:,:) = 0.D0 
	CALL SCHOLLU(P0,T,nx,NULLITY,EPS,IFAIL)
     
C  OMC = T' OM T + I 
      DO 20 I=1,nx
      DO 20 J=1,nx
20	COM(I,J) = SUM(T(:,I)*OM(:,J))

      OMC(:,:) = 0.D0
      DO 30 I=1,nx
      OMC(I,I) = 1.D0
      DO 30 J=1,I
	OMC(I,J) = SUM(COM(I,:)*T(:,J)) + OMC(I,J)
30    OMC(J,I) = OMC(I,J)

C  OMC = inv(T' OM T + I )
	COM(1:nx,:) = OMC(:,:)
	IFAIL = -1
	CALL F01ADF(nx, COM, nx+1, IFAIL) 
	DO 40 I=1,nx
	DO 40 J=1,I
	OMC(I,J) = COM(I+1,J)
40	OMC(J,I) = OMC(I,J)           

	COM(1:nx,:) = OMC(:,:) ! CARE OMC IS INVERTED
	IFAIL=-1
	CALL F03ABF(COM(1:nx,1:nx),nx,nx,DETV,W,IFAIL)

C  AUX = -.5*m' OM m  + MU'*m +.5*(mu-om*m)'*T(T'omT+I)^-1*T'*(mu-om*m)
      AUX = SUM(MU(:)*X0(:)) ! MU'*m
      DO 50 I=1,nx	    
50	C(I,1) = SUM(X0(:)*OM(:,I))  ! OM*m
	AUX = AUX - .5D0*SUM(C(:,1)*X0(:)) ! -.5*m' OM m
	  
	DO 60 I=1,nx
60	COM(I,1) = SUM((MU(:)-C(:,1))*T(:,I))  !(mu-om*m)'*T

C  .5*(mu-om*m)'*T(T'omT+I)^-1*T'*(mu-om*m)
      DO 70 I=1,nx
      DO 70 J=1,nx
70	AUX = AUX + .5D0*COM(I,1)*COM(J,1)*OMC(I,J)

	LEMMA4 = .5D0*DLOG(DETV) + AUX  	 	 	   
      
	RETURN
	END
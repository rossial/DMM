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
C ------------------------------------------------------------
	DOUBLE PRECISION FUNCTION LEMMA4(OM,MU,X0,P0,nx)

! INPUT
	INTEGER nx
	DOUBLE PRECISION OM(nx,nx),MU(nx),X0(nx),P0(nx,nx)

! LOCALS
	INTEGER I,J,NULLITY,IFAIL
      DOUBLE PRECISION COM(nx+1,nx),OMC(nx,nx),C(nx,nx),T(nx,nx),
	1 DETV,AUX,EPS  ! W(nx),LAM(nx),WORK(3*nx)

! EXTERNAL SUBROUTINES
      EXTERNAL SCHOLLU,DPOTRF,DPOTRI

	EPS    = 1.D-14
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
C	CALL F01ADF(nx, COM, nx+1, IFAIL)
      CALL DPOTRF('L',nx,COM(1:nx,1:nx),nx,IFAIL) ! COM = L*L'
      DETV = 1.D0 ! det(L)
      DO I =1,nx
       DETV = DETV*COM(I,I)
      ENDDO
      CALL DPOTRI('L',nx,COM(1:nx,1:nx),nx,IFAIL) ! COM = VV^-1

	DO 40 I=1,nx
      OMC(I,I) = COM(I,I)
	DO 40 J=1,I-1
	OMC(I,J) = COM(I,J)
40	OMC(J,I) = OMC(I,J)

C	COM(1:nx,:) = OMC(:,:) ! CARE OMC IS INVERTED
C	IFAIL=-1
C	CALL F03ABF(COM(1:nx,1:nx),nx,nx,DETV,W,IFAIL)

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

c	LEMMA4 = .5D0*DLOG(DETV) + AUX
      LEMMA4 = -1.D0*DLOG(DETV) + AUX

	RETURN
	END

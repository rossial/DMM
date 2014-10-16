C -------------------------------------------------------------
C LYAP solves the Lyapunov equation Ps-F*Ps*F'=RR
C where RR and Ps are symmetric matrices.
C Developed by __GFORTRAN__ team
C Recoded in Fortran by A.Rossi, C.Planas and G.Fiorentini
C
C Copyright (C) 2006-2012 Dynare Team
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
C -------------------------------------------------------------
	SUBROUTINE LYAP(nx,nu,compt,F,R,Ps)
C INPUT
	INTEGER nx,nu
	DOUBLE PRECISION compt,F(nx,nx),R(nx,nu)
C OUTPUT
	DOUBLE PRECISION Ps(nx,nx)
C LOCALS
      INTEGER I,J,IFAIL,LWORK,IPIV(nx**2)
	LOGICAL BWORK(nx),SELECT
      DOUBLE PRECISION, ALLOCATABLE:: WORK(:),RR(:,:),WR(:),WI(:),
	1 Z(:,:),T(:,:),ZR(:,:),ZRZ(:,:),Q(:,:),WR1(:),Z1(:)
C EXTERNAL SUBROUTINES
      EXTERNAL DGETRF,DGETRI,DGEES,SELECT

C RR = R*R'
	ALLOCATE(RR(nx,nx))
	DO 10 I = 1,nx
	RR(I,I) = SUM(R(I,:)*R(I,:))
	DO 10 J = 1,I-1
	RR(I,J) = SUM(R(I,:)*R(J,:))
10	RR(J,I) = RR(I,J)

	IF (nx.EQ.1) THEN
	 Ps(1,1) = RR(1,1)/(1.D0-F(1,1)**2)
	 DEALLOCATE(RR)
	 GOTO 7777
	ENDIF

	LWORK = 3*nx
	ALLOCATE(T(nx,nx),WORK(LWORK),WR(nx),WI(nx),Z(nx,nx),ZRZ(nx,nx),
	1         ZR(nx,nx),Q(2*nx,2*nx),WR1(nx),Z1(2*nx))

      T(1:nx,1:nx) = F(1:nx,1:nx)
	IFAIL = -1
C	CALL F02EAF('V',nx,T,nx,WR,WI,Z,nx,WORK,LWORK,IFAIL) ! F = ZTZ'
      I = 0
      CALL DGEES('V','N',SELECT,nx,T,nx,I,WR,WI,Z,nx,WORK,LWORK,
     #           BWORK,IFAIL)
	DO I = 1,nx
	IF (WI(I)**2+WR(I)**2.GE.1.D0) THEN
#if defined(MEX)
       CALL mexErrMsgTxt('\nLYAPUNOV SUBROUTINE: Some parameters out of\nstationary region. Check hyptheta in namelist prior.\nProgram aborting')
#elif defined(__GFORTRAN__)
       WRITE(*,*) ' '
       WRITE(*,*) ' LYAPUNOV SUBROUTINE: Some parameters out of '
       WRITE(*,*) ' stationary region. Check hyptheta in namelist '
       WRITE(*,*) ' prior.'
       WRITE(*,*) ' Program aborting'
#else
	 TYPE *, ' '
	 TYPE *, ' LYAPUNOV SUBROUTINE: Some parameters out of '
	 TYPE *, ' stationary region. Check hyptheta in namelist prior.'
	 TYPE *, ' Program aborting'
#endif
	 PAUSE
	 STOP
	ENDIF
	ENDDO

C ZRZ = Z'*RR*Z  (B in Matlab)
	DO 20 I = 1,nx
	DO 20 J = 1,nx
20	ZR(I,J) = SUM(Z(:,I)*RR(:,J))  ! Z'*RR

	DO 30 I = 1,nx
	DO 30 J = 1,nx
30	ZRZ(I,J) = SUM(ZR(I,:)*Z(:,J))

	I = nx
	WR(:)   = 0.D0  !(c  in Matalab)
	WR1(:)  = 0.D0  !(c1 in Matalab)
	Ps(:,:) = 0.D0  !(x  in Matalab)
	DO WHILE (I.GE.2)
	 IF (DABS(T(I,I-1)).LT.compt) THEN
	  IF (I.LT.nx) THEN
c        c = T(1:i,:)*(x(:,i+1:end)*T(i,i+1:end)') + ...
c            T(i,i)*T(1:i,i+1:end)*x(i+1:end,i);
	   	DO 40 J = 1,nx
40		WI(J) = SUM(Ps(J,I+1:nx)*T(I,I+1:nx))
          DO 50 J = 1,I
50		WR(J) = SUM(T(J,1:nx)*WI(1:nx))
     +          + T(I,I)*SUM(T(J,I+1:nx)*Ps(I+1:nx,I))

	  ENDIF
c       q = eye(i)-T(1:i,1:i)*T(i,i);
c       x(1:i,i) = q\(B(1:i,i)+c); x = inv(q)*(B(1:i,i)+c)
c       x(i,1:i-1) = x(1:i-1,i)';
c       i = i - 1;
	  Q(1:I,1:I) = -T(I,I)*T(1:I,1:I)
	  DO 60 J = 1,I
60	  Q(J,J) = Q(J,J) + 1.D0

        IFAIL = -1
	  ZR(1:I,1:I) = Q(1:I,1:I)
c	  CALL F07ADF(I,I,ZR(1:I,1:I),I,IPIV(1:I),IFAIL)
c	  CALL F07AJF(I,ZR(1:I,1:I),I,IPIV(1:I),WORK(1:64*I),64*I,IFAIL)
        CALL DGETRF(I,I,ZR(1:I,1:I),I,IPIV(1:I),IFAIL)
	  CALL DGETRI(I,ZR(1:I,1:I),I,IPIV(1:I),WORK(1:64*I),64*I,IFAIL)
	  DO 70 J = 1,I
70	  Ps(J,I) = SUM(ZR(J,1:I)*(ZRZ(1:I,I)+WR(1:I)))
	  Ps(I,1:I-1) = Ps(1:I-1,I)
	  I = I - 1
	 ELSE
	  IF (I.LT.nx) THEN
c        c = T(1:i,:)*(x(:,i+1:end)*T(i,i+1:end)') + ...
c            T(i,i)*T(1:i,i+1:end)*x(i+1:end,i) + ...
c            T(i,i-1)*T(1:i,i+1:end)*x(i+1:end,i-1);
	   	DO 90 J = 1,nx
90		WI(J) = SUM(Ps(J,I+1:nx)*T(I,I+1:nx))
          DO 100 J = 1,I
100		WR(J) = SUM(T(J,1:nx)*WI(1:nx))
     +          + T(I,I)*SUM(T(J,I+1:nx)*Ps(I+1:nx,I))
     +          + T(I,I-1)*SUM(T(J,I+1:nx)*Ps(I+1:nx,I-1))
c         c1 = T(1:i,:)*(x(:,i+1:end)*T(i-1,i+1:end)') + ...
c              T(i-1,i-1)*T(1:i,i+1:end)*x(i+1:end,i-1) + ...
c              T(i-1,i)*T(1:i,i+1:end)*x(i+1:end,i);
	   	DO 110 J = 1,nx
110		WI(J) = SUM(Ps(J,I+1:nx)*T(I-1,I+1:nx))

          DO 120 J = 1,I
120		WR1(J) = SUM(T(J,1:nx)*WI(1:nx))
     +           + T(I-1,I-1)*SUM(T(J,I+1:nx)*Ps(I+1:nx,I-1))
     +           + T(I-1,I)*SUM(T(J,I+1:nx)*Ps(I+1:nx,I))
        ENDIF
c       q = [  eye(i)-T(1:i,1:i)*T(i,i) ,  -T(1:i,1:i)*T(i,i-1) ; ...
c            -T(1:i,1:i)*T(i-1,i)     ,   eye(i)-T(1:i,1:i)*T(i-1,i-1) ];
	  Q(1:I,1:I)         = -T(I,I)*T(1:I,1:I)
	  Q(I+1:2*I,I+1:2*I) = -T(I-1,I-1)*T(1:I,1:I)
	  Q(1:I,I+1:2*I)     = -T(I,I-1)*T(1:I,1:I)
	  Q(I+1:2*I,1:I)     = -T(I-1,I)*T(1:I,1:I)
	  DO 130 J = 1,2*I
130	  Q(J,J) = Q(J,J) + 1.D0
c       z =  q\[ B(1:i,i)+c ; B(1:i,i-1) + c1 ];
        IFAIL = -1
C	  CALL F07ADF(2*I,2*I,Q(1:2*I,1:2*I),2*I,IPIV(1:2*I),IFAIL)
C	  CALL F07AJF(2*I,Q(1:2*I,1:2*I),2*I,IPIV(1:2*I),WORK(1:64*2*I),
C     #	          64*2*I,IFAIL)
	  CALL DGETRF(2*I,2*I,Q(1:2*I,1:2*I),2*I,IPIV(1:2*I),IFAIL)
	  CALL DGETRI(2*I,Q(1:2*I,1:2*I),2*I,IPIV(1:2*I),WORK(1:64*2*I),
     # 64*2*I,IFAIL)

	  DO 140 J = 1,2*I
140	  Z1(J) = SUM(Q(J,1:I)*(ZRZ(1:I,I) + WR(1:I)))
     +        + SUM(Q(J,I+1:2*I)*(ZRZ(1:I,I-1) + WR1(1:I)))
c       x(1:i,i)     = z(1:i);
c       x(1:i,i-1)   = z(i+1:end);
c       x(i,1:i-1)   = x(1:i-1,i)';
c       x(i-1,1:i-2) = x(1:i-2,i-1)';
     	  Ps(1:I,I)     = Z1(1:i)
	  Ps(1:I,I-1)   = Z1(I+1:2*I)
	  Ps(I,1:I-1)   = Ps(1:I-1,I)
	  Ps(I-1,1:I-2) = Ps(1:I-2,I-1)
	  I = I - 2
	 ENDIF
	END DO
c if i == 1
c    c = T(1,:)*(x(:,2:end)*T(1,2:end)') + T(1,1)*T(1,2:end)*x(2:end,1);
c    x(1,1) = (B(1,1)+c)/(1-T(1,1)*T(1,1));
c end
	IF (I.EQ.1) THEN
	 DO 150 J =1,nx
	 WI(J) = SUM(Ps(J,2:nx)*T(1,2:nx))
150    WR(1) = SUM(T(1,1:nx)*WI(1:nx))
     +       + T(1,1)*SUM(T(1,2:nx)*Ps(2:nx,1))
	 Ps(1,1) = (ZRZ(1,1)+WR(1))/(1.D0-T(1,1)**2)
	ENDIF
c x = U(:,:)*x*U(:,:)';
	DO 160 I = 1,nx
	DO 160 J = 1,nx
160	ZR(I,J) = SUM(Z(I,1:nx)*Ps(1:nx,J))

	DO 170 I = 1,nx
	Ps(I,I) = SUM(ZR(I,1:nx)*Z(I,1:nx))
	DO 170 J = 1,I-1
	Ps(I,J) = SUM(ZR(I,1:nx)*Z(J,1:nx))
170	Ps(J,I) = Ps(I,J)

	DEALLOCATE(WORK,RR,WR,WI,Z,T,ZRZ,ZR,Q,WR1,Z1)

7777	RETURN
      END

C For DEEGS - not used
      LOGICAL FUNCTION SELECT(A,B)
      DOUBLE PRECISION A,B
      RETURN
      END


c      SUBROUTINE PROVA(A)
c      DOUBLE PRECISION A
c      DIMENSION A( * )
c      A(4) = 1.D0
c      RETURN
c      END

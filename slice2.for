C -------------------------------------------------------------
C SLICE2 (no missing values) implements the SINGLE-VARIABLE SLICE SAMPLING 
C in Neal (2003), Slice Sampling, Annals of Statistics 31, 705-67 
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
C -------------------------------------------------------------
	SUBROUTINE SLICE2(it,nobs,d,ny,nz,nx,nu,ns,nt,S,yk,
	1                  theta,thetaprior,tipo,pdll,NEVAL,XSIM)
C INPUT
	INTEGER it,nobs,d(2),ny,nz,nx,nu,ns(6),nt,S(nobs,6)
	DOUBLE PRECISION yk(nobs,ny+nz),theta(nt),thetaprior(4)
	CHARACTER*2 tipo
	CHARACTER*1 fittizia
	POINTER (pdll,fittizia)

C OUTPUT
	INTEGER NEVAL
	DOUBLE PRECISION XSIM
C LOCALS
	INTEGER M,J,K,OK
	DOUBLE PRECISION XOLD,XLB,XUB
	DOUBLE PRECISION FXOLD,U,Z,L,R,W,FXL,FXR,FXSIM
	DOUBLE PRECISION genunf,PTHETA2

	NEVAL = 0
	XOLD  = theta(it)
	XLB   = thetaprior(3)
	XUB   = thetaprior(4)

C -------------------------------------------------------
C 1. DRAW Z = ln[f(X0)] - EXP(1) where EXP(1)=-ln(U(0,1))
C    THIS DEFINES THE SLICE S={x: z < ln(f(x))}
C -------------------------------------------------------
	theta(it) = XOLD
	FXOLD=PTHETA2(it,nobs,d,ny,nz,nx,nu,ns,nt,S,yk,
	1              theta,thetaprior,tipo,pdll)
	NEVAL = NEVAL + 1
      U = genunf(0.D0,1.D0)
	Z = FXOLD + DLOG(U)
	
C -------------------------------------------------------------
C 2. FIND I=(L,R) AROUND X0 THAT CONTAINS S AS MUCH AS POSSIBLE 
C    STEPPING-OUT PROCEDURE
C    W = an estimate of the scale of SC    
C    M = Limit on steps (-1 = +INF)
C -------------------------------------------------------------
	M = -1
	W = max((XUB-XLB)/10.0,1.D0) 
C	U = G05CAF(U)
      U = genunf(0.D0,1.D0)
	L = XOLD - W*U
	R = XOLD + W - W*U ! L + W
	IF (M.EQ.-1) THEN
	 DO 100 WHILE(L.GT.XLB)
        theta(it) = L
 	  FXL=PTHETA2(it,nobs,d,ny,nz,nx,nu,ns,nt,S,yk,
	1              theta,thetaprior,tipo,pdll)
	  NEVAL = NEVAL + 1
	  IF (FXL.LE.Z) GOTO 110  
100	   L = L - W
110     DO 200 WHILE(R.LT.XUB)
        theta(it) = R
 	  FXR=PTHETA2(it,nobs,d,ny,nz,nx,nu,ns,nt,S,yk,
	1              theta,thetaprior,tipo,pdll)
	  NEVAL = NEVAL + 1
	  IF (FXR.LE.Z) GOTO 210  
200	 R = R + W
210    CONTINUE
	ELSE
C	 U = G05CAF(U)
       U = genunf(0.D0,1.D0)
	 J = M*U
	 K = M-1-J
	 DO 300 WHILE (J.GT.0)  
	  IF (L.LE.XLB) GOTO 310
	  theta(it) = L
 	  FXL=PTHETA2(it,nobs,d,ny,nz,nx,nu,ns,nt,S,yk,
	1              theta,thetaprior,tipo,pdll)
	  NEVAL = NEVAL + 1  
	  IF (FXL.LE.Z) GOTO 310
	  L = L - W
300     J = J - 1
310    CONTINUE
	  DO 400 WHILE (K.GT.0)  
	  IF (R.GE.XLB) GOTO 410
	  theta(it) = R
 	  FXR=PTHETA2(it,nobs,d,ny,nz,nx,nu,ns,nt,S,yk,
	1              theta,thetaprior,tipo,pdll)
	  NEVAL = NEVAL + 1  
	  IF (FXR.LE.Z) GOTO 410
	  R = R + W
400     K = K - 1
410    CONTINUE
	ENDIF
	IF (L.LT.XLB) L = XLB
	IF (R.GT.XUB) R = XUB
	
C ------------------------------------------------------
C 3. SAMPLING FROM THE SET A = (I INTERSECT S) = (LA,RA)
C ------------------------------------------------------
	OK = 0
	DO 500 WHILE (OK.EQ.0)
C	 U = G05CAF(U)
       U = genunf(0.D0,1.D0)
	 XSIM = L + U*(R - L)
	 theta(it) = XSIM
	 FXSIM=PTHETA2(it,nobs,d,ny,nz,nx,nu,ns,nt,S,yk,
	1               theta,thetaprior,tipo,pdll)
	 NEVAL = NEVAL + 1  
	 IF (FXSIM.GE.Z) OK = 1
	 IF(XSIM.GT.XOLD) THEN
	  R = XSIM
	 ELSE
	  L = XSIM
	 ENDIF 
500	CONTINUE	

	theta(it) = XOLD     

	RETURN
	END
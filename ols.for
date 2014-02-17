C --------------------------------------------------
C OLS reurns the ordinary least square estimates of 
C a linear regression
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C      
C  INPUT: matrix of regressors   X (N x K),
C         vector of observations Y (N x 1)
C         N number of observation
C		K number of regressors
C
C  OUTPUT: BETA    = model parameters 
C          SEB     = standard error of parameters
C          SIGMA   = var covar matrix of parameters
C          RES     = model residuals
C          VA      = variance of residuals
C          IFAULT  = OUTPUT, ERROR INDICATOR:
C                      1 IF N < 1
C                      2 IF X'X IS NOT +VE SEMI-DEFINITE
C                      0 OTHERWISE
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
C ------------------------------------------------
      SUBROUTINE OLS(Y,X,N,K,BETA,SEB,SIGMA,RES,VA,IFAULT)
	INTEGER N,K
	DOUBLE PRECISION Y(N),X(N,K),BETA(K),SEB(K),SIGMA(K,K)
	DOUBLE PRECISION RES(N),VA
	INTEGER I,J,IL,IFAULT,NULLTY
	DOUBLE PRECISION XX(K,K),invXX(K,K),LTRXX(K*(K+1)/2),W(K),
     1 pro(K,N),RMAX

C Compute X'X 	
      DO I = 1,K
	 XX(I,I) = SUM(X(1:N,I)*X(1:N,I))
	 DO J = I+1,K
	  XX(I,J) = SUM(X(1:N,I)*X(1:N,J))
	  XX(J,I) = XX(I,J) 
 	 ENDDO
	ENDDO

C Compute inv(X'X) 	
	DO 10 i=1,K
	DO 10 j=1,i 
10    LTRXX(i*(i+1)/2-i+j)=XX(i,j)

      CALL SYMINV(LTRXX,K,LTRXX,W,NULLTY,IFAULT,RMAX)
      DO i = 1,K
	 invXX(i,i) = LTRXX(i*(i+1)/2)
	 DO j = 1,i-1
	  invXX(i,j) = LTRXX(i*(i+1)/2-i+j)
	  invXX(j,i) = invXX(i,j)
	 ENDDO
	ENDDO	

C beta=inv(x'*x)*x'*y;
      DO I = 1,K
	 DO J = 1,N
	  pro(I,J) = SUM(invXX(I,:)*X(J,:))
 	 ENDDO
	ENDDO
      DO I = 1,K
	 beta(I) = SUM(pro(I,:)*Y(:))
 	ENDDO
	 	   	  	   
C Residuals: res=y-x*beta
	DO I = 1,N
	 RES(I) = Y(I)-SUM(X(I,:)*beta(:))
 	ENDDO
	   
C Variance of residuals: va=(res'res)/n
	VA = SUM(RES(:)*RES(:))/DFLOAT(N)

C Var-Covar matri:x  sigma=va*inv(x'*x);
      SIGMA = VA*invXX 

C Coefficient Standard errors 	   
	DO 20 i=1,K
20	SEB(i) = dsqrt(SIGMA(i,i))

      RETURN
	END
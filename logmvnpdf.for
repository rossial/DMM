C ------------------------------------------------------------------------ 
C LOGMVNPDF returns the Multivariate Normal pdf with parameters mu, SIG, 
C evaluated at x.
C mu (px1),SIG (pxp), x (px1)
C Bauwens et al. (1999): "Bayesian Inference in Dynamic Econometric 
C models", Oxford University Press, page 298
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
C -------------------------------------------------------------------------- 
	DOUBLE PRECISION FUNCTION logmvnpdf(X,mu,SIG,p)

C INPUT
	INTEGER p
      DOUBLE PRECISION X(p),mu(p),SIG(p,p)
C LOCALS
	INTEGER I,J,IFAIL
	DOUBLE PRECISION,ALLOCATABLE:: ISIG(:,:),COM(:,:),WKSPCE(:)
      DOUBLE PRECISION PI,H,DET,KER
	DATA PI/3.141592653589793D0/,H/-.5D0/	

      ALLOCATE(ISIG(p,p),COM(p+1,p),WKSPCE(p))
      
C INVERT SIG 
      COM(1:p,1:p) = SIG(1:p,1:p)
	IFAIL = -1
	CALL F01ADF(p,COM(1:p+1,1:p),p+1,IFAIL) 
	DO 10 I=1,p
	DO 10 J=1,I
	ISIG(I,J) = COM(I+1,J)
10	ISIG(J,I) = ISIG(I,J)    

C DETERMINANT of SIG	
	COM(1:p,1:p) = SIG(1:p,1:p)
	IFAIL = -1
      CALL F03ABF(COM(1:p,1:p),p,p,DET,WKSPCE,IFAIL) 

C QUADRATIC FORM (log)
	KER = 0.D0
	DO 20 I=1,p
	DO 20 J=1,p
20	KER = KER + (X(I)-mu(I))*ISIG(I,J)*(X(J)-mu(J))

C LOG PDF
	logmvnpdf = H*p*dlog(2.*PI)+H*dlog(DET)+H*KER
	
      DEALLOCATE(ISIG,COM,WKSPCE)
	
      RETURN
	END
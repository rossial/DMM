C ------------------------------------------------------------------------
C MVNPDF returns the multivariate Normal pdf with parameters mu, SIG,
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
C --------------------------------------------------------------------------
	DOUBLE PRECISION FUNCTION mvnpdf(X,mu,SIG,p)

C INPUT
	INTEGER p
	DOUBLE PRECISION X(p),mu(p),SIG(p,p)
C LOCALS
	INTEGER I,J,IFAIL
	DOUBLE PRECISION,ALLOCATABLE:: ISIG(:,:),COM(:,:)
	DOUBLE PRECISION PI,H,DET,KER
	DATA PI/3.141592653589793D0/,H/-.5D0/
C EXTERNAL SUBROUTINES
      EXTERNAL DPOTRF,DPOTRI

      ALLOCATE(ISIG(p,p),COM(p+1,p))
      COM(1:p,1:p) = SIG(1:p,1:p)
	IFAIL = -1
C	CALL F01ADF(p,COM(1:p+1,1:p),p+1,IFAIL)
      CALL DPOTRF('L',p,COM(1:p,1:p),p,IFAIL) ! COM = L*L'
      DET = 1.D0  ! det(L)
      DO I =1,p
       DET = DET*COM(I,I)
      ENDDO
      CALL DPOTRI('L',p,COM(1:p,1:p),p,IFAIL) ! COM = VV^-1

	DO 10 I=1,p
      ISIG(I,I) = COM(I,I)
	DO 10 J=1,I-1
	ISIG(I,J) = COM(I,J)
10	ISIG(J,I) = ISIG(I,J)

C DETERMINANT of SIG
C	COM(1:p,1:p) = SIG(1:p,1:p)
C	IFAIL = -1
C     CALL F03ABF(COM(1:p,1:p),p,p,DET,WKSPCE,IFAIL)

C QUADRATIC FORM (log)
	KER = 0.D0
	DO 20 I=1,p
	DO 20 J=1,p
20	KER = KER + (X(I)-mu(I))*ISIG(I,J)*(X(J)-mu(J))

C PDF
	mvnpdf = (2.*PI)**(H*p)*DET**(-1.D0)*DEXP(H*KER)

      DEALLOCATE(ISIG,COM)
	RETURN
	END

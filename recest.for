C -----------------------------------------------
C RECEST performs recursive estimation of the 
C MEAN, SD, and COVARIANCE MATRIX
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
C -----------------------------------------------
	SUBROUTINE RECEST(NX,N,X,MED0,SIG0,SK0,MED,SIG,SK)
C INPUT
	INTEGER NX,N
	DOUBLE PRECISION X(NX),MED0(NX),SK0(NX),SIG0(NX,NX)
C OUTPUT
	DOUBLE PRECISION MED(NX),SK(NX),SIG(NX,NX)
C LOCALS
	INTEGER I,J
	DOUBLE PRECISION SSX,SSXY,S3X

C MEAN
	MED = (N*MED0+X)/DFLOAT(N+1)
C VAR	and SKEW 
	DO I = 1,NX
	 SSX = SIG0(I,I) + MED0(I)**2
	 S3X = SK0(I)*SIG0(I,I)**1.5D0+3.D0*MED0(I)*SSX-2.D0*MED0(I)**3 
	 SSX = (N*SSX+X(I)**2)/DFLOAT(N+1)
	 SIG(I,I) = SSX-MED(I)**2
	 S3X = (S3X*N+X(I)**3)/DFLOAT(N+1)
	 DO J = 1,I-1
	  SSXY = SIG0(I,J) + MED0(I)*MED0(J)
	  SIG(I,J) = (N*SSXY+X(I)*X(J))/DFLOAT(N+1)-MED(I)*MED(J)
	  SIG(J,I) = SIG(I,J)
       ENDDO
	 SK(I) = SIG(I,I)**(-1.5D0)*(S3X-3.D0*MED(I)*SSX+2.D0*MED(I)**3)
	ENDDO     
	RETURN
	END

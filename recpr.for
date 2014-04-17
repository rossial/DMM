C ------------------------------------------------------------------------
C RECPR compute transition and marginal probabilities for
C the adaptive MH block-sampler (see Fiorentini,Planas
C and Rossi, Statistics and Computing 2014)
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
C ------------------------------------------------------------------------
	SUBROUTINE RECPR(N,NS,nobs,S,SW,PM,PTR)
#ifdef __GFORTRAN__
      USE dynare
#endif
! INPUT
	INTEGER N,NS,nobs
	INTEGER S(nobs),SW(2*nobs)
! INPUT/OUTPUT
	DOUBLE PRECISION PM(nobs,NS),PTR(nobs,NS,NS)
! LOCALS
	INTEGER I,J,IT
#ifdef __GFORTRAN__
      INTEGER tmp0, tmp1
#endif
C Use previous run
      SW(1:nobs)        = SW(nobs+1:2*nobs)
	SW(nobs+1:2*nobs) = S(1:nobs)

C Transition probs Pr[S(t)|S(t-1),Y]
	DO 10 IT = 2,nobs
	DO 10 J  = 1,NS
	DO 9  I  = 1,NS-1
#ifdef __GFORTRAN__
       tmp0 = SW(IT).EQ.I
       tmp1 = SW(IT-1).EQ.J
9	PTR(IT,I,J) = (N*PTR(IT,I,J)*PM(IT-1,J)
     #            + tmp0*tmp1)
     #            / (N*PM(IT-1,J)+ABS(tmp1))
#else
9	PTR(IT,I,J) = (N*PTR(IT,I,J)*PM(IT-1,J)
     #            + (SW(IT).EQ.I)*(SW(IT-1).EQ.J))
     #            / (N*PM(IT-1,J)+ABS(SW(IT-1).EQ.J))
#endif
10	PTR(IT,NS,J) = 1.D0-SUM(PTR(IT,1:NS-1,J))

	DO 30 IT = 1,nobs
	DO 20 I  = 1,NS-1
#ifdef __GFORTRAN__
       tmp0 = SW(IT).EQ.I
20	PM(IT,I) = (N*PM(IT,I)
     #         + ABS(tmp0))/DFLOAT(N+1)
#else
20	PM(IT,I) = (N*PM(IT,I)
     #         + ABS(SW(IT).EQ.I))/DFLOAT(N+1)
#endif
30	PM(IT,NS) = 1.D0 - SUM(PM(IT,1:NS-1))

	RETURN
	END

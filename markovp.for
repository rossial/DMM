C ---------------------------------------------------------------------------
C MARKOVP computes the conditional Log-Probability of a Block of S variables:
C log{P[S(t+1),...,S(t+h)|S(1),...,S(t),S(t+h+1),...,S(T)]}
C HFIX = block length
C NH   = 1 first block; 
C NH   > 1 and NH < nobs/HFIX middle block
C NH   = nob/HFIX last block
C S    = {S(t),S(t+1),...,S(t+h),S(t+h+1)} 
C S(t) takes values in {1,2,...,N}
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
C ---------------------------------------------------------------------------
	DOUBLE PRECISION FUNCTION MARKOVP(PMAT,PE,N,HFIX,NH,nobs,S)

C INPUT 
	INTEGER N,HFIX,NH,S(HFIX+2),nobs
	DOUBLE PRECISION PMAT(N,N,HFIX+1),PE(N)

C LOCALS
	INTEGER IT

	IF (NH.EQ.1) THEN

	 MARKOVP = DLOG(PMAT(S(HFIX+2),S(2),HFIX))+DLOG(PE(S(2))) 
     +	     - DLOG(PE(S(HFIX+2)))
	 DO 10 IT=2,HFIX
10	 MARKOVP = MARKOVP + DLOG(PMAT(S(IT+1),S(IT),1)) 
     +         + DLOG(PMAT(S(HFIX+2),S(IT+1),HFIX+1-IT))
     +         - DLOG(PMAT(S(HFIX+2),S(IT),HFIX+2-IT))

	ELSEIF ((NH.GT.1).AND.(NH.LT.nobs/HFIX)) THEN

	 MARKOVP = 0.D0
	 DO 20 IT=1,HFIX
20	 MARKOVP = MARKOVP + DLOG(PMAT(S(IT+1),S(IT),1))
     +         + DLOG(PMAT(S(HFIX+2),S(IT+1),HFIX+1-IT))
     +         - DLOG(PMAT(S(HFIX+2),S(IT),HFIX+2-IT))

	ELSE

	 MARKOVP = 0.D0
	 DO 30 IT=1,HFIX
30	 MARKOVP = MARKOVP + DLOG(PMAT(S(IT+1),S(IT),1))
	
	ENDIF
	RETURN
	END
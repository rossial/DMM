C ------------------------------------------------------------
C INT2SEQ converts integers into the S-sequence
C int   = integer to be converted
C nv    = # of S variables
C INFOS = info for S-var
C SEQ   = S-sequence
C IS    = map to c, H, G, etc
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
C -----------------------------------------------------------------------
	SUBROUTINE INT2SEQ(int,nv,INFOS,SEQ,IS)
C INPUT	
	INTEGER int,nv,INFOS(9,6)
C OUTPUT
	INTEGER SEQ(nv),IS(6)
C LOCALS
	INTEGER M,i,j,k,ns(nv)

	IS(:)  = 1
	SEQ(:) = 1
	ns(1:nv) = INFOS(8,1:nv)
	j = int
	i = nv
	DO WHILE (i.GT.1)
        M = PRODUCT(ns(nv-i+2:nv))
        DO 10 k = 1,ns(nv-i+1)
         IF (j.LE.k*M) THEN
           SEQ(nv-i+1) = k
           GOTO 11   
         ENDIF
10      CONTINUE
11      j = j-(k-1)*M
        i = i - 1
	ENDDO
	SEQ(nv) = j
	DO 20 i = 1,nv
	 DO 20 j = 1,INFOS(1,i)
20	 IS(INFOS(j+1,i)) = SEQ(i)
	RETURN
	END

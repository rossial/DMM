C -------------------------------------------------------------
C FINFINPUT finds input values in the inputfile
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C      
C INPUT:
C   STR1 
C   STR2
C
C OPUTPUT:
C   NUM the value to be recovered 
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
	SUBROUTINE FINDINPUT(STR1,N1,STR2,N2,NUM)
C INPUT
	INTEGER N1,N2
	CHARACTER STR1*N1,STR2*N2	
C OUTPUT
	DOUBLE PRECISION NUM
C LOCALS
      INTEGER IPOS,IPOS2
	CHARACTER*1 :: OSTR
	CHARACTER*300 ERRMSG

	NUM   = 0.D0
	OSTR  = ' '
	ERRMSG = 'INPUT ERROR: '// STR2 //' is not set' 
	IPOS = INDEX(TRIM(STR1),STR2)
	IF (IPOS.EQ.0) THEN 
	 WRITE(*,*) ERRMSG
	 PAUSE
	 STOP
      ELSE	 
	 IPOS = IPOS+2
	 DO WHILE ((OSTR.NE.'=').AND.(OSTR.NE.'$'))
	  READ(STR1(IPOS:IPOS),*) OSTR
	  IPOS = IPOS+1
	 END DO
	 IPOS2 = IPOS
	 DO WHILE ((OSTR.NE.';').AND.(OSTR.NE.'$'))
	  READ(STR1(IPOS2:IPOS2),*) OSTR
	  IPOS2 = IPOS2 + 1
	 END DO
	 IF (IPOS2-1.GT.IPOS) THEN
	  READ(STR1(IPOS:IPOS2-2),'(F20.10)') NUM
	  READ(STR1(IPOS:IPOS2-2),*) NUM	  
	 ELSE	  
	  WRITE(*,*) ERRMSG
	  PAUSE
	  STOP
	 ENDIF	   		 
	ENDIF 

      RETURN
	END
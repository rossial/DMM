C ------------------------------------------------------------
C INITRAND initializes the random number generator:
C INITIAL = 0 FIXED SEED, otherwise SEED set according to clock 
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
C ------------------------------------------------------------
      SUBROUTINE INITRAND(INITIAL,DATE_ITIME)
C Input
      INTEGER INITIAL,DATE_ITIME(8)

      IF(INITIAL.NE.0) THEN
       CALL setall(DATE_ITIME(7),DATE_ITIME(8))      
      ELSE
       CALL setall(12345,54321)
      ENDIF
      RETURN
      END

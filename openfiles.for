C --------------------------------------------------------------------------
C OPENFILES opens files for writing the DMM output       
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
C --------------------------------------------------------------------------
      SUBROUTINE OPENFILES(ESTIMATION,SEED,NV,NF,NMIS,SIMULATION,
	1                     MARGLIK,PATH,NMLNAME)
      INTEGER SEED,nv,NF,NMIS
	CHARACTER*1 SIMULATION,MARGLIK
	CHARACTER*2 ESTIMATION
	CHARACTER*3 CSEED
	CHARACTER*200 NMLNAME,PATH,FILEOUT 

	IF ((estimation.EQ.'ML').OR.(estimation.EQ.'ml').OR.
     &    (estimation.EQ.'Ml').OR.(estimation.EQ.'mL')) THEN

	 FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//'.PAR'
	 OPEN(9, FILE = FILEOUT, ACCESS='SEQUENTIAL') 
      
	 FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//'.UNB'
	 OPEN(10,FILE = FILEOUT, ACCESS='SEQUENTIAL')

       IF (nv.GT.0) THEN
	  FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//'.DIS'
 	  OPEN(11,FILE = FILEOUT, ACCESS='SEQUENTIAL')
       ENDIF 
       
     	 FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//'.INN'
	 OPEN(12, FILE = FILEOUT, ACCESS='SEQUENTIAL')
	
	ELSE

	 IF (SEED.LE.9) THEN 
	  WRITE(CSEED,'(I1)') SEED
	 ELSEIF ((SEED.GE.10).AND.(SEED.LE.99)) THEN 
	  WRITE(CSEED,'(I2)') SEED
	 ELSEIF ((SEED.GE.100).AND.(SEED.LE.999)) THEN 
	  WRITE(CSEED,'(I3)') SEED
	 ENDIF
	
	 IF (nv.GT.0) THEN
	  FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//TRIM(CSEED)//'.DIS'
 	  OPEN(11,FILE = FILEOUT, ACCESS='SEQUENTIAL')
	 ENDIF 
		
	 FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//TRIM(CSEED)//'.PAR'
	 OPEN(9, FILE = FILEOUT, ACCESS='SEQUENTIAL') 
      
	 FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//TRIM(CSEED)//'.UNB'
	 OPEN(10,FILE = FILEOUT, ACCESS='SEQUENTIAL')
	
	 IF ((SIMULATION.EQ.'N').OR.(SIMULATION.EQ.'n')) THEN
	  FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//TRIM(CSEED)//'.INN'
	  OPEN(12, FILE = FILEOUT, ACCESS='SEQUENTIAL')
	
	  IF (nf.GT.0) THEN
         FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//TRIM(CSEED)//'.FST'
	   OPEN(13,FILE = FILEOUT, ACCESS='SEQUENTIAL')
	  ENDIF
	
	  IF (nmis.GT.0) THEN
         FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//TRIM(CSEED)//'.MIS'
	   OPEN(14,FILE = FILEOUT, ACCESS='SEQUENTIAL')
	  ENDIF
	 
	  IF ((MARGLIK.EQ.'Y').OR.(MARGLIK.EQ.'y')) THEN
         FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//TRIM(CSEED)//'.ML'
	   OPEN(15,FILE = FILEOUT, ACCESS='SEQUENTIAL')
	  ENDIF
	
	 ELSE
	
	  FILEOUT = TRIM(PATH)//TRIM(NMLNAME)//TRIM(CSEED)//'.DAT'
	  OPEN(15,FILE = FILEOUT, ACCESS='SEQUENTIAL')
	
	 ENDIF

	ENDIF

      RETURN
	END

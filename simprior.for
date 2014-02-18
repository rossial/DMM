C --------------------------------------------------------------------
C SIMPRIOR SIMULATES theta FROM THE PRIOR DISTRIBUTION 
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
C --------------------------------------------------------------------	
	SUBROUTINE SIMPRIOR(estimation,nt,thetaprior,tipo,ntf,INDT,theta)
C INPUT
	INTEGER nt
	DOUBLE PRECISION thetaprior(nt,4)
	CHARACTER*2 tipo(nt),estimation
C OUTPUT
	INTEGER ntf,INDT(nt+2)
	DOUBLE PRECISION theta(nt)
C LOCALS
	INTEGER I,IFAIL
	DOUBLE PRECISION LB,UB,PLB,PUB
	DOUBLE PRECISION S15ABF,TNORMI,gengam,genbet

	INDT(:) = 0
	ntf     = 0
	DO I = 1,nt
	 IF (thetaprior(I,3).LT.thetaprior(I,4)) THEN
	  ntf       = ntf + 1
	  INDT(ntf) = I
	  IF (estimation.EQ.'BA') THEN
         IF (tipo(I).EQ.'IG') THEN
	    IFAIL = -1
C   	    CALL G05FFF(thetaprior(I,2)/2.D0,2.D0/thetaprior(I,1),1,
C    1                theta(I),IFAIL)
          theta(I) = 1.D0
     #             / gengam(thetaprior(I,1)/2.D0,thetaprior(I,2)/2.D0) 
	    IF( (theta(I).LT.thetaprior(I,3)).OR.
     #        (theta(I).GT.thetaprior(I,4)) ) THEN
	      theta(I) = (thetaprior(I,4)+thetaprior(I,3))/2.D0 
	    ENDIF  
         ELSEIF (tipo(I).EQ.'NT') THEN  
	    LB = (thetaprior(I,3)-thetaprior(I,1))/DSQRT(thetaprior(I,2))
	    UB = (thetaprior(I,4)-thetaprior(I,1))/DSQRT(thetaprior(I,2))
	    PLB = S15ABF(LB,IFAIL)
	    PUB = S15ABF(UB,IFAIL)
	    theta(I) = TNORMI(PLB,PUB)
	    theta(I) = thetaprior(I,1)+theta(I)*DSQRT(thetaprior(I,2))	 
	   ELSEIF (tipo(I).EQ.'BE') THEN 
	    IFAIL = -1
C	    CALL G05FEF(thetaprior(I,1),thetaprior(I,2),1,theta(I),IFAIL)
          theta(I) = genbet(thetaprior(I,1),thetaprior(I,2))
	    theta(I) = theta(I)*(thetaprior(I,4)-thetaprior(I,3)) 
     +             + thetaprior(I,3)       
         ENDIF         
        ELSE
         theta(I) = thetaprior(I,1)
        ENDIF    
	 ELSE
	  theta(I) = thetaprior(I,3)  
	 ENDIF
	ENDDO
	INDT(nt+2) = ntf  ! # OF FREE PARS	
	RETURN
	END
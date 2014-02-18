C -------------------------------------------------------------------
C PRIOR COMPUTES THE LOG-VALUE PRIOR pdf EVALUATED AT THETA 
C ACCORDING to 'TIPO'.
C Developed by A.Rossi, C.Planas and G.Fiorentini     
C
C TIPO: 'BE' = Beta over (a,b) 
C       'IG' = Inverted Gamma, parameterization as in Bauwens et al. 
C       'NT' = Truncated Normal(mean,variance)
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
C -------------------------------------------------------------------
	DOUBLE PRECISION FUNCTION PRIOR(theta,thetaprior,tipo)  	
C INPUT
	DOUBLE PRECISION theta,thetaprior(4)
	CHARACTER*2 tipo

C LOCALS
      INTEGER IFAIL
	DOUBLE PRECISION XMIN,XMAX,AUX,EPS,P,Q
	DOUBLE PRECISION S14ABF,S15ABF,PI
	EXTERNAL S14BAF
	DATA PI/3.141592653589793D0/,EPS/1.D-14/	
	
	IF (tipo.EQ.'IG') THEN
        IFAIL= -1
        PRIOR = -S14ABF(.5D0*thetaprior(2),IFAIL)
     +	    - .5D0*thetaprior(2)*DLOG(2.D0/thetaprior(1))
     +        - (.5D0*thetaprior(2)+1.D0)*DLOG(theta)
     +        - .5D0*thetaprior(1)/theta 
	  XMAX = thetaprior(1)/(2.D0*thetaprior(4))
        CALL S14BAF(thetaprior(2)/2.D0,XMAX,EPS,P,Q,IFAIL) 
        PRIOR = PRIOR -DLOG(Q)
      ELSEIF (tipo.EQ.'NT') THEN  
	  XMIN  = (thetaprior(3)-thetaprior(1))/DSQRT(thetaprior(2))
	  XMAX  = (thetaprior(4)-thetaprior(1))/DSQRT(thetaprior(2))
	  IFAIL = -1
	  AUX   = S15ABF(XMAX,IFAIL)-S15ABF(XMIN,IFAIL)
	  PRIOR = -.5D0*DLOG(thetaprior(2)) 
     +	    -  .5D0*(theta-thetaprior(1))**2/thetaprior(2)
     +        - DLOG(AUX) -.5D0*DLOG(2.D0*PI)
	ELSEIF (tipo.EQ.'BE') THEN 
        AUX = S14ABF(thetaprior(1) + thetaprior(2),IFAIL)
     +	  - S14ABF(thetaprior(1),IFAIL)
     +      - S14ABF(thetaprior(2),IFAIL)
        XMIN = (theta-thetaprior(3))/(thetaprior(4)-thetaprior(3))                                 
        PRIOR = AUX + (thetaprior(1)-1.D0)*DLOG(XMIN) 
     +        + (thetaprior(2)-1.D0)*DLOG(1.D0-XMIN)
     +        - DLOG(thetaprior(4)-thetaprior(3))
      ENDIF
      RETURN
      END
C --------------------------------------------------------------------------------      
C INVNORMCDF Computes the inverse of the N(0,1) cdf according to 
C the algorithm shown in Wichura, M.J. (1988).
C Algorithm AS 241: The Percentage Points of the Normal Distribution.
C Applied Statistics, 37, 477-484.
C
C Copyright (C) 2002 Przemyslaw Sliwa and Jason H. Stover.
C
C Recoded in Fortran by A. Rossi
C Copyright (C) 2014 European Commission
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
C --------------------------------------------------------------------------------      
      DOUBLE PRECISION FUNCTION INVNORMCDF(P)
C INPUT
      DOUBLE PRECISION P
C LOCALS
      DOUBLE PRECISION r,x,pp,dP
C EXTERNAL FUNCIONS      
      DOUBLE PRECISION small,intermediate,tail 
      
      IF (P.GE.1.D0) THEN
       INVNORMCDF = 1.D100
       RETURN
      ENDIF    
      
      IF (P.LE.0.D0) THEN
       INVNORMCDF = -1.D100
       RETURN
      ENDIF    
      
      dP = P - 0.5D0
      IF (DABS(dP).LE.0.425D0) THEN
       INVNORMCDF = small(dP)   
       RETURN
      ENDIF
      
      IF (P.LT.0.5D0) THEN
       pp = P
      ELSE
       pp = 1.D0-P   
      ENDIF 
      
      r = dsqrt(-dlog(pp))
      IF (r.LE.5.D0) THEN
       x = intermediate(r)       
      ELSE
       x = tail(r)
      ENDIF

      IF (P.LT.0.5D0) THEN 
       INVNORMCDF = -x
      ELSE
       INVNORMCDF = x
      ENDIF            
      
      RETURN
      END        

C *************************************************************      
      double precision function small(q)
      double precision q,r,a(8),b(8)
      double precision rat_eval
      data a/3.387132872796366608,133.14166789178437745,
     * 1971.5909503065514427, 13731.693765509461125,
     * 45921.953931549871457, 67265.770927008700853,
     * 33430.575583588128105, 2509.0809287301226727/
      data b/ 1.D0,42.313330701600911252,
     * 687.1870074920579083, 5394.1960214247511077,
     * 21213.794301586595867, 39307.89580009271061,
     * 28729.085735721942674, 5226.495278852854561/
  
      r = .180625D0 - q * q
      small = q * rat_eval(a, 8, b, 8, r)
  
      return
      end     
C *************************************************************            
      double precision function intermediate(r)
      double precision r,a(8),b(8)
      double precision rat_eval
      data a/1.42343711074968357734, 4.6303378461565452959,
     * 5.7694972214606914055, 3.64784832476320460504,
     * 1.27045825245236838258, 0.24178072517745061177,
     * 0.0227238449892691845833, 7.7454501427834140764d-4/
      data b/ 1.D0, 2.05319162663775882187,
     * 1.6763848301838038494, 0.68976733498510000455,
     * 0.14810397642748007459, 0.0151986665636164571966,
     * 5.475938084995344946d-4, 1.05075007164441684324d-9/
  
      intermediate = rat_eval(a, 8, b, 8, r-1.6D0)
  
      return
      end     
C *************************************************************                  
      double precision function tail(r)
      double precision r,a(8),b(8)
      double precision rat_eval
      data a/6.6579046435011037772, 5.4637849111641143699,
     * 1.7848265399172913358, 0.29656057182850489123,
     * 0.026532189526576123093, 0.0012426609473880784386,
     * 2.71155556874348757815d-5, 2.01033439929228813265d-7/
      data b/ 1.D0, 0.59983220655588793769,
     * 0.13692988092273580531, 0.0148753612908506148525,
     * 7.868691311456132591d-4, 1.8463183175100546818d-5,
     * 1.4215117583164458887d-7, 2.04426310338993978564d-15/
  
      tail = rat_eval(a, 8, b, 8, (r - 5.D0))
      
      return
      end     
C *************************************************************                              
      double precision function rat_eval(a,na,b,nb,x)
      integer na,nb,i
      double precision a(na),b(nb),x,u,v
      
      u = a(na)      
      do i=na,2,-1
       u = x*u+a(i-1)    
      enddo
      
      v = b(nb)
      do i=nb,2,-1
       v = x*v+b(i-1)    
      enddo
      rat_eval = u/v
            
      return
      end  
C *************************************************************                  
      
      
      
      
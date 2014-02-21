C -------------------------------------------------------------
C PRIORDIR COMPUTES THE LOG PRIOR DISTRIBUTION EVALUATED AT psi
C FOR THE DIRICHLET pdf 
C f(psi(1),..,psi(N-1);a1,...,aN)=1/B(a)*prod_i=1^N psi(i)^(ai-1)
C B(a) = prod_i=1^N G(ai)/G(a0), a0 = a1+...+aN 
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
C --------------------------------------------------------------
	DOUBLE PRECISION FUNCTION PRIORDIR(psi,psiprior,N)  	
C INPUT
	INTEGER N
	DOUBLE PRECISION psi(N-1),psiprior(N)

C LOCALS
      INTEGER I
      
C EXTERNAL FUNCTIONS      
	DOUBLE PRECISION GAMMLN
	
	PRIORDIR = gammln(SUM(psiprior(1:N)))
      DO I = 1,N-1
	 PRIORDIR = PRIORDIR - gammln(psiprior(I))
     #          + (psiprior(I)-1.D0)*DLOG(psi(I))
	ENDDO
	PRIORDIR = PRIORDIR - gammln(psiprior(N))
     #         + (psiprior(N)-1.D0)*DLOG(1.D0-SUM(psi(1:N-1)))

      RETURN
      END
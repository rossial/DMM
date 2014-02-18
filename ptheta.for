C ----------------------------------------------------------------
C PTHETA COMPUTES THE LOG POSTERIOR 
C P(theta1(it)|theta1(~it),S,Y)
C Developed by A.Rossi, C.Planas and G.Fiorentini             
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
C ----------------------------------------------------------------
	DOUBLE PRECISION FUNCTION PTHETA(it,nobs,d,ny,nz,nx,nu,ns,nt,
	1 S,yk,IYK,theta,thetaprior,tipo,pdll)

	USE dfwin
	INTERFACE
	 SUBROUTINE DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R) 
	 INTEGER ny,nz,nx,nu,ns(6),nt
	 DOUBLE PRECISION theta(nt)
	 DOUBLE PRECISION c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6))
	 END SUBROUTINE
	END INTERFACE
	CHARACTER*1 fittizia
	POINTER (pdll,fittizia)  ! ASSOCIATE  pointer pdll alla DLL ad una varibile fittizia
	POINTER (pdesign,DESIGN) ! IMPORTANT associo il puntatore pdesign alla Interface definita
	
C INPUT
	INTEGER it,nobs,d(2),ny,nz,nx,nu,ns(6),nt,S(nobs,6),
	1 IYK(nobs,ny+1)
	DOUBLE PRECISION yk(nobs,ny+nz),theta(nt),thetaprior(4)
	CHARACTER*2 tipo 

C LOCALS
      DOUBLE PRECISION,ALLOCATABLE:: c(:,:,:),H(:,:,:),G(:,:,:),a(:,:),
     1 F(:,:,:),R(:,:,:),LIKE(:),XT(:,:),PT(:,:,:),Xdd(:,:),Pdd(:,:,:)	
      DOUBLE PRECISION PRIOR
      
      ALLOCATE(c(ny,max(nz,1),ns(1)),H(ny,nx,ns(2)),G(ny,nu,ns(3)),
     1 a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6)),LIKE(nobs),
	2 XT(0:nobs,nx),PT(0:nobs,nx,nx),Xdd(max(d(1),1),nx),
     3 Pdd(max(d(1),1),nx,nx))	      
	
C computes the log-posterior
	pdesign = getprocaddress(pdll, "design_"C)
	CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
	CALL IKF(d,ny,nz,nx,nu,ns,S(1:max(d(1),1),1:6),
	1         yk(1:max(d(1),1),1:ny+nz),IYK(1:max(d(1),1),1:ny+1),
     2         c,H,G,a,F,R,Xdd,Pdd,LIKE(1:max(d(1),1)))
	XT(d(1),1:nx)      = Xdd(max(d(1),1),1:nx)   
	PT(d(1),1:nx,1:nx) = Pdd(max(d(1),1),1:nx,1:nx)
	CALL KF(nobs,d,ny,nz,nx,nu,ns,S,yk,IYK,c,H,G,a,F,R,XT,PT,LIKE)
	PTHETA = SUM(LIKE) + PRIOR(theta(it),thetaprior,tipo) 
      
      DEALLOCATE(c,H,G,a,F,R,LIKE,XT,PT,Xdd,Pdd)
	RETURN
	END
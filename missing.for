C --------------------------------------------------------------------
C MISSING Simulates missing observations
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
C --------------------------------------------------------------------
	SUBROUTINE MISSING(yk,ny,nz,nx,nu,ns,nt,nmis,theta,S,STATE,pdll,ykmis)
#if !defined(__GFORTRAN__)
	USE dfwin
#endif
	INTERFACE
	 SUBROUTINE DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)
	 INTEGER ny,nz,nx,nu,ns(6),nt
	 DOUBLE PRECISION theta(nt)
	 DOUBLE PRECISION c(ny,max(1,nz),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)),R(nx,nu,ns(6))
	 END SUBROUTINE
	END INTERFACE
	CHARACTER*1 fittizia
	POINTER (pdll,fittizia)
	POINTER (pdesign,DESIGN)

C INPUT
	INTEGER ny,nz,nx,nu,nt,ns(6),nmis
	DOUBLE PRECISION yk(ny+nz),theta(nt),STATE(nx)

C OUTPUT
	DOUBLE PRECISION ykmis(nmis)
C LOCALS
      INTEGER S(6),I,J,K,IFAIL,NIYK(ny)
	DOUBLE PRECISION U(nu)
      DOUBLE PRECISION gennor
	DOUBLE PRECISION,ALLOCATABLE::R(:,:,:),c(:,:,:),H(:,:,:),
	1 G(:,:,:),a(:,:),F(:,:,:)

	ALLOCATE(R(nx,nu,ns(6)),c(ny,max(nz,1),ns(1)),H(ny,nx,ns(2)),
	1 G(ny,nu,ns(3)),a(nx,ns(4)),F(nx,nx,ns(5)))
#ifdef __GFORTRAN__
      pdesign = getprocaddress(pdll, "design_")
#else
	pdesign = getprocaddress(pdll, "design_"C)
#endif
	CALL DESIGN(ny,nz,nx,nu,ns,nt,theta,c,H,G,a,F,R)

C NIYK = not(IYK)
	K = 0
      DO J = 1,ny
         IF(yk(J).EQ.-99999.D0) THEN
            K = K+1
            NIYK(K) = J
         ENDIF
      ENDDO

C SAMPLING U
	IFAIL = -1
	U(1:nu) = 0.D0
	DO 20 I = 1,nu
c	CALL G05EAF(U(I),1,1.D0,1,1.D-14,WORKU,3,IFAIL)
c20   CALL G05EZF(U(I),1,WORKU,3,IFAIL)
20    U(I) = gennor(0.D0,1.D0)


C DRAW yk ~ f(yk|x,S,zk,theta)
	DO 30 I = 1,nmis
30    ykmis(I) = SUM(c(NIYK(I),1:nz,S(1))*yk(ny+1:ny+nz))
     +         + SUM(H(NIYK(I),1:nx,S(2))*STATE(1:nx))
     +         + SUM(G(NIYK(I),1:nu,S(3))*U(1:nu))

      DEALLOCATE (R,c,H,G,a,F)

	RETURN
	END

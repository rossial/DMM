C --------------------------------------------------------------------
C CONFUN computes non linear constrains for E04UCF
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
C -------------------------------------------------------------
	SUBROUTINE CONFUN(MODE,NCNLNS,NPAR,LDCJ,NEEDC,CHI,CC,CJAC,
     ,	              NNN,IU,U)
C Input
	INTEGER MODE,NPAR,NCNLNS,LDCJ,NEEDC(NCNLNS),NNN
	INTEGER*8 IU(72)
	DOUBLE PRECISION U(IU(1)*(2*IU(4)+IU(5)+7)+3*IU(8)+2),CHI(NPAR)
C Output
	DOUBLE PRECISION CC(LDCJ),CJAC(LDCJ,NPAR)

	CC(:)     = 0.D0
	CJAC(:,:) = 0.D0

	RETURN
	END
